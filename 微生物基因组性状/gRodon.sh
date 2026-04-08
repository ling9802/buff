# 下面的先放在服务器中跑，生成需要的fna文件，可能也需要相对丰度的文件来计算加权（可选，建议），用sh脚本直接bash运行即可
#!/bin/bash

# ==============================
# 参数设置
# ==============================
BASE_DIR=/home/Yangly/JunLing/WG_metagenomics
READS_DIR=$BASE_DIR/02_cleandata
OUTDIR=$BASE_DIR/gRodon_results_test
THREADS=60

# ==============================
# 跳过样本（当前为空）
# ==============================
SKIP_SAMPLES=()

# ==============================
# 循环处理所有 contig
# ==============================
for CONTIG in $BASE_DIR/*.contig.fa
do
    SAMPLE=$(basename $CONTIG .contig.fa)

    # ==============================
    # 跳过逻辑（严格匹配）
    # ==============================
    for skip in "${SKIP_SAMPLES[@]}"; do
        if [ "$SAMPLE" == "$skip" ]; then
            echo "⚠ Skipping $SAMPLE"
            continue 2
        fi
    done

    R1=${READS_DIR}/${SAMPLE}.R1.clean.fastq
    R2=${READS_DIR}/${SAMPLE}.R2.clean.fastq
    SAMPLE_OUT=${OUTDIR}/${SAMPLE}

    mkdir -p $SAMPLE_OUT

    echo "=============================="
    echo "=== Processing $SAMPLE ==="
    echo "=============================="

    # ==============================
    # 检查输入
    # ==============================
    if [ ! -f "$CONTIG" ]; then
        echo "❌ Missing contig"
        continue
    fi

    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
        echo "❌ Missing reads"
        continue
    fi

    # ==============================
    # 1. 构建 bowtie2 index
    # ==============================
    if [ ! -f "$SAMPLE_OUT/${SAMPLE}_index.1.bt2" ]; then
        echo "→ Building index"
        bowtie2-build $CONTIG $SAMPLE_OUT/${SAMPLE}_index
    else
        echo "→ Index exists, skip"
    fi

    # ==============================
    # 2. Mapping + BAM
    # ==============================
    if [ ! -f "$SAMPLE_OUT/${SAMPLE}.sorted.bam" ]; then
        echo "→ Mapping"

        bowtie2 -x $SAMPLE_OUT/${SAMPLE}_index \
                -1 $R1 -2 $R2 \
                -S $SAMPLE_OUT/${SAMPLE}.sam \
                -p $THREADS \
                --very-sensitive \
                --no-unal

        echo "→ Converting & sorting BAM"

        samtools view -@ $THREADS -bS $SAMPLE_OUT/${SAMPLE}.sam | \
        samtools sort -@ $THREADS -o $SAMPLE_OUT/${SAMPLE}.sorted.bam

        samtools index $SAMPLE_OUT/${SAMPLE}.sorted.bam
        rm $SAMPLE_OUT/${SAMPLE}.sam
    else
        echo "→ BAM exists, skip mapping"
    fi

    # ==============================
    # 3. Prokka 注释
    # ==============================
    PROKKA_DIR=$SAMPLE_OUT/prokka

    if [ ! -f "$PROKKA_DIR/${SAMPLE}.gff" ]; then
        echo "→ Running Prokka"

        conda run -n prokka_env prokka $CONTIG \
            --outdir $PROKKA_DIR \
            --prefix $SAMPLE \
            --metagenome \
            --cpus $THREADS \
            --force
    else
        echo "→ Prokka exists, skip"
    fi

    # ==============================
    # 4. featureCounts
    # ==============================
    GFF=$PROKKA_DIR/${SAMPLE}.gff
    BAM=$SAMPLE_OUT/${SAMPLE}.sorted.bam
    OUT=$SAMPLE_OUT/${SAMPLE}_counts.txt

    if [ ! -f "$GFF" ]; then
        echo "❌ Missing GFF, skip featureCounts"
        continue
    fi

    if [ ! -f "$OUT" ]; then
        echo "→ Running featureCounts"

        featureCounts \
            -a $GFF \
            -o $OUT \
            -t CDS \
            -g ID \
            $BAM \
            -T $THREADS \
            -p
    else
        echo "→ counts exists, skip"
    fi

    echo "=== DONE: $SAMPLE ==="

done

echo "🎉 ALL SAMPLES FINISHED"


###### 后续在R中进行 #######
library(Biostrings)
library(gRodon)

# =====================================================================
# 1. 文件路径设置 (请替换为你 MacBook 上的实际路径)
# =====================================================================
path_to_metagenome <- "/Volumes/JunLU/AOF_1/AOF_1.ffn"
path_to_coverage <- "/Volumes/JunLU/AOF_1/AOF_1_counts.txt"

# =====================================================================
# 2. 处理基因序列与识别高表达基因 (严格参考官方范例)
# =====================================================================
genes <- readDNAStringSet(path_to_metagenome)

# 【关键步骤】：Prokka 的 .ffn 文件头通常包含基因功能注释（如 ">PROKKA_00001 50S ribosomal protein L2"）
# 必须在截断名字之前，根据官方建议找到核糖体蛋白作为 highly_expressed
highly_expressed <- grepl("ribosomal protein", names(genes), ignore.case = TRUE)

# 确保在正确提取核糖体标记后，清理序列名以匹配 counts 表格的 ID
names(genes) <- gsub(" .*", "", names(genes))

# =====================================================================
# 3. 处理丰度与计算覆盖度 (适配 featureCounts 格式)
# =====================================================================
# featureCounts 的输出第一行通常是命令注释，使用 skip=1 跳过
counts_df <- read.table(path_to_coverage, header = TRUE, skip = 1, stringsAsFactors = FALSE)

# 提取基因 ID、长度和计数 (假设只有1个样本，Counts在第7列)
gene_ids <- counts_df$Geneid
gene_lengths <- counts_df$Length
gene_counts <- counts_df[, 7]

# 计算相对覆盖度 (Relative Coverage)
# 这一步相当于官方示例中的 meandepth，消除了基因长度偏好
relative_depths <- gene_counts / gene_lengths

# 给覆盖度向量赋名
names(relative_depths) <- gene_ids

# =====================================================================
# 4. 数据对齐与排序 (对应官方范例中的 names(depths) <- ...)
# =====================================================================
# 找出共有的序列 ID，防止匹配错误
common_ids <- intersect(names(genes), names(relative_depths))

# 过滤并按相同顺序排列所有变量
genes_filtered <- genes[names(genes) %in% common_ids]
highly_expressed_filtered <- highly_expressed[names(genes) %in% common_ids]

# 严格按照基因序列的顺序提取深度，这与官方的 depths[gsub(...)] 逻辑一致
depth_of_coverage <- relative_depths[names(genes_filtered)]
head(depth_of_coverage)

# (可选) 过滤掉丰度为 0 的基因，以加速计算并减少低信噪比影响
valid_idx <- depth_of_coverage > 0
genes_filtered <- genes_filtered[valid_idx]
highly_expressed_filtered <- highly_expressed_filtered[valid_idx]
depth_of_coverage <- depth_of_coverage[valid_idx]

# =====================================================================
# 5. 运行加权宏基因组模式 (严格参考官方范例)
# =====================================================================
result <- predictGrowth(genes_filtered, 
                        highly_expressed_filtered, 
                        mode = "metagenome_v2", 
                        depth_of_coverage = depth_of_coverage)

print(result)

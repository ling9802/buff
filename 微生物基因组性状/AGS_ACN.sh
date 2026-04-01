#### 运行MicrobeCensus-计算平均基因组大小（AGS） ####
conda activate mc_env
export TMPDIR=/home/Yangly/JunLing/TMPDIR/

run_microbe_census.py -t 100 WSF_1.R1.clean.fastq,WSF_1.R2.clean.fastq MicrobeCensus_results/WSF_1.txt


#### 从cleandata中进行 ####
conda activate mc_env
mkdir -p AGS_ACN

## 将R1和R2端拼接起来，并且随机提取出2M的序列，参考的原作者
SAMPLE="WOF_3"
vsearch --fastq_mergepairs ${SAMPLE}.R1.clean.fastq \
        --reverse ${SAMPLE}.R2.clean.fastq \
        --fastqout AGS_ACN/${SAMPLE}.merged.fastq \
        --fastqout_notmerged_fwd AGS_ACN/${SAMPLE}.unmerged_R1.fastq \
        --fastqout_notmerged_rev AGS_ACN/${SAMPLE}.unmerged_R2.fastq
        
seqtk sample -s100 AGS_ACN/${SAMPLE}.merged.fastq 2000000 | seqtk seq -A > AGS_ACN/${SAMPLE}_sampled_2M.fna

# 可选预测蛋白序列
SAMPLE="WSF_3"
prodigal -i data/${SAMPLE}_sampled_2M.fna \
            -a data/${SAMPLE}_sampled_2M.faa \
            -p meta


for fna in data/*_sampled_2M.fna; do
    # 提取样本名 (比如从 data/WSF_1_sampled_2M.fna 提取出 WSF_1)
    # ${fna%_sampled_2M.fna} 表示去掉右边的后缀
    # ${fna##*/} 表示去掉左边的路径前缀 data/
    temp=${fna%_sampled_2M.fna}
    SAMPLE=${temp##*/}

    echo "---------------------------------------"
    echo "正在处理样本: ${SAMPLE}"
    echo "输入文件: ${fna}"
    
    # 3. 运行加速版 Prodigal (推荐使用我们之前说的 -t 16 加速)
    # 如果你没装 prodigal-gv，就用普通的 prodigal
    prodigal -i "$fna" \
             -a "data/${SAMPLE}_sampled_2M.faa" \
             -p meta
             
    echo "样本 ${SAMPLE} 注释完成！"
done

echo "======================================="
echo "所有 FNA 文件已成功转化为 FAA 蛋白序列！"


## 下面的代码在macbook上运行的，直接在终端的文件夹中运行
### 定义脚本名称 ###
run_ags=/Users/lemon9802/AGS-and-ACN-tools/run_ags.sh
run_acn=/Users/lemon9802/AGS-and-ACN-tools/run_acn.sh
run_ags_fgsrs=/Users/lemon9802/AGS-and-ACN-tools/run_ags_fgsrs.sh

### 运行脚本 - AGS - ACN ###
sample="WSF_1"

"${run_ags_fgsrs}" data/${sample}_sampled_2M.fna ${sample}_ags_output \
  --min_length 100 \
  --sample_name "${sample}" \
  --verbose t \
  --overwrite t \
  --nslots 2 \
  --save_complementary_data t

# 利用已经预测好的蛋白序列直接得到AGS，但是ACN不行
"${run_ags}" \
  data/${sample}_sampled_2M.fna \
  data/${sample}_sampled_2M.faa \
  ${sample}_ags_output2 \
  --sample_name "${sample}" \
  --verbose t \
  --overwrite t \
  --nslots 2 \
  --save_complementary_data t

"${run_acn}" \
  ${sample}_ags_output/${sample}_FBL.fna \
  ${sample}_ags_output/${sample}_ags.tsv \
  ${sample}_acn_output \
  --sample_name "${sample}" \
  --evalue 1e-10 \
  --min_identity 80 \
  --min_length 30 \
  --verbose t \
  --overwrite t \
  --nslots 2 \
  --save_complementary_data t

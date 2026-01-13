#!/bin/bash
# 宏基因组SPAdes批量组装脚本（新手友好版）
# 作者：VX公众号：GutMicroLab
# ===================== 只需要改这3个变量！=====================
clean_data_dir="../02_cleandata"   # 干净数据目录（对应上面的02_cleandata）
result_dir="./spades_result"       # 组装结果输出目录（对应上面的spades_result）
sample_list="../samplelist"        # 样本列表路径（和质控时的samplelist一致）
# ==============================================================

# 1. 检查目录是否存在，不存在就自动创建（避免报错）
if [ ! -d $clean_data_dir ]; then
    echo "干净数据目录 $clean_data_dir 没找到！请检查路径"
    exit 1  # 路径错了就退出脚本，避免白跑
fi

if [ ! -d $result_dir ]; then
    mkdir -p $result_dir  # 自动创建结果目录
    echo "已创建组装结果目录：$result_dir"
fi

# 2. 循环处理每个样本（核心部分）
for sample in `cat $sample_list`; do
    # 记录开始时间和进度（日志里能看到）
    start_time=$(date "+%Y-%m-%d %H:%M:%S")
    echo "开始组装样本 $sample，时间：$start_time" >> $result_dir/spades_assembly.log 2>&1

    # 检查当前样本的干净数据是否存在（避免数据缺失导致报错）
    r1_file="$clean_data_dir/${sample}.R1.clean.fastq"
    r2_file="$clean_data_dir/${sample}.R2.clean.fastq"
    if [ ! -f $r1_file ] || [ ! -f $r2_file ]; then
        echo "样本 $sample 缺少干净数据（R1/R2文件不存在），跳过组装" >> $result_dir/spades_assembly.log 2>&1
        continue  # 跳过这个样本，继续处理下一个
    fi

    # 执行SPAdes组装（关键参数逐行解释）
    spades.py \
        --meta \
        --only-assembler \
        -1 $r1_file \
        -2 $r2_file \
        --threads 120 \
        --memory 300 \
        -o $result_dir/${sample}_spades \
        >> $result_dir/spades_assembly.log 2>&1

    # 检查组装是否成功
    if [ $? -eq 0 ]; then
        end_time=$(date "+%Y-%m-%d %H:%M:%S")
        echo "样本 $sample 组装完成，时间：$end_time" >> $result_dir/spades_assembly.log 2>&1
    else
        echo "样本 $sample 组装失败！请查看日志" >> $result_dir/spades_assembly.log 2>&1
        continue  # 失败就跳过，不影响其他样本
    fi
done

# 组装全部完成后提示
echo "所有样本组装任务结束！结果在 $result_dir，日志在 $result_dir/spades_assembly.log"
#!/bin/bash
# MEGAHIT 批量组装脚本（meta-sensitive + 功能注释友好）

clean_data_dir="../02_cleandata"
result_dir="./megahit_result"
sample_list="../samplelist"

# 创建结果目录
mkdir -p "$result_dir"

# 循环组装
while read sample; do
    echo "开始组装样本 $sample：$(date)" >> "$result_dir/megahit_assembly.log" 2>&1

    megahit \
        -1 "$clean_data_dir/${sample}.R1.clean.fastq" \
        -2 "$clean_data_dir/${sample}.R2.clean.fastq" \
        -o "$result_dir/${sample}_megahit" \
        --continue \
        --presets meta-sensitive \
        -t 100 \
        --memory 300 \
        >> "$result_dir/megahit_assembly.log" 2>&1

    if [ $? -eq 0 ]; then
        echo "样本 $sample 组装完成：$(date)" >> "$result_dir/megahit_assembly.log" 2>&1
    else
        echo "样本 $sample 组装失败：$(date)" >> "$result_dir/megahit_assembly.log" 2>&1
    fi
done < "$sample_list"

echo "MEGAHIT 组装全部结束！结果文件在 $result_dir"
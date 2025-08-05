##### 复制以下python代码作为一个py脚本，命名为run_md_analysis.py
#!/usr/bin/env python3
import os
import subprocess
import argparse
import numpy as np
from collections import defaultdict
import random


def run_md_on_sample(faa_path, out_root):
    sample_name = os.path.basename(faa_path).replace('.faa','')
    sample_outdir = os.path.join(out_root, sample_name)
    os.makedirs(sample_outdir, exist_ok=True)
    print(f"\n>>> Processing sample: {sample_name}")

    # Clean up any existing intermediate dirs/files
    for item in ['DB', 'DB_h', 'DB.index', 'DB.lookup', 'Clu', 'Clu.tsv', 'res', 'res.m8', 'tmp']:
        path = os.path.join(sample_outdir, item)
        if os.path.isdir(path):
            subprocess.run(['rm', '-rf', path])
        elif os.path.exists(path):
            os.remove(path)

    try:
        # MMseqs2 steps
        subprocess.run(["mmseqs", "createdb", faa_path, f"{sample_outdir}/DB"], check=True)
        subprocess.run(["mmseqs", "cluster", f"{sample_outdir}/DB", f"{sample_outdir}/Clu", f"{sample_outdir}/tmp"], check=True)
        subprocess.run(["mmseqs", "createtsv", f"{sample_outdir}/DB", f"{sample_outdir}/DB", f"{sample_outdir}/Clu", f"{sample_outdir}/Clu.tsv"], check=True)

        # Ensure old res/tmp are removed before search
        for d in ['res', 'tmp']:
            p = os.path.join(sample_outdir, d)
            if os.path.isdir(p): subprocess.run(['rm','-rf',p])

        subprocess.run(["mmseqs", "search", f"{sample_outdir}/DB", f"{sample_outdir}/DB", f"{sample_outdir}/res", f"{sample_outdir}/tmp"], check=True)
        result = subprocess.run(["mmseqs", "convertalis", f"{sample_outdir}/DB", f"{sample_outdir}/DB", f"{sample_outdir}/res", f"{sample_outdir}/res.m8"], cwd=sample_outdir)
        if result.returncode != 0:
            print(f"⚠️ convertalis failed for {sample_name}, skipping.")
            return

        # Parse clustering
        clu_file = os.path.join(sample_outdir, 'Clu.tsv')
        clusters = defaultdict(list)
        with open(clu_file) as f:
            for l in f:
                q,t = l.strip().split('\t',1)
                clusters[q].append(t)

        # Parse alignments
        aln_file = os.path.join(sample_outdir, 'res.m8')
        pairsim = {}
        with open(aln_file) as f:
            for l in f:
                q,t,pid = l.split('\t')[0:3]
                pairsim[(q,t)] = float(pid)

        # Compute distances
        clustdist = defaultdict(list)
        for q, members in clusters.items():
            for m in members:
                if (q,m) in pairsim:
                    clustdist[q].append(1 - pairsim[(q,m)])

        # Diversity metrics
        O = sum(len(v) for v in clusters.values())
        P = len(clustdist)
        CDlens = [len(v) for v in clustdist.values()]
        total_links = np.sum(CDlens)
        H = [ (l/total_links)*np.log(l/total_links) for l in CDlens ]
        S = [ (l/total_links)**2 for l in CDlens ]
        MD = [ 1 + (np.sum(d)/len(d)) for d in clustdist.values() ]

        # Write summary
        summary = os.path.join(sample_outdir, 'diversity_summary.txt')
        with open(summary, 'w') as f:
            f.write(f"Sample: {sample_name}\n")
            f.write(f"Total protein-encoding genes: {O}\n")
            f.write(f"Protein Richness: {P}\n")
            f.write(f"Shannon Diversity: {(-np.sum(H)):.4f}\n")
            f.write(f"Simpson Evenness: {(1-np.sum(S)):.4f}\n")
            f.write(f"Log10 Protein Dissimilarity: {np.log10(np.sum(MD)):.4f}\n")
            f.write(f"Metagenomic Diversity Index: {(1/O)*np.sum(MD):.4f}\n")

        print(f"✅ Finished {sample_name}, summary at {summary}")

    except subprocess.CalledProcessError as e:
        print(f"❌ Error processing {sample_name}: {e}")


def run_batch(faa_list_file, out_root):
    with open(faa_list_file) as f:
        for line in f:
            faa = line.strip()
            if faa.endswith('.faa') and os.path.exists(faa):
                run_md_on_sample(faa, out_root)
            else:
                print(f"❌ Skipping invalid path: {faa}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Batch MD diversity analysis')
    parser.add_argument('-l','--faa-list', required=True, help='List of .faa file paths')
    parser.add_argument('-o','--out-dir',   required=True, help='Root output directory')
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    run_batch(args.faa_list, args.out_dir)

####### 截止到107 line


######### 运行python脚本，-o可以指定输出位置，faa_list.txt 文件按照每行放置faa（蛋白文件）的文件路径即可
python run_md_analysis.py   -l faa_list.txt   -o /home/Yangly/JunLing/ISSM_data_2025/metagenome/MD_analysis/results_xu

# 也可后台运行脚本，log文件默认输出到当前文件夹中
nohup python run_md_analysis.py \
  -l faa_list.txt \
  -o /home/Yangly/JunLing/ISSM_data_2025/metagenome/MD_analysis/results_xu &

######### 提取跑出来的diversity_summary.txt放到另外一个文件夹中
### 将以下新建成一个bash文件即可，然后进行bash xxxx.sh即可
#!/bin/bash
# 设置路径
input_dir="/home/Yangly/JunLing/ISSM_data_2025/metagenome/MD_analysis/results_xu"
output_dir="/home/Yangly/JunLing/ISSM_data_2025/metagenome/MD_analysis/summary_xu"

# 创建输出文件夹（如果不存在）
mkdir -p "$output_dir"

# 遍历所有子文件夹
for folder in "$input_dir"/*/; do
    # 获取当前子文件夹的名称（不含路径）
    folder_name=$(basename "$folder")
    
    # 原始文件路径
    src_file="$folder/diversity_summary.txt"

    # 目标路径，重命名为子文件夹名.txt
    dest_file="$output_dir/${folder_name}.txt"

    # 如果文件存在，则复制并重命名
    if [[ -f "$src_file" ]]; then
        cp "$src_file" "$dest_file"
        echo "提取：$folder_name"
    else
        echo "跳过（文件不存在）：$folder_name"
    fi
done

######### 将得到的所有的diversity_summary.txt文件，整理合并到一个文件中
### 将以下新建成一个bash文件即可，然后进行bash xxxx.sh即可

# 定义变量，方便修改
INPUT_DIR="/home/Yangly/JunLing/ISSM_data_2025/metagenome/MD_analysis/summary_xu"
OUTFILE="/home/Yangly/JunLing/ISSM_data_2025/metagenome/MD_analysis/summary_xu/summary_all.tsv"

# 1. 输出表头（会覆盖已存在的同名文件）
echo -e "Sample\tTotal_protein_encoding_genes\tProtein_Richness\tShannon_Diversity\tSimpson_Evenness\tLog10_Protein_Dissimilarity\tMetagenomic_Diversity_Index" \
  > "${OUTFILE}"

# 2. 循环遍历 .txt 文件，抽取字段并追加
for f in "${INPUT_DIR}"/*.txt; do
  awk -F': ' 'BEGIN{ OFS="\t" }
    /^Sample:/                       { sample = $2 }
    /^Total protein-encoding genes:/ { tp = $2 }
    /^Protein Richness:/             { pr = $2 }
    /^Shannon Diversity:/            { sd = $2 }
    /^Simpson Evenness:/             { se = $2 }
    /^Log10 Protein Dissimilarity:/  { ld = $2 }
    /^Metagenomic Diversity Index:/  { mdi = $2 }
    END {
      print sample, tp, pr, sd, se, ld, mdi
    }
  ' "$f" >> "${OUTFILE}"
done

# 3. 查看结果前几行
head -n 5 "${OUTFILE}"

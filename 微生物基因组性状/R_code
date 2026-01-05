##### GC content （%）分析 #####
# Calculate GC content (可以分析contigs和bins的数据)

# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("Biostrings")
library(Biostrings)
# install.packages("metaCluster",
#                 repos = c("https://diprosinha.r-universe.dev",
#                           "https://cran.r-project.org"))
library(metaCluster)
library(seqinr)

library(seqinr)

filenames <- list.files("/Volumes/JL/4.assemble/fa_files", pattern = "*.fa", full.names = TRUE)
gc_list <- list()

for (f in filenames) {
  cat("Processing file:", f, "\n")
  
  # 不使用 as.string，让 seqinr 默认读为字符向量
  seqs <- read.fasta(file = f, seqtype = "DNA")
  
  gc_values <- sapply(seqs, GC)
  
  gc_list[[basename(f)]] <- gc_values
}

gc_list[[1]]

# 创建空数据框存储每个文件的平均GC
gc_avg_df <- data.frame(
  file = character(),
  mean_GC = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_along(gc_list)) {
  gc_values <- gc_list[[i]]       # 当前文件的GC值向量
  file_name <- names(gc_list)[i]  # 文件名
  
  if (is.null(file_name)) {
    file_name <- paste0("file_", i)
  }
  
  # 计算平均GC
  mean_gc <- mean(gc_values, na.rm = TRUE)
  
  # 添加到数据框
  gc_avg_df <- rbind(gc_avg_df, data.frame(file = file_name, mean_GC = mean_gc))
}
gc_avg_df

# 保存为CSV
write.csv(gc_avg_df, "/Volumes/JL/4.assemble/GC_content_average_per_file.csv", row.names = FALSE)

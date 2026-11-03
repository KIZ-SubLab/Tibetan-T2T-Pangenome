library(data.table)

# ================= 配置参数 =================
input_file      <- "00-TRGT_result/Merged.tsv"
output_file     <- "trds_batch_results.tsv" 
trds_threshold  <- 0
num_splits      <- 50  
tmp_dir         <- "tmp_workdir" 

if(!dir.exists(tmp_dir)) dir.create(tmp_dir)

# ================= 1. 读取并切分原始数据 =================
cat("正在读取大文件并切分成", num_splits, "份...\n")
df <- fread(input_file, header = FALSE, sep = "\t")
total_rows <- nrow(df)

# 计算每份的行数
rows_per_split <- ceiling(total_rows / num_splits)

for (i in 1:num_splits) {
  start_row <- (i - 1) * rows_per_split + 1
  end_row   <- min(i * rows_per_split, total_rows)
  
  if (start_row > total_rows) break
  
  # 提取子集并保存为临时输入文件
  sub_df <- df[start_row:end_row, ]
  fwrite(sub_df, file = file.path(tmp_dir, paste0("split_in_", i, ".tsv")), sep = "\t", col.names = FALSE)
}
rm(df); gc() # 释放内存

# ================= 2. 动态生成子 R 脚本并同时启动 =================
cat("正在生成子脚本并并发启动...\n")

for (i in 1:num_splits) {
  sub_in  <- file.path(tmp_dir, paste0("split_in_", i, ".tsv"))
  sub_out <- file.path(tmp_dir, paste0("split_out_", i, ".tsv"))
  sub_r_script <- file.path(tmp_dir, paste0("task_", i, ".R"))
  
  script_content <- paste0('
library(waddR)
library(data.table)

input_file <- "', sub_in, '"
output_file <- "', sub_out, '"
trds_threshold <- ', trds_threshold, '

cat("正在读取临时 TSV 矩阵数据...\\n")
df <- fread(input_file, header = FALSE, sep = "\\t")

group1_cols <- 2:71
group2_cols <- 72:ncol(df)

calculate_trds_for_row <- function(row_index, data_frame) {
  current_row <- data_frame[row_index, ]
  site_id <- as.character(current_row[[1]]) 

  g1_vectors <- round(as.numeric(current_row[, ..group1_cols]))
  g2_vectors <- round(as.numeric(current_row[, ..group2_cols]))
  
  g1_vectors <- g1_vectors[!is.na(g1_vectors)]
  g2_vectors <- g2_vectors[!is.na(g2_vectors)]

  if (length(g1_vectors) == 0 || length(g2_vectors) == 0) {
    return(data.frame(TR_ID = site_id, TRDS = NA, stringsAsFactors = FALSE))
  }
  score <- waddR::wasserstein_metric(g1_vectors, g2_vectors, p = 2)
  print(score)

  return(data.frame(TR_ID = site_id, TRDS = score, stringsAsFactors = FALSE))
}

total_rows <- nrow(df)
results_list <- vector("list", total_rows)

system.time({
  for (i in 1:total_rows) {
    results_list[[i]] <- calculate_trds_for_row(i, df)
  }
})

final_results <- do.call(rbind, results_list)
final_results <- final_results[!is.na(final_results$TRDS), ]
final_results <- final_results[order(-final_results$TRDS), ]

filtered_results <- final_results[final_results$TRDS >= trds_threshold, ]
fwrite(filtered_results, file = output_file, sep = "\\t", row.names = FALSE)
')
  
  # 写入临时 R 脚本文件
  writeLines(script_content, sub_r_script)

  log_file <- file.path(tmp_dir, paste0("task_", i, ".log"))
  
  if (.Platform$OS.type == "windows") {
    system(paste0("start /B Rscript ", sub_r_script, " > ", log_file, " 2>&1"), wait = FALSE)
  } else {
    system(paste0("Rscript ", sub_r_script, " > ", log_file, " 2>&1 &"), wait = FALSE)
  }
}

# ================= 3. 监控等待所有任务结束 =================

while (TRUE) {
  Sys.sleep(5) 
  completed_files <- list.files(tmp_dir, pattern = "split_out_.*\\.tsv")
  cat(paste0("\r当前已完成子任务: ", length(completed_files), " / ", num_splits))
  
  if (length(completed_files) == num_splits) {
    cat("\n所有子任务计算完成！\n")
    break
  }
}

# ================= 4. 合并所有临时结果并过滤排序 =================
cat("正在合并 50 个分块的结果...\n")
out_files <- list.files(tmp_dir, pattern = "split_out_.*\\.tsv", full.names = TRUE)

# 用 data.table 高效读取并合并
all_results <- rbindlist(lapply(out_files, fread))

cat("正在进行最终的清洗、排序和筛选...\n")
all_results <- all_results[!is.na(TRDS)]
all_results <- all_results[TRDS >= trds_threshold]
all_results <- all_results[order(-TRDS)]

# 写出最终结果
fwrite(all_results, file = output_file, sep = "\t", row.names = FALSE)

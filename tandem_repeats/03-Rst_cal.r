library(data.table)
library(parallel)

# ==========================================
# 1. 读取数据与严格的类型转换
# ==========================================
cat("正在读取输入文件...\n")
raw_data <- fread("/home/ffzhou/APG/APG_xizang/TRGT/Rst_workflow/00-TRGT_result/pSTR.TRGT.cn.tsv")

locus_names <- raw_data[[1]]
matrix_data <- as.matrix(raw_data[, -1, with = FALSE])
matrix_data[matrix_data == "."] <- NA
class(matrix_data) <- "numeric"

# 拆分人群 (保持原始矩阵结构)
popA_mat <- matrix_data[, 1:70]
popB_mat <- matrix_data[, 71:ncol(matrix_data)]
n_loci <- length(locus_names)

# ==========================================
# 2. 核心计算函数（包含下采样与100次迭代取平均）
# ==========================================
calc_rst_subsample <- function(i) {
  # 提取当前位点原始数据
  a_orig <- popA_mat[i, ]
  b_orig <- popB_mat[i, ]
  
  # 设定下采样目标：与PopA个体数等量 (70)
  target_n <- length(a_orig) 
  
  # 设定抽样重复次数（推荐100次以保证稳定）
  n_bootstrap <- 100
  
  # 用来存放每次抽样成功计算的 Rst 和 状态
  rst_iterations <- numeric()
  he_a_iterations <- numeric()
  he_b_iterations <- numeric()
  
  # 为了不让单次抽样因为缺失值而全盘皆输，我们循环抽样 n_bootstrap 次
  # 每次抽样时，对 PopA 和 PopB 均抽取等量样本(70个)，然后计算
  # 设置局部随机种子，防止并行计算时多核获取相同随机数
  set.seed(12345 + i)
  
  for (b in 1:n_bootstrap) {
    # 随机抽取 index 
    idx_a <- sample(1:length(a_orig), target_n, replace = FALSE)
    idx_b <- sample(1:length(b_orig), target_n, replace = FALSE)
    
    a <- a_orig[idx_a]
    b <- b_orig[idx_b]
    
    a_clean <- a[!is.na(a)]
    b_clean <- b[!is.na(b)]
    
    # 过滤缺失率 > 50%
    if (length(a_clean) / target_n < 0.5 || length(b_clean) / target_n < 0.5) {
      next
    }
    
    # 计算两个群体整体的 He
    combined_clean <- c(a_clean, b_clean)
    freq <- table(combined_clean) / length(combined_clean)
    he <- 1 - sum(freq^2)
    
    # 严格过滤：整体杂合度 > 0.1
    if (he <= 0.1) {
      next
    }
    # 计算 Rst
    v_a <- var(a_clean)
    v_b <- var(b_clean)
    sw <- mean(c(v_a, v_b))
    st <- var(c(a_clean, b_clean))
    
    if (is.na(st) || st == 0) {
      rst <- 0
    } else {
      rst <- (st - sw) / st
      if (is.na(rst) || is.nan(rst) || rst < 0) rst <- 0
    }
    
    # 记录单次迭代结果
    rst_iterations <- c(rst_iterations, rst)
  }
  
  # 判断最终状态：如果100次抽样里有超过一半成功通过QC，则保留该位点
  if (length(rst_iterations) < (n_bootstrap * 0.5)) {
    # 无法通过QC，推测是低杂合度或高缺失位点
    return(list(Rst = NA_real_, Status = "Filtered_LowHe_or_Missing"))
  } else {
    # 通过QC，对多次抽样结果取平均值
    return(list(Rst = mean(rst_iterations), Status = "Calculated"))
  }
}

# ==========================================
# 3. 多核并行计算
# ==========================================
num_cores <- min(48, detectCores() - 1) 
cat("开始等量下采样迭代计算（100次重复），共", n_loci, "个位点...\n")

results_list <- mclapply(1:n_loci, calc_rst_subsample, mc.cores = num_cores)

rst_values <- sapply(results_list, function(x) x$Rst)
status_values <- sapply(results_list, function(x) x$Status)

# ==========================================
# 4. 构建结果表并导出
# ==========================================
output_df <- data.table(
  Locus = locus_names,
  Rst = rst_values,
  QC_Status = status_values
)

filtered_df <- output_df[QC_Status == "Calculated"]

fwrite(output_df, "./strict_he_subsample_all.txt", sep = "\t", quote = FALSE)
fwrite(filtered_df, "./strict_he_subsample_filtered.txt", sep = "\t", quote = FALSE)

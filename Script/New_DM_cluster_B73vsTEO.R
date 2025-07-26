# ----------------------------
# 数据准备
# ----------------------------
# 加载必要的包
library(tibble)
library(dplyr)
library(openxlsx)  # 用于写入Excel
library(effsize)   # 用于计算效应量（可选）

# 读取数据
SM_BT <- readRDS(file.path(getwd(), "Git_data/BT_combined_MEM_objects.rds"))

# 提取数据并转换
B73_SM_neg <- SM_BT[[1]] %>% rownames_to_column(var = "label")
B73_SM_pos <- SM_BT[[2]] %>% rownames_to_column(var = "label")
Teo_SM_neg <- SM_BT[[3]] %>% rownames_to_column(var = "label")
Teo_SM_pos <- SM_BT[[4]] %>% rownames_to_column(var = "label")

# 合并数据
b73_merged <- inner_join(B73_SM_neg, B73_SM_pos, by = c("label", "Sample", "Cluster"))
teo_merged <- inner_join(Teo_SM_neg, Teo_SM_pos, by = c("label", "Sample", "Cluster"))
combined_df <- bind_rows(b73_merged, teo_merged)

# 清理环境
rm(SM_BT, B73_SM_neg, B73_SM_pos, Teo_SM_neg, Teo_SM_pos, b73_merged, teo_merged)

# 设置Sample为因子
combined_df$Sample <- factor(combined_df$Sample, levels = c("B73", "Teo"))

# ----------------------------
# 统计分析：Wilcoxon检验 + 中位数差异
# ----------------------------
# 按Cluster分组
byClust <- split(combined_df, combined_df$Cluster)
cl_names <- names(byClust)

# 初始化进度条
pb <- txtProgressBar(min = 0, max = length(cl_names), style = 3)

# 对每个Cluster进行检验
res_list <- lapply(seq_along(cl_names), function(i) {
  setTxtProgressBar(pb, i)  # 更新进度条
  cl <- cl_names[i]
  df <- byClust[[cl]]
  mzs <- names(df)[!(names(df) %in% c("label", "Sample", "Cluster"))]  # 排除非数值列
  
  # 对每个代谢物进行检验
  res <- lapply(mzs, function(mz) {
    x <- df[df$Sample == "B73", mz]
    y <- df[df$Sample == "Teo", mz]
    
    # 跳过无效数据
    if (length(x) == 0 || length(y) == 0 || all(is.na(x)) || all(is.na(y))) {
      return(data.frame(
        mz = mz,
        p.value = NA_real_,
        median_diff = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    
    # Wilcoxon检验
    p <- wilcox.test(x, y, exact = FALSE)$p.value
    
    # 计算中位数差（B73 - Teo）
    median_diff <- median(x, na.rm = TRUE) - median(y, na.rm = TRUE)
    
    # 可选：计算效应量（Cliff's Delta）
    delta <- cliff.delta(x, y)$estimate
    
    data.frame(
      mz = mz,
      p.value = p,
      median_diff = median_diff,
      cliff_delta = delta,  # 可选
      stringsAsFactors = FALSE
    )
  })
  
  # 合并当前Cluster的结果
  res_df <- do.call(rbind, res)
  res_df$Cluster <- cl
  res_df
})

close(pb)  # 关闭进度条

# 合并所有Cluster的结果
unpaired_wilcox_fast <- do.call(rbind, res_list)

# FDR校正
unpaired_wilcox_fast$p.adj <- p.adjust(unpaired_wilcox_fast$p.value, method = "fdr")

# ----------------------------
# 结果筛选与标记方向
# ----------------------------
# 筛选显著结果（p.adj < 0.001）
sig <- unpaired_wilcox_fast[unpaired_wilcox_fast$p.adj < 0.001 & !is.na(unpaired_wilcox_fast$p.adj), ]

# 标记表达方向
sig$Direction <- ifelse(
  sig$median_diff > 0, 
  "B73_Higher", 
  ifelse(sig$median_diff < 0, "Teo_Higher", "No_Difference")
)

# 按Cluster和p值排序
sig <- sig[order(sig$Cluster, sig$p.adj), ]

# ----------------------------
# 结果导出
# ----------------------------
# 写入Excel文件
write.xlsx(
  sig,
  file.path(getwd(), "Git_data/KEGG/Diff_metabolites_with_direction.xlsx"),
  sheetName = "Significant_Metabolites",
  rowNames = FALSE
)

# 按Cluster分组导出（可选）
# sig_list <- split(sig, sig$Cluster)
# write.xlsx(
#   sig_list,
#   file.path(getwd(), "Git_data/KEGG/Diff_metabolites_by_cluster.xlsx")
# )

# ----------------------------
# 结果预览
# ----------------------------
cat("显著差异代谢物数量：", nrow(sig), "\n")
print(head(sig))

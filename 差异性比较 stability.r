# contributor: Daliang Ning
rm(list = ls())
library(vegan)
library(dplyr)

setwd("/Volumes/xiaowang/zhoushan/新分析/diversity/差异性分析 组成稳定性")


phyto_stab <- read.csv("phyto_stabm_reduced.csv", header = TRUE, sep = ",", fill = TRUE)
prok_stab <- read.csv("prok_stabm_reduced.csv", header = TRUE, sep = ",", fill = TRUE)

phyto_stab$type <- "phyto"
prok_stab$type <- "prok"

# 整合数据
combined_data <- rbind(phyto_stab, prok_stab)

# 检查列名是否一致
if (!all(colnames(phyto_stab) == colnames(prok_stab))) {
  stop("Column names of phyto_stab and prok_stab are not identical.")
}

# 移除不需要的列（例如X列）
combined_data <- combined_data %>% select(-X)

# 进行t检验或Wilcoxon秩和检验
zeta_columns <- colnames(combined_data)[grepl("^Zeta", colnames(combined_data))]

results <- list()

for (zeta_col in zeta_columns) {
  # 提取当前Zeta列的数据
  zeta_data <- combined_data %>% select(zeta_col, type)
  
  # 移除NA值
  zeta_data <- na.omit(zeta_data)
  
  # 分组
  phyto_group <- zeta_data %>% filter(type == "phyto") %>% pull(zeta_col)
  prok_group <- zeta_data %>% filter(type == "prok") %>% pull(zeta_col)
  
  # 检验
  test_result <- tryCatch({
    wilcox.test(phyto_group, prok_group)
  }, error = function(e) {
    return(list(p.value = NA))
  })
  
  results[[zeta_col]] <- test_result$p.value
}

# 显示结果
results_df <- data.frame(Zeta_Column = names(results), P_Value = unlist(results))

results_df$Adjusted_P <- p.adjust(results_df$P_Value, method = "fdr")  # FDR调整
print(results_df)
write.csv(results_df,file = "sig_stab.csv")

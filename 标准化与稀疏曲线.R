rm(list = ls())
setwd("/yourworkpath")

#-----------------------------------------------------
library(dplyr)
library(ggplot2)
library(stringr)
library(ggplot2)
library(cowplot)
library(grid)
library(dplyr)
library(vegan)
library(corrplot)
library(readxl)

#Input-----------------------------------------------------
#ASVtable <- read.table("ASVtable_raw.txt", header = TRUE, sep = "\t", row.names = 1,fill = TRUE)
ASVtable <- read.table("taxonomy.txt", header = TRUE, sep = "\t", row.names = 1,fill = TRUE)

#Discard
exclude_categories <- c("D_4__Mitochondria", "D_3__Chloroplast")
ASVtable <- ASVtable %>%
  rowwise() %>%  # 对每一行进行操作
  filter(!any(c_across(c(Kingdom, Phylum, Order, Family, Genus, Class)) %in% exclude_categories))  # According to your database
numeric_cols <- names(ASVtable)[sapply(ASVtable, is.numeric)]
ASVtable1 <- ASVtable %>%
  mutate(across(where(is.numeric), ~ifelse(. < 10, 0, .))) #According to your own data
ASVtable2 <- ASVtable1[rowSums(ASVtable1[, numeric_cols], na.rm = TRUE) > 0, ]

ASVtable3 <- ASVtable2%>% select (numeric_cols)


library(vegan)	
library(ggplot2)	
library(doBy)	#用于分组统计
library(ggalt)	#用于绘制拟合曲线


##定义函数
#计算多种 Alpha 多样性指数，结果返回至向量
#各子函数的用法详见 http://blog.sciencenet.cn/blog-3406804-1179983.html
alpha_index <- function(x, method = 'chao1', tree = NULL, base = exp(1)) {
  if (method == 'richness') result <- rowSums(x > 0)	#丰富度指数
  else if (method == 'chao1') result <- estimateR(x)[2, ]	#Chao1 指数
  else if (method == 'ace') result <- estimateR(x)[4, ]	#ACE 指数
  else if (method == 'shannon') result <- diversity(x, index = 'shannon', base = base)	#Shannon 指数
  else if (method == 'simpson') result <- diversity(x, index = 'simpson')	#Gini-Simpson 指数
  else if (method == 'pielou') result <- diversity(x, index = 'shannon', base = base) / log(estimateR(x)[1, ], base)	#Pielou 均匀度
  else if (method == 'gc') result <- 1 - rowSums(x == 1) / rowSums(x)	#goods_coverage
  else if (method == 'pd' & !is.null(tree)) {	#PD_whole_tree
    pd <- pd(x, tree, include.root = FALSE)
    result <- pd[ ,1]
    names(result) <- rownames(pd)
  }
  result
}

#根据抽样步长（step），统计每个稀释梯度下的 Alpha 多样性指数，结果返回至列表
alpha_curves <- function(x, step, method = 'shannon', rare = NULL, tree = NULL, base = exp(1)) {#chao1
  x_nrow <- nrow(x)
  if (is.null(rare)) rare <- rowSums(x) else rare <- rep(rare, x_nrow)
  alpha_rare <- list()
  
  for (i in 1:x_nrow) {
    step_num <- seq(0, rare[i], step)
    if (max(step_num) < rare[i]) step_num <- c(step_num, rare[i])
    
    alpha_rare_i <- NULL
    for (step_num_n in step_num) alpha_rare_i <- c(alpha_rare_i, alpha_index(x = rrarefy(x[i, ], step_num_n), method = method, tree = tree, base = base))
    names(alpha_rare_i) <- step_num
    alpha_rare <- c(alpha_rare, list(alpha_rare_i))
  }
  
  names(alpha_rare) <- rownames(x)
  alpha_rare
}

#otu <- read.delim('otu_table.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
#otu <- read.csv("yt.csv", row.names = 1)
otu <- t(ASVtable3)
#以 2000 步长（step=2000）为例统计
chao1_curves <- alpha_curves(otu, step = 1000, method = 'shannon') #步长可以根据实际情况调整#chao1

#获得 ggplot2 作图文件
plot_chao1 <- data.frame()
for (i in names(chao1_curves)) {
  chao1_curves_i <- (chao1_curves[[i]])
  chao1_curves_i <- data.frame(rare = names(chao1_curves_i), alpha = chao1_curves_i, sample = i, stringsAsFactors = FALSE)
  plot_chao1 <- rbind(plot_chao1, chao1_curves_i)
}

rownames(plot_chao1) <- NULL
plot_chao1$rare <- as.numeric(plot_chao1$rare)
plot_chao1$alpha <- as.numeric(plot_chao1$alpha)

#ggplot2 作图
yt_chao1 <- ggplot(plot_chao1, aes(rare, alpha, color = sample)) +
  geom_line() +
  labs(x = 'Number of sequences', y = 'Shannon index', color = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = min(rowSums(otu)), linetype = 2) +
  scale_x_continuous(breaks = seq(0, 200000, 20000), labels = as.character(seq(0, 200000, 20000)))+
  scale_y_continuous(breaks = seq(0, 6, 1),
                     labels = as.character(seq(0, 6, 1))) +
  coord_cartesian(ylim = c(0, 6))  # 限制y轴的范围为0到5
yt_chao1#这一步是为了让你知道什么时候电脑运算完成上方所有任务
ggsave("yt_chao1.pdf", yt_chao1, width = 16, height = 5)
ggsave("prokaryotes_chao1.svg", yt_chao1, width = 16, height = 5)
min_value <- min(rowSums(otu))
print(min_value)

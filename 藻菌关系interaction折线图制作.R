rm(list = ls())
setwd("/Volumes/xiaowang/zhoushan/新分析/diversity/16S18SV91_interaction/network/PR2_0.6 0.05/edge/bipartitate/ridgeline")
#https://mp.weixin.qq.com/s/5osR3humOEehvpGG2Xxx1A
library(ggplot2)
library(ggridges)
library(reshape2)
library(ggsci)
library(dplyr)
library(tidyr)
library(tibble)
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(reshape2)
library(ggprism)
library(ggalt)
library(ggplot2)

table1 <- read.table("连接数量表.txt", header = TRUE, sep = "\t", row.names = 1, fill = TRUE)

table <- table1[grepl("D_3__", rownames(table1)), ]

row_sums <- rowSums(table, na.rm = TRUE)

col_sums <- colSums(table, na.rm = TRUE)

relative_abundance <- table / col_sums
relative_abundance <- sweep(table, 2, col_sums, "/")
test <- relative_abundance %>% 
  arrange(desc(row_sums)) %>% 
  head(100)

# 根据总和排序并选取总和最高的10行
top_10 <- relative_abundance %>% 
  arrange(desc(row_sums)) %>% 
  head(25)

top_10 <- rownames_to_column(top_10, var = "RowNames")


long_format <- top_10 %>%
  gather(month, value, -RowNames)

table <- rownames_to_column(table, var = "RowNames")


top_10_long <- pivot_longer(top_10, cols = -RowNames, names_to = "Sample", values_to = "Value")

samples <- c("S2111", "S2201", "S2205", "S2207", "S2209", "S2211", "S2301", "S2303")#, "S2306"

# 转换为因子，并确保因子水平按特定顺序排列

mapping <- c(S2111 =1, S2201=2, S2205=3, S2207=4, S2209=5, S2211=6, S2301=7, S2303=8)#, S2306=9
top_10_long$sample <- mapping[top_10_long$Sample]


top_10_long <- top_10_long %>%
  separate(RowNames, into = c("Type","Prokaryote" ), sep = " D_3__")

top_10_long$Prokaryote <- gsub("^D_3__", "", top_10_long$Prokaryote)
top_10_long <- top_10_long %>%
  mutate(Value = Value * 10)

top_10_long <- top_10_long %>%
  filter(!(Type %in% c("Mamiellophyceae", "Bacillariophyceae", "Cryptophyceae:nucl")))


peakcol = c("#9ECAE1", "#2171B5", "#999999", "#FF7F0E")  # 添加一个橙色，具有良好的对比度


p2c <- ggplot(top_10_long, aes(x = sample, y = Prokaryote, height = Value, color = Type,fill = Type)) +
  geom_ridgeline(stat = "identity", scale = 1.5,  show.legend = T, alpha = 0.2) + # 设置透明度
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = peakcol, limits = c('Dinophyceae', 'Syndiniales',"Mediophyceae","Cryptophyceae"))+
  scale_color_manual(values = peakcol, limits = c('Dinophyceae', 'Syndiniales', "Mediophyceae", "Cryptophyceae")) +  # 设置边框颜色
  scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") 
samples <- c("2111-2201", "2201-2205", "2205-2207", "2207-2209", 
             "2209-2211", "2211-2301", "2301-2303", "2303-2306","cc")



p2c1 <- p2c +
  scale_x_continuous(expand = c(0, 0), breaks = 1:9, labels = samples) +
  theme_minimal() +
  theme(
   # panel.background = element_rect(fill = "white", color = "black", size = 1), # 添加外边框
    panel.grid.major.x = element_line(color = "lightgray", linetype = "dashed"), # 添加垂直灰色虚线
    panel.grid.minor.x = element_blank(), # 移除次要垂直网格线
    panel.grid.major.y = element_line(color = "lightgray"), # 添加水平灰色实线
    axis.title.x = element_text(), 
    axis.text.x = element_text(),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt") # 增加图形四周的边距
  ) +
  labs(x = '', y = 'Links')
#p2c1 <-p2c+ 
  #scale_x_continuous(expand = c(0, 0))+#修改x轴刻度，这个也可以在theme()中输入axis.ticks.X=element_blank()取消x轴刻度
  #theme_minimal()+ #预设主题
  #theme(axis.title.x=element_blank(),axis.text.x.bottom=element_blank())+#取消设定的x轴标题和文本
  #labs(x = 'Samping month', y = 'Interactions')#设置x轴和y轴的标题
p2c1


  #scale_x_continuous(expand = c(0, 0), breaks = 1:9, labels = samples) + # 确保x轴显示并设置标签

  #theme(axis.title.x = element_text(), axis.text.x = element_text()) + # 确保x轴标题和文本可见



ggsave("Interactions.svg", p2c1, width = 8, height = 10)
ggsave("Interactions.png", p2c1, width = 8, height = 10)






rm(list = ls())
setwd("/Volumes/xiaowang/zhoushan/新分析/diversity/16S/NMDS/GAM_nmds")
library(openxlsx)
library(vegan)
library(ggplot2)
df <- read.table("sequences_rarefied_20240222.txt", header = TRUE, sep = "\t",fill = TRUE)
#排除分类学数据
df1 <- df[, -c(1:8)]#原来是8 为了匹配GAM 去掉前两个月数据
df3 <- df[, -c(1:16)]

df2=t(df1)
dfGroup = read.delim("sample-metadata_modifiedNMDS.txt",row.names = 1)

distance <- vegdist(df2, method = 'bray')
#NMDS排序分析，k = 2预设两个排序轴
nmds <- metaMDS(distance, k = 2)

#dfNmds<-metaMDS(df,distance="bray",k = 2)
data = data.frame(nmds$points)

data <- merge(data, dfGroup, by = "row.names", all = TRUE)  # "row.names"指定了合并依据为行名
#data <- data[-c(1:8), ]
# 重命名合并后的第一列，它现在包含了原始的行名
# 假设merged_data是合并后的数据框，且第一列是想要设置为行名的列
rownames(data) <- data[,1]  # 将第一列的值设置为行名
data <- data[,-1]  # 移除第一列
data$Order <- c(1, 3, 2, 4, 1, 2, 4, 3, 1, 6, 3, 2, 4, 5, 1, 2, 3, 5, 4, 6, 2, 3, 1, 6, 5, 4, 1, 4, 3, 2, 6, 5, 1,3, 4, 5, 2, 6, 1, 6, 2, 3, 5, 4, 1, 3, 2, 6, 4, 5)
data <- data[order(data$Group), ]#先按分组排序
#1, 3, 2, 4, 1, 2, 4, 3, 
data <- data[order(data$Group, data$Order), ]

#ZHIXING
#p <- ggplot(data, aes(x = MDS1, y = MDS2, color = Group, group = Group, fill = Group)) +
#  geom_point(size = 1.5) +
#  stat_ellipse(aes(x = MDS1, y = MDS2, fill = Group, group = Group, color = Group), 
#               level = 0.95,    # 设定置信度水平，这里是95%
#               type = "norm",    # 椭圆类型，"norm"适用于正态分布的数据
#               alpha = 0.3, 
#               linetype = "longdash", 
#               linewidth = 0.5) +
#  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"), 
#        panel.grid = element_blank(), 
#        panel.background = element_rect(fill = 'white')) +
#  labs(subtitle = paste("stress=", round(nmds$stress, 3), sep=""))
#
#p

#设置随机种子
set.seed(123)
#基于bray-curtis距离进行PERMANOVA分析
adonis <-  adonis2(df2 ~ Group, data = dfGroup, permutations = 999, method = "bray")
#基于bray-curtis距离进行anosim分析
anosim = anosim(df2, dfGroup$Group, permutations = 999, distance = "bray")

print(adonis)
print(anosim)

adonis_text <- paste(paste("Adonis  =", round(adonis$R2, 2)), "***")[1]
anosim_text <- paste(paste("Anosim  =", round(anosim$statistic, 2)), "***")


#添加GAM
plot(nmds, type = "n")  # 'type = "n"'创建一个只有坐标轴的空图，适合作为后续添加元素的基础
df3 <- df[, -c(1:16)]

#df3  <- df3[ , -4]


df4=t(df3)
distance_gam <- vegdist(df4, method = 'bray')
#NMDS排序分析，k = 2预设两个排序轴
nmds_gam <- metaMDS(distance_gam, k = 2)
#dfNmds<-metaMDS(df,distance="bray",k = 2)
data_gam = data.frame(nmds_gam$points)
data_gam <- merge(data_gam, dfGroup, by = "row.names", all = TRUE)  # "row.names"指定了合并依据为行名


data_gam <- data_gam[-c(1:8), ]
#data_gam <- data_gam[-4, ]


env <- read.table("totalenvs202402.txt", header = TRUE, sep = "\t",fill = TRUE,row.names = 1)
#env_T <- env[-3, ]
env_T <- env
#好像不是数值类型？
#env_T$T <- as.numeric(env_T$T)

#str(env_T)
fit <- ordisurf(nmds_gam ~ T, env_T, col = "blue", add = TRUE, select = FALSE, method = "GCV.Cp")

summary(fit)

ordi.grid <- fit$grid
#合并数据
ordi.mite <- expand.grid(x = ordi.grid$x, y = ordi.grid$y)
ordi.mite$z <- as.vector(ordi.grid$z)
ordi.mite.na <- data.frame(na.omit(ordi.mite))




p2 <- ggplot() +
  geom_point(data = data, 
             aes(x = MDS1, y = MDS2, group = Group, fill = Group,shape = location), # 添加shape映射fill = Group,
             size = 3,pch = 21) + #color = 'transparent', pch = 21, 
  scale_shape_manual(values = c("farming" = 16, "beach" = 1)) +
  geom_polygon(data = data, aes(x = MDS1, y = MDS2, fill = Group), 
               alpha = 0.3, linetype = "longdash", size = 0.5)+
  stat_contour(data = ordi.mite.na, 
               aes(x = x, y = y, z = z, color = ..level..), size = 1.5) +
  scale_colour_continuous(high = '#066502', low = '#c4fd70') +
 labs(x = 'NMDS1', y = 'NMDS2', color = 'T', fill = "Group", shape = "Location") + # 更新图例标题
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"), 
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'white')) +
  labs(                     # 在副标题处添加stress
    subtitle = paste(paste("Stress=", round(nmds$stress, 3), sep = " ", adonis_text), anosim_text),
  )+
  annotate("text", x = Inf, y = -1, label = "Deviance explained = 96.1% ***", 
           hjust = 1.05, vjust = 2, size = 4, color = "black", angle = 0)


 # guides(color = guide_legend(override.aes = list(shape = NA))) # 如果需要，保留这行来调整图例显示
p2
#  geom_text(aes(label = row.names(df2)), vjust = -1, hjust = 0.5, size = 3, color = "black") +  # 添加站位名标签
ggsave("NMDS 16S_T.svg", p2, height = 3, width = 5)
ggsave("NMDS 16S_T.png", p2, height = 3, width = 5)

write.xlsx(data, "16S_NMDS_data_T.xlsx",rowNames = TRUE)

#ggtitle(paste(paste(stress_text, adonis_text), anosim_text))
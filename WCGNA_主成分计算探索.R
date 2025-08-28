library(ggplot2)
library(reshape2)
library(forcats)
library(igraph)
library(tidyverse)
library(WGCNA)
#BiocManager::install("preprocessCore")
library(ggClusterNet)
library(viridis)
library(dplyr)
library(stringr)
library(ggnewscale)
library(readxl)
rm(list = ls())
setwd("/Volumes/xiaowang/zhoushan/新分析/diversity/16S18SV91_interaction/WGCNA/ASV_new")

pro16S <- read.table("sequences_rarefied_20240222.txt", header = TRUE, sep = "\t" )
pro16S <- pro16S[pro16S[, 1] != "", ]
taxonomy <- pro16S[, 1:8]
taxonomy <- taxonomy[taxonomy[, 1] != "", ]
rownames(taxonomy) <- taxonomy[, 1]

rownames(pro16S) <- pro16S[, 1]
pro16S <- pro16S[, -(1:8)]
otu_presence16S <- pro16S > 0  # 创建一个逻辑矩阵，标记OTU是否在每个样本中出现（计数大于0）
otu_count16S <- rowSums(otu_presence16S)  # 对每一行（OTU）计算True（即出现）的总数
pro16S_1 <- pro16S[otu_count16S >=6, ]# #(ncol(pro16S) / 3)
#pro16S_1 <-pro16S
pro16S_t <- t(pro16S_1)
pro16S_df <- as.data.frame(pro16S_t)

phyto18S <- read.table("PR2_phyto_ASV_ENV.txt", header = TRUE, sep = "\t" , row.names = 1)

num_cols <- ncol(phyto18S)
phyto_taxonomy <- phyto18S %>%
  select(1:7)

# 删除空字符串所在的行＆# 设置第一列为行名，taxonomy完成
phyto_taxonomy <- phyto_taxonomy[phyto_taxonomy[, 1] != "", ]

phyto_taxonomy <- phyto_taxonomy %>%
  tibble::rownames_to_column(var = "ASV") %>%
  select(ASV, everything())
rownames(phyto_taxonomy) <- phyto_taxonomy$ASV

#后面需要根据每个月月份来进行asvtable的筛选

phyto_asvtable <- phyto18S %>%
  select(8:46) 


#

combined_tax <- rbind(taxonomy, phyto_taxonomy)


otu_presence18S <- phyto_asvtable > 0  # 创建一个逻辑矩阵，标记OTU是否在每个样本中出现（计数大于0）
otu_count18S <- rowSums(otu_presence18S)  # 对每一行（OTU）计算True（即出现）的总数
phyto18S_1 <- phyto_asvtable[otu_count18S >= 6, ]#(ncol(phyto_asvtable) / 3)

#phyto18S_1 <- phyto_asvtable

phyto18S_t <- t(phyto18S_1)
phyto18S_df <- as.data.frame(phyto18S_t)

# 将行名转换为数据框的一个列
pro16S_df$ID <- rownames(pro16S_df)
phyto18S_df$ID <- rownames(phyto18S_df)


# 根据新创建的ID列合并数据
merged_data <- merge(pro16S_df, phyto18S_df, by="ID", all=TRUE)
# 将ID列设置为行名
rownames(merged_data) <- merged_data$ID
# 从数据框中移除ID列
merged_data$ID <- NULL

write_csv(merged_data, "merged_data.csv")

library(WGCNA)

datExpr <- as.data.frame(merged_data) # 转置数据，如果需要

# 选择一个合适的软阈值
powers = c(c(1:10), seq(from = 12, to=30, by=2))

#powers <- c(c(1:10), 
#seq(from = 12, 
#    to = 30,
#    by = 2))
# 使用pickSoftThreshold函数，计算每个power值对应的scale free topology模型拟合程度和平均连接度，
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# 查看sft$powerEstimate来确定最佳的软阈值

par(mfrow = c(1,2))
cex1 = 0.9
# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
# 网络越符合无标度特征 (non-scale)
pdf("my_softthreshold plot.pdf")#最好手动保存

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="turquoise")
# 筛选标准。R-square=0.85
abline(h=0.85,col="red")

# Soft threshold与平均连通性
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
     cex=cex1, col="turquoise")

dev.off()


# 构建无尺度网络
softPower <- 8 # 例如，这里使用6作为软阈值
adjacency <- adjacency(datExpr, power = softPower)
dimnamesAdjacency <- dimnames(adjacency)
# 转换为拓扑重叠矩阵 它考虑了网络中各节点的共享邻居。
TOM <- TOMsimilarity(adjacency)
dimnames(TOM) <- dimnamesAdjacency
dissTOM <- 1-TOM

# 基因树和模块检测 deepSplit = 1
geneTree <- hclust(as.dist(dissTOM), method = "average")
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 1, pamRespectsDendro = FALSE)

# 将模块颜色分配给基因树
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

dynamicColors2 = as.data.frame(dynamicColors)


# 绘制模块检测结果 这边最终是划分出了10个模块
svg("tree.svg", width = 7, height = 7)
p <- plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                         dendroLabels = FALSE, hang = 0.03,
                         addGuide = TRUE, guideHang = 0.05)

# 关闭SVG设备
dev.off()

#moduleLabels_MB = net_MB$colors

#moduleColors_MB = labels2colors(net_MB$colors)

#MEs_MB = net_MB$MEs;
#geneTree_MB = net_MB$dendrograms[[1]];
#save(MEs_MB, moduleLabels_MB, moduleColors_MB, geneTree_MB, 
#     file = "/Users/anni_djurhuus/Documents/Projects/MBON/m2w/google_drive/data/MEGAN6/network_data/MB_network_construction.RData")


# 使用图形设备函数保存图形

library(stringr)

moduleNetworkThreshold2 <- TOM
#moduleNetworkThreshold[moduleNetworkThreshold < 0.1] <- 0
moduleGraph2 <- graph.adjacency(moduleNetworkThreshold2, mode = "undirected", weighted = TRUE ,diag = FALSE)# weighted是NULL还是TRUE？


#sub_net_layout4 <- layout_in_circle(moduleGraph2)
sub_net_layout4 <- layout_with_fr(moduleGraph2, niter=999,grid = 'nogrid')
data2 = as.data.frame(sub_net_layout4)
data2$OTU = igraph::get.vertex.attribute(moduleGraph2)$name
colnames(data2) = c("X1","X2","elements")

row.names(data2) = data2$elements
dat2 = data.frame(
  OTU = V(moduleGraph2)$name,
  X1 = data2$X1,
  X2 = data2$X2
)


list(data2, dat2, moduleGraph2)
node2 <- data.frame(name = dat2[,1])
row.names(node2) <- dat2[,1]
node2$ModuleColor <- dynamicColors
nodes2 = nodeadd(plotcord =node2,otu_table = merged_data,tax_table = combined_tax)
head(nodes2)
#-----计算边#--------

temm1 = moduleNetworkThreshold2 %>%
  tidyfst::mat_df() %>%
  dplyr::filter(row != col) %>%
  dplyr::rename(OTU_1 = row, OTU_2 = col, weight = value) %>%
  dplyr::filter(weight != 0)

temm2 = temm1 %>% dplyr::left_join(data2,by = c("OTU_1" = "elements")) %>%
  dplyr::rename(Y1 = X2)
#head(tem2)
temm3 = data2 %>%
  dplyr::rename(Y2 = X2,X2 = X1) %>%
  dplyr::right_join(temm2,by = c("elements" = "OTU_2")) %>%
  dplyr::rename(OTU_2 = elements)

edge3 = temm3 %>%
  dplyr::mutate(
    cor = ifelse(weight > 0,"+","-")
  )
colnames(edge)[8] = "cor"


temm2 = dat2 %>% 
  dplyr::select(OTU) %>%
  dplyr::right_join(edge3,by =c("OTU" = "OTU_1" ) ) %>%
  dplyr::rename(OTU_1 = OTU)
head(temm2)

temm3 = dat2 %>% 
  dplyr::select(OTU) %>%
  dplyr::right_join(edge3,by =c("OTU" = "OTU_2" ) ) %>%
  dplyr::rename(OTU_2 = OTU)
head(temm3)

temm4 = temm2 %>%inner_join(temm3)
#head(tem4)

temm5 <- temm4 %>%
  filter((str_detect(OTU_1, "phyto") & !str_detect(OTU_2, "phyto")) | 
           (!str_detect(OTU_1, "phyto") & str_detect(OTU_2, "phyto")))
#temm6 <- temm5 %>%
#  filter(weight > 0.1)

temm6 <- temm5 %>%
  inner_join(node2, by = c("OTU_1" = "name")) %>%
  rename(Module_1 = ModuleColor) %>%
  # 再次内部连接来为OTU_2找到模块
  inner_join(node2, by = c("OTU_2" = "name")) %>%
  rename(Module_2 = ModuleColor) %>%
  # 保留那些两个OTU属于同一模块的边
  filter(Module_1 == Module_2)

#csv_filename7 <- sprintf("edge%s.csv", month_pattern)
#write.table(tem4, file=csv_filename7, sep=",", quote=FALSE)

fixed_Classes <- c("D_2__Alveolata", "D_2__Stramenopiles", "D_2__Cryptomonadales","D_2__Alphaproteobacteria","D_2__Nitrososphaeria","D_2__Bacteroidia","D_2__Gammaproteobacteria")

phylum1_colors = viridis(length(fixed_Classes) + 1)#top_phyla替换为fixed_Classes
phylum_colors <- c("#9C27B0","#FFFF00","#FFA500", "#4CAF50","#FF0000", "#35B779FF","#0000FF","#757575")
phylum_colors[length(phylum_colors)] <- "#757575"  # 将最后一个颜色设置为 "#C1C1C1"# 为排名前7的门和 "Others" 分配颜色
phylum_color_map <- setNames(phylum_colors, c(fixed_Classes, "Others"))
# 更新 nodes 数据框，为排名前7的门和 "Others" 分配颜色-这边其实分类也得统一，而不应该是排名前7的物种
nodes2$Color <- ifelse(nodes2$Class %in% fixed_Classes, phylum_color_map[nodes$Class], "#757575")

nodes2 <- nodes2 %>%
  mutate(phylum_select = ifelse(Class %in% fixed_Classes, Class, "Others"))
# 创建包含Class和对应颜色的数据框
phylum_color <- data.frame(Class = c(fixed_Classes, "Others"), Color = phylum_colors)
# 映射颜色到 ggplot
col2 = phylum_color$Color
names(col2) = phylum_color$Class

#提取分类注释结果
select_nodes2 <- nodes2 %>%
  select(name, Class, Class)
nodes_dat2 <- data2 %>%
  left_join(select_nodes2, by = c("elements" = "name"))
nodes_dat2$ModuleColor <- dynamicColors
nodes_dat2 <- nodes_dat2 %>%
  mutate(kingdoM = if_else(str_detect(elements, "phyto"), "Eukaryote", "Prokaryote"))

# 获取唯一值
unique_values <- unique(temm6$Module_2)
unique_values1 <- unique(nodes_dat2$ModuleColor)

p1 <- ggplot() + 
  geom_curve(aes(x = X1, y = Y1, xend = X2, yend = Y2, color = Module_2),
             data = temm6, size = 0.2, alpha = 0.05, curvature = 0.2) +
  scale_colour_manual(values = c("green" = "#00FF00", "blue" = "#0096c7","brown"="#C2A878","grey"= "grey","turquoise" ="turquoise","yellow"="#f4d35e",
                                 "red"="red","pink"="pink", "black"="black","magenta"="magenta")) 
p2 = p1 +
  new_scale_color() +
  geom_point(aes(X1, X2, color = ModuleColor,shape = kingdoM), data = nodes_dat2, size = 0.4,alpha = 0.3) +
  scale_colour_manual(values = c("green" = "#00FF00", "blue" = "#0096c7","brown"="#C2A878","grey"= "grey","turquoise" ="turquoise","yellow"="#f4d35e",
                                 "red"="red","pink"="pink", "black"="black","magenta"="magenta"))  + 
  scale_shape_manual(values = c("Eukaryote" = 17, "Prokaryote" = 16)) +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +#这一步可以更改背景颜色和边框
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
p2

ggsave("./cs2.svg",p2,width = 4,height = 5)#原6 5
ggsave("./cs2.png",p2,width = 4,height = 5)
write.table(nodes_dat2, file = "相关节点.txt", sep = "\t", quote = FALSE, row.names = FALSE)


prokaryote_df <- filter(nodes_dat2, kingdoM == "Prokaryote")

eukaryote_df <- filter(nodes_dat2, kingdoM == "Eukaryote")
write.table(prokaryote_df, file = "相关节点Prokaryote.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(eukaryote_df, file = "相关节点Eukaryote.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(phyto18S_finalout, file = "节点Eukaryote.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(pro16S_finalout, file = "节点Prokaryote.txt", sep = "\t", quote = FALSE, row.names = FALSE)


datExpr_MB <- datExpr
nGenes_MB = ncol(datExpr)
nSamples_MB = nrow(datExpr)
# Recalculate MEs with color labels
#MEs0_MB = moduleEigengenes(datExpr_MB, moduleColors_MB)$eigengenes
##MEs_MB = orderMEs(MEs0_MB)
#moduleTraitCor_MB = cor(MEs_MB, environmentalData, use = "p");
#moduleTraitPvalue_MB = corPvalueStudent(moduleTraitCor_MB, nSamples_MB);

#textMatrix_MB =  paste(signif(moduleTraitCor_MB, 2), "\n(",
   #                    signif(moduleTraitPvalue_MB, 1), ")", sep = "");
#dim(textMatrix_MB) = dim(moduleTraitCor_MB)
#par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
#svg("heatmap.svg", width = 7, height = 7)
#labeledHeatmap(Matrix = moduleTraitCor_MB,
          #     xLabels = names(environmentalData),
         #      yLabels = names(MEs_MB),
        #       ySymbols = names(MEs_MB),
       #        colorLabels = FALSE,
      #         colors = blueWhiteRed(50),
     #          textMatrix = textMatrix_MB,
    #           setStdMargins = FALSE,
   #            cex.text = 0.5,
  #             zlim = c(-1,1),
 #              main = paste("Module-trait relationships"))

#dev.off()



MEs <- moduleEigengenes(datExpr, colors = dynamicColors)$eigengenes
write.table(MEs, file = "MEs.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

# 绘制模块间相关性图
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap",
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = TRUE,
                      xLabelsAngle = 90)

#network heatmap plot
TOMplot(dissTOM, MEs$dendrograms, dynamicColors,
        main = "Network heatmap plot")

# 假设您有一个名为environmentalData的DataFrame，其行与datExpr的行对应 undone
environmentalData  <- read.table("env_LN_doc.txt", header = TRUE, sep = "\t", row.names = 1 )

MEs1 <- MEs[-c(1:10),]
MEs1 <-MEs1[-2,]
datExpr1 <- datExpr[-c(1:10),-12,]
datExpr1<-datExpr1[-2,]

moduleTraitCorrelation <- cor(MEs1, environmentalData, use = "p")

geneTraitSignificance <- cor(datExpr1, environmentalData, use = "p") 
library(pheatmap)

# 计算距离矩阵，这里假设你想要基于模块与环境因子的相关性
#row_dist <- dist(moduleTraitCorrelation)
#col_dist <- dist(t(moduleTraitCorrelation))
#dev.off()
# 使用pheatmap绘制热图
#pheatmap(moduleTraitCorrelation,
#         clustering_distance_rows = row_dist,
#         clustering_distance_cols = col_dist,
#         cluster_rows = TRUE,  # 如果你想要基于距离聚类行
#         cluster_cols = TRUE,  # 如果你想要基于距离聚类列
#         show_rownames = TRUE,  # 根据需要调整
#         show_colnames = TRUE,  # 根据需要调整
#         angle_col = 45,  # 调整列名显示的角度
#         main = "Species vs Environmental Factors Correlation Heatmap")
#导出moduleTraitCorrelation
#write.table(moduleTraitCorrelation, file = "moduleTraitCorrelation.txt", sep = "\t", row.names = TRUE, col.names = NA)

# 选择与特定环境因子（如温度）相关性最高的模块 undone
#selectedModule <- names(MEs)[which.max(abs(moduleTraitCorrelation[, "NO3"]))]这个代码有点问题





# 假设datExpr的每列代表一个taxon，environmentalData包含环境因子
correlationMatrix <- cor(datExpr1, environmentalData, use = "complete.obs")
pValueMatrix <- matrix(nrow = ncol(datExpr1), ncol = ncol(environmentalData))

for (i in 1:ncol(datExpr1)) {
  for (j in 1:ncol(environmentalData)) {
    testResult <- cor.test(datExpr1[, i], environmentalData[, j], method = "pearson")
    pValueMatrix[i, j] <- testResult$p.value
  }
}


# 将相关性矩阵和P值矩阵合并为长格式的数据框
longData <- melt(correlationMatrix)
longData$pValue <- melt(pValueMatrix)$value

longData[is.na(longData)] <- 0
# 绘制热图
ggplot(longData, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  geom_text(aes(label = ifelse(pValue < 0.05, "*", "")), color = "black") +
  theme_minimal() +
  xlab("Taxon") +
  ylab("Environmental Factor")

row_hclust <- hclust(dist(correlationMatrix), method = "complete")
col_hclust <- hclust(dist(t(correlationMatrix)), method = "complete")
# 使用聚类结果排序
sorted_correlationMatrix <- correlationMatrix[row_hclust$order, col_hclust$order]
sorted_pValueMatrix <- pValueMatrix[row_hclust$order, col_hclust$order]
ggsave("my_plot.pdf", plot = last_plot(), device = "pdf")





# 初始化p值矩阵
pValueMatrix <- matrix(nrow = ncol(MEs1), ncol = ncol(environmentalData))

# 计算相关性和p值
for (i in 1:ncol(MEs1)) {
  for (j in 1:ncol(environmentalData)) {
    testResult <- cor.test(MEs1[, i], environmentalData[, j], method = "pearson")
    pValueMatrix[i, j] <- testResult$p.value
  }
}

# 转换为长格式并准备绘图数据
longData <- reshape2::melt(cor(MEs1, environmentalData, use = "complete.obs"))
longData$pValue <- reshape2::melt(pValueMatrix)$value

# 绘制热图并添加显著性标注
penv<- ggplot(longData, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  geom_text(aes(label = ifelse(pValue < 0.05, "*", "")), color = "black") +
  theme_minimal() +
  xlab("Module Eigengene") +
  ylab("Environmental Factor")
penv

library(ggplot2)

penv <- ggplot(longData, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#07BEB8", high = "#EF5B5B", mid = "white", midpoint = 0, 
                       name = "Correlation", limits = c(-max(abs(longData$value)), max(abs(longData$value)))) +
  geom_text(aes(label = ifelse(pValue < 0.05, "*", "")), color = "black", size = 3) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  xlab("") +
  ylab("") +
  ggtitle("")

print(penv)

ggsave("env_LN_doc.svg",plot = penv, width = 10, height = 10)# 原来10 10
ggsave("env_LN_doc.png",plot = penv, width = 10, height = 10)
# 假设dynamicColors是分配给基因的模块颜色
moduleColors <- labels2colors(dynamicMods)
names(moduleColors) <- names(datExpr)  # 确保基因或者taxon的名称是对应的

#子节点的分析
#子节点的分析
#子节点的分析
#子节点的分析
#子节点的分析
# 提取该模块的基因/物种 undone
moduleGenes <- names(datExpr)[dynamicColors == "turquoise"]
# 提取 turquoise 模块基因的表达数据
module_expression_data <- datExpr[, moduleGenes]
# 查看数据的基本统计信息
summary(module_expression_data)
# 绘制箱线图查看 turquoise 模块基因表达的分布
boxplot(module_expression_data, las = 2, main = "Expression of Turquoise Module Genes",
        ylab = "Expression", xlab = "Genes")
# 计算模块特征基因
module_eigengenes <- moduleEigengenes(datExpr, dynamicColors)$eigengenes

# 提取 turquoise 模块的特征基因表达
turquoise_eigengene <- module_eigengenes[["MEturquoise"]]
list(turquoise_eigengene)

# 使用TOM矩阵作图 undone
moduleNetwork <- TOM[moduleGenes, moduleGenes]
moduleNetworkThreshold <- moduleNetwork
#moduleNetworkThreshold[moduleNetworkThreshold < 0.1] <- 0
moduleGraph <- graph.adjacency(moduleNetworkThreshold, mode = "undirected", weighted = TRUE ,diag = FALSE)# weighted是NULL还是TRUE？


sub_net_layout3 <- layout_with_fr(moduleGraph, niter=999,grid = 'nogrid')


#layout_circles(g)
data = as.data.frame(sub_net_layout3)
data$OTU = igraph::get.vertex.attribute(moduleGraph)$name
colnames(data) = c("X1","X2","elements")

row.names(data) = data$elements
dat = data.frame(
  OTU = V(moduleGraph)$name,
  X1 = data$X1,
  X2 = data$X2
)

# 如果这段代码不是在函数中，您不需要使用return
# 直接列出对象即可，它们会在环境中保持可用状态
list(data, dat, moduleGraph)
node <- data.frame(name = dat[,1])
row.names(node) <- dat[,1]
#导入nodes,tax and otu _ here is "merged_data"  combined_tax

nodes = nodeadd(plotcord =node,otu_table = merged_data,tax_table = combined_tax)
head(nodes)
#-----计算边#--------

tem1 = moduleNetworkThreshold %>%
  tidyfst::mat_df() %>%
  dplyr::filter(row != col) %>%
  dplyr::rename(OTU_1 = row, OTU_2 = col, weight = value) %>%
  dplyr::filter(weight != 0)

tem2 = tem1 %>% dplyr::left_join(data,by = c("OTU_1" = "elements")) %>%
  dplyr::rename(Y1 = X2)
head(tem2)
tem3 = data %>%
  dplyr::rename(Y2 = X2,X2 = X1) %>%
  dplyr::right_join(tem2,by = c("elements" = "OTU_2")) %>%
  dplyr::rename(OTU_2 = elements)

edge = tem3 %>%
  dplyr::mutate(
    cor = ifelse(weight > 0,"+","-")
  )
colnames(edge)[8] = "cor"


tem2 = dat %>% 
  dplyr::select(OTU) %>%
  dplyr::right_join(edge,by =c("OTU" = "OTU_1" ) ) %>%
  dplyr::rename(OTU_1 = OTU)
head(tem2)

tem3 = dat %>% 
  dplyr::select(OTU) %>%
  dplyr::right_join(edge,by =c("OTU" = "OTU_2" ) ) %>%
  dplyr::rename(OTU_2 = OTU)
head(tem3)

tem4 = tem2 %>%inner_join(tem3)
#head(tem4)

tem5 <- tem4 %>%
  filter((str_detect(OTU_1, "phyto") & !str_detect(OTU_2, "phyto")) | 
           (!str_detect(OTU_1, "phyto") & str_detect(OTU_2, "phyto")))

#csv_filename7 <- sprintf("edge%s.csv", month_pattern)
#write.table(tem4, file=csv_filename7, sep=",", quote=FALSE)

fixed_Classes <- c("D_2__Alveolata", "D_2__Stramenopiles", "D_2__Cryptomonadales","D_2__Alphaproteobacteria","D_2__Nitrososphaeria","D_2__Bacteroidia","D_2__Gammaproteobacteria")

phylum1_colors = viridis(length(fixed_Classes) + 1)#top_phyla替换为fixed_Classes
phylum_colors <- c("#9C27B0","#FFFF00","#FFA500", "#4CAF50","#FF0000", "#35B779FF","#0000FF","#757575")
phylum_colors[length(phylum_colors)] <- "#757575"  # 将最后一个颜色设置为 "#C1C1C1"# 为排名前7的门和 "Others" 分配颜色
phylum_color_map <- setNames(phylum_colors, c(fixed_Classes, "Others"))
# 更新 nodes 数据框，为排名前7的门和 "Others" 分配颜色-这边其实分类也得统一，而不应该是排名前7的物种
nodes$Color <- ifelse(nodes$Class %in% fixed_Classes, phylum_color_map[nodes$Class], "#757575")

nodes <- nodes %>%
  mutate(phylum_select = ifelse(Class %in% fixed_Classes, Class, "Others"))
# 创建包含Class和对应颜色的数据框
phylum_color <- data.frame(Class = c(fixed_Classes, "Others"), Color = phylum_colors)
# 映射颜色到 ggplot
col2 = phylum_color$Color
names(col2) = phylum_color$Class

#提取分类注释结果
select_nodes <- nodes %>%
  select(name, Class, Class)
nodes_dat <- data %>%
  left_join(select_nodes, by = c("elements" = "name"))

nodes <- nodes %>%
  rename(PHYLUm = Class, Class = phylum_select)


nodes_dat <- nodes_dat %>%
  mutate(kingdoM = if_else(str_detect(elements, "phyto"), "Eukaryote", "Prokaryote"))


p1 <- ggplot() + 
  geom_curve(aes(x = X1, y = Y1, xend = X2, yend = Y2),
             data = tem5, size = 0.3, alpha = 0.2, curvature = 0.2, color = "turquoise") 

p2 = p1 +
  new_scale_color() +
  geom_point(aes(X1, X2, color = Class,shape = kingdoM), data = nodes_dat, size = 5) +
  scale_colour_manual(values = col2) +
  scale_shape_manual(values = c("Eukaryote" = 17, "Prokaryote" = 16)) +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +#这一步可以更改背景颜色和边框
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
p2
ggsave("./cs1.png",p2,width = 16,height = 14)
ggsave("./cs1.svg",p2,width = 16,height = 14)#聚类实在是过于离谱，感觉可以试着用整体网络做一个？但是对于单个网络来说还是没啥必要——使用igraph自己的可视化进行制作即可？
length(V(moduleGraph));length(E(moduleGraph))  
#还有节点的大小需要调整，在这时候需要知道一些degree

#-----------------------MEs的可视化

result <- data.frame()
result1 <- data.frame()
# 循环遍历每个模块列
#MEblue <- c("MEblue")
for (module in colnames(MEs)) {#colnames(MEs)
  # 创建临时数据框，仅保留当前模块列
  temp_df <- data.frame(
    station = rownames(MEs),
    Expression = MEs[[module]]
  )
  
  # 提取第二到第五个字符生成 Month 列
  temp_df$Month <- substr(temp_df$station, 2, 5)
  
  # 按 Month 分组，计算平均值和标准差
  summary_df <- temp_df %>%
    group_by(Month) %>%
    summarize(
      Mean = mean(Expression, na.rm = TRUE),
      SD = sd(Expression, na.rm = TRUE)
    )
  
  # 添加模块名称列，标识数据来源
  temp_df$Module <- module
  summary_df$Module <- module
  
  # 将结果合并到结果数据框
  result <- rbind(result, summary_df)
  result1 <- rbind(result1,temp_df)
}

module_error<-ggplot(result, aes(x = Month, y = Mean, group = Module)) +
  geom_line(aes(color = Module), size = 1) +  # 添加折线图，按 Module 分配颜色
  geom_point(aes(color = Module), size = 2) +  # 添加点，按 Module 分配颜色
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, color = "black") +  # 添加误差棒
  facet_wrap(~ Module, scales = "free_y", ncol = 2) +  # 按模块分组
  scale_y_continuous(breaks = seq(-1, 1, by = 0.2)) +  # 设置 y 轴间隔为 0.2
  scale_color_manual(values = c("MEgreen" = "#00FF00", "MEblue" = "#0096c7", "MEbrown" = "#E76F51", 
                                "MEgrey" = "grey", "MEturquoise" = "turquoise", "MEyellow" = "#f4d35e",
                                "MEred" = "red", "MEpink" = "pink", "MEblack" = "black")) +
  scale_fill_manual(values = c("MEgreen" = "#00FF00", "MEblue" = "#0096c7", "MEbrown" = "#E76F51", 
                               "MEgrey" = "grey", "MEturquoise" = "turquoise", "MEyellow" = "#f4d35e",
                               "MEred" = "red", "MEpink" = "pink", "MEblack" = "black")) +
  labs(x = "Month", y = "Temporal ASVs profile") +
  theme_minimal() +
  theme(
    legend.position = "none",  # 去掉图例
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转 x 轴标签
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)  # 添加黑色边框
  )
module_error
ggsave("module_boxplot.svg",plot = module_error,width = 4,height=10 )
ggsave("module_boxplot.jpg",plot = module_error,width = 4,height=10 )


mean_values <- result1 %>%
  group_by(Module, Month) %>%
  summarize(Mean = mean(Expression, na.rm = TRUE), .groups = 'drop')
module_boxplot <-ggplot(result1, aes(x = Month, y = Expression, fill = Module)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # 添加箱线图，不显示离群点，设置透明度
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +  # 
  geom_line(data = mean_values, aes(x = Month, y = Mean, group = Module, color = Module), size = 0.5) +  # 添加均值连接线
  geom_point(data = mean_values, aes(x = Month, y = Mean, group = Module, color = Module), size = 0.5) + # 添加均值点
  facet_wrap(~ Module, scales = "free_y", ncol = 1) +  # 按模块分组
  scale_y_continuous(breaks = seq(-0.2, 0.6, by = 0.2)) +  # 设置 y 轴范围和间隔
  labs(x = "Month", y = "Temporal ASVs profile", title = "") +
  scale_color_manual(values = c("MEgreen" = "#00FF00", "MEblue" = "#0096c7", "MEbrown" = "#E76F51", 
                                "MEgrey" = "grey", "MEturquoise" = "turquoise", "MEyellow" = "#f4d35e",
                                "MEred" = "red", "MEpink" = "pink", "MEblack" = "black")) +
  scale_fill_manual(values = c("MEgreen" = "#00FF00", "MEblue" = "#0096c7", "MEbrown" = "#E76F51", 
                               "MEgrey" = "grey", "MEturquoise" = "turquoise", "MEyellow" = "#f4d35e",
                               "MEred" = "red", "MEpink" = "pink", "MEblack" = "black")) +
  theme_minimal() +
  theme(
    legend.position = "none",  # 去掉图例
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转 x 轴标签
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)  # 添加黑色边框
  )
module_boxplot
ggsave("module_boxplot.svg",plot = module_boxplot,width = 4,height=10 )#4 10
ggsave("module_boxplot.jpg",plot = module_boxplot,width = 4,height=10 )


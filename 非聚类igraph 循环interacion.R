#try to use SparCC network  igraph方法
library(phyloseq)
library(igraph)
library(network)
library(sna)
library(tidyverse)
library(ggClusterNet)
library(readxl)
library(readr)
library(tidyfst)
library(ggnewscale)
library(viridis)
library(vegan)
library(ggforce)
library(svglite)
library(cowplot)


#devtools::install_github("eliocamp/ggnewscale")


rm(list = ls())
setwd("E:/Zhoushan eDNA/新分析/diversity/16S18SV91_interaction/network")
metadata <- read.table("metadata.txt", header = TRUE, sep = "\t", row.names = 1 )
asvtable1 <- read.table("筛选后的sliva18Sphytoplankton.txt", header = TRUE, sep = "\t" , row.names = 1)
num_cols <- ncol(asvtable1)
taxonomy1 <- asvtable1[, (num_cols - 6):num_cols]

# 删除空字符串所在的行＆# 设置第一列为行名，taxonomy完成
taxonomy1 <- taxonomy1[taxonomy1[, 1] != "", ]
rownames(taxonomy1) <- rownames(asvtable1)



#后面需要根据每个月月份来进行asvtable的筛选
asvtable1 <- asvtable1[,  1:(num_cols - 7)]
taxonomy1$ASV <- rownames(taxonomy1)

# 将新的ASV列移动到第一列位置
taxonomy1 <- taxonomy1[, c("ASV", setdiff(names(taxonomy1), "ASV"))]


asvtable <- read.table("sequences_rarefied_20240222.txt", header = TRUE, sep = "\t" )
taxonomy <- asvtable[, 1:8]

# 删除空字符串所在的行＆# 设置第一列为行名，taxonomy完成
taxonomy <- taxonomy[taxonomy[, 1] != "", ]
rownames(taxonomy) <- taxonomy[, 1]



#后面需要根据每个月月份来进行asvtable的筛选
asvtable <- asvtable[asvtable[, 1] != "", ]
rownames(asvtable) <- asvtable[, 1]
asvtable <- asvtable[, -(1:8)]

combined_asvtable <- rbind(asvtable, asvtable1)
asvtable <- combined_asvtable
combined_tax <- rbind(taxonomy, taxonomy1)
taxonomy <- combined_tax
#这一步也可以是去除小于0.1%的ASV样本
rowSums <- rowSums(asvtable)

# 将每个OTU计数除以其对应样本的总计数，转换为相对丰度的数据框relativeAbundance
relativeAbundance <- asvtable / rowSums
#对ASV表格进行排除+抽平
# 计算每个OTU出现的样本数
#relativeAbundance[relativeAbundance < 0.001] <- 0

row_totals <- rowSums(relativeAbundance)

# 保留总和大于0的行
relativeAbundance1 <- relativeAbundance[row_totals > 0, ]
plots_list <- list()
months <- c("2111", "2201", "2205", "2207", "2209", "2211", "2301", "2303", "2306")
for (month_pattern in months) {
# 筛选出至少在一半样本中出现的OTUs。这个没有问题，但是这个要分月份后再进行分析

# 构造正则表达式，以匹配列名中包含该月份及其后的站位的列
# 这里的正则表达式意味着寻找以"S2201_"开始，后面跟任意字符（代表站位）的字符串
pattern = paste0("S", month_pattern, "_\\d+")

# 使用grep函数搜索匹配模式的列名
selected_columns <- grep(pattern, colnames(relativeAbundance1), value = TRUE)

# 使用筛选出的列名来索引relativeAbundance1数据框，保留这些列
filtered_data <- relativeAbundance1[, selected_columns]
otu_presence <- filtered_data > 0  # 创建一个逻辑矩阵，标记OTU是否在每个样本中出现（计数大于0）
otu_count <- rowSums(otu_presence)  # 对每一行（OTU）计算True（即出现）的总数
relativeAbundance2 <- filtered_data[otu_count >= (ncol(filtered_data) / 2), ]
otus_to_keep <- rownames(relativeAbundance2)
# 筛选taxonomy，只保留存在于otus_to_keep中的行
filtered_taxonomy <- taxonomy[rownames(taxonomy) %in% otus_to_keep, ]
filtered_ASV <- asvtable[rownames(asvtable)%in% otus_to_keep, ]
filtered_ASV1 <- filtered_ASV[, selected_columns]


# 计算每个物种的总丰度
species_total_abundance <- rowSums(asvtable)

# 计算所有物种的总丰度
total_abundance <- sum(species_total_abundance)

# 计算每个物种的相对丰度
sp.ra <- species_total_abundance / total_abundance
sp.ra <- sp.ra[names(sp.ra) %in% otus_to_keep]
#数据整合
ps = phyloseq(sample_data(metadata),
              otu_table(as.matrix(filtered_ASV1), taxa_are_rows=TRUE),
              tax_table(as.matrix(filtered_taxonomy))#,
              # phy_tree(tree),
              # refseq(rep)
)

num_rows = nrow(asvtable)
#-提取丰度最高的指定数量的otu进行构建网络

#----------计算相关#----
#result = cor_Big_micro(ps = ps,
#                       N = num_rows,
#                       r.threshold=0.6,
#                       p.threshold=0.05,
#                       method = "pearson"
#                       )
#低于0.6的都变成0了，但是样本量少的话cor很容易非0即1？
#--提取相关矩阵

result= corMicro(ps=ps,
                 N=num_rows,
                 r.threshold=0.6,
                 p.threshold = 0.05,
                 method.scale = "rela",
                 method = "sparcc",
                 R = 20,
                 ncpus = 1
)


cor = result[[1]]
dim(cor)
head(cor)

igraph <-  graph.adjacency(cor, weighted = NULL, mode = 'undirected')
igraph <- simplify(igraph)#删除孤立节点
length(V(igraph));length(E(igraph))  
sub_net_layout <- layout_with_fr(igraph, niter=999,grid = 'nogrid')


#layout_circles(g)
data = as.data.frame(sub_net_layout)
data$OTU = igraph::get.vertex.attribute(igraph)$name
colnames(data) = c("X1","X2","elements")

row.names(data) = data$elements
dat = data.frame(
  OTU = V(igraph)$name,
  X1 = data$X1,
  X2 = data$X2
)

# 如果这段代码不是在函数中，您不需要使用return
# 直接列出对象即可，它们会在环境中保持可用状态
list(data, dat, igraph)
node <- data.frame(name = dat[,1])
row.names(node) <- dat[,1]

otu_table = as.data.frame(t(vegan_otu(ps)))
tax_table = as.data.frame(vegan_tax(ps))
nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
head(nodes)
#-----计算边#--------

tem1 = cor %>%
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
head(tem4)


csv_filename7 <- sprintf("edge%s.csv", month_pattern)
write.table(tem4, file=csv_filename7, sep=",", quote=FALSE)




#
#phylum_abundance = nodes %>% 
#  group_by(Class) %>% 
#  summarise(TotalMean = sum(mean, na.rm = TRUE)) %>%
#  arrange(desc(TotalMean))
#top_phyla = head(phylum_abundance, 7)$Class# 选取丰度排名前7的门
# 使用 dplyr 包中的 mutate 函数来创建新列
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

#生成全部分类学颜色
#unique_phylum <- unique(nodes$Class)
#phylum_colors <- viridis(length(unique_phylum))# 创建包含Class和对应颜色的数据框
#phylum_color <- data.frame(Class = unique_phylum, Color = phylum_colors)# 为NA值分配灰色
#phylum_color$Color[is.na(phylum_color$Class)] <- "#C1C1C1"

#col1 = phylum_color$Color
#names(col1) = phylum_color$Class

#提取分类注释结果
select_nodes <- nodes %>%
  select(name, Class, Class)
nodes_dat <- data %>%
  left_join(select_nodes, by = c("elements" = "name"))

nodes <- nodes %>%
  rename(PHYLUm = Class, Class = phylum_select)

center_x <- mean(range(tem4$X1))
center_y <- mean(range(tem4$Y1))

# 计算最大距离作为半径
max_radius <- max(sqrt((tem4$X1 - center_x)^2 + (tem4$Y1 - center_y)^2))


library(ggnewscale)

#p1 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = color),
#                              data = tem4, size = 1) +
#  scale_colour_manual(values = col0) 
p1 <- ggplot() + 
  geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2, color = cor), 
               data = tem4, size = 0.1,alpha = 0.2) +
  scale_colour_manual(values = c("+" = "orange", "-" = "grey"),
                      labels = c("Negative", "Positive"))
p1
p2 = p1 +
  new_scale_color() +
  geom_point(aes(X1, X2, color = Class), data = nodes_dat, size = 3) +
  scale_colour_manual(values = col2) +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white", colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

length(V(igraph));length(E(igraph)) #这边为止igraph都是对的 268个节点 151条边 
#degree_nodes <- igraph::degree(igraph)#似乎存在孤立节点

#  geom_circle(aes(x0 = center_x, y0 = center_y, r = max_radius+1), color = "black", linetype = "solid", size = 1) + # 添加圆形轮廓
# ggsave("./cs1.pdf",p1,width = 16,height = 14)
#以下是nodes颜色等通过model来注释。
#p2 = p1 +
#  new_scale_color() +
#  geom_point(aes(X1, X2,color =model), data = dat,size = 3) +
#  scale_colour_manual(values = col) +
#  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
#  theme(panel.background = element_blank()) +
#  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
#  theme(legend.background = element_rect(colour = NA)) +
# theme(panel.background = element_rect(fill = "white",  colour = NA)) +
#theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
#p2
non_zero_count <- sum(cor != 0)

# 打印结果
print(non_zero_count)
# 计算所有值都是0的行数
all_zero_rows <- sum(rowSums(cor == 0) == ncol(cor))

# 打印结果
print(all_zero_rows) #发现有两行是彻底是0，所以有两个节点是孤立节点

#节点模块化和可视化计算 这边的node和之前的node变了
result4 = nodeEdge(cor = cor)
#提取变文件
edge = result4[[1]]
#--提取节点文件
node = result4[[2]]
igraph  = igraph::graph_from_data_frame(edge, directed = FALSE, vertices = node)
degree_nodes <- igraph::degree(igraph)

length(V(igraph));length(E(igraph))

nodes_dat$Degree <- degree_nodes[match(nodes_dat$elements, V(igraph)$name)]

# 根据度数范围设置节点大小
nodes_dat$Size <- ifelse(nodes_dat$Degree <= 20, 2,
                         ifelse(nodes_dat$Degree <= 40, 4, 6))
p2 <- p2 +
  geom_point(aes(X1, X2, color = Class, size = Size), data = nodes_dat) +
  scale_size_continuous(range = c(2, 6), 
                        breaks = c(2, 4, 6), 
                        labels = c("0-20", "20-40", ">40")) +
  scale_colour_manual(values = col2)
p2
png_filename <- sprintf("network%s.png", month_pattern)
svg_filename <- sprintf("network%s.svg", month_pattern)
ggsave(png_filename,p2,width = 16,height = 14)
ggsave(svg_filename,p2,width = 16,height = 14)

plots_list[[length(plots_list) + 1]] <- p2#保存至组图中

res = ZiPiPlot(igraph = igraph,method = "cluster_walktrap")#这边是与后续的网络节点属性计算同步，那个改不了
p <- res[[1]]
taxa.roles <- res[[2]]
write.table(taxa.roles, file="keytaxa2306.csv", sep=",", quote=FALSE)
p
png_filename2 <- sprintf("key%s.png", month_pattern)
svg_filename2 <- sprintf("key%s.svg", month_pattern)
ggsave(png_filename2,p,width = 16,height = 14)
ggsave(svg_filename2,p,width = 16,height = 14)
#

# 升级后包含的网络属性更多
daaat = net_properties.2(igraph,n.hub = T)
head(daaat,n = 16)

#节点性质计算
nodepro = node_properties(igraph)
head(nodepro)
#关键ASV的挑选
hub = hub_score(igraph)$vector %>%
  sort(decreasing = TRUE) %>%
  head(5) %>%
  as.data.frame()

colnames(hub) = "hub_sca"

hub_plot <- ggplot(hub) +
  geom_bar(aes(x = hub_sca,y = reorder(row.names(hub),hub_sca)),stat = "identity",fill = "#4DAF4A")

png_filename3 <- sprintf("hub_plot%s.png", month_pattern)
svg_filename3 <- sprintf("hub_plot%s.svg", month_pattern)
ggsave(png_filename3, plot = hub_plot, width = 10, height = 8)
ggsave(svg_filename3, plot = hub_plot, width = 10, height = 8)
#导出网络属性
#daaat
csv_filename3 <- sprintf("node property%s.csv", month_pattern)
txt_filename3 <- sprintf("node property%s.txt", month_pattern)
write.table(nodepro, file=csv_filename3, sep=",", quote=FALSE)
write.table(nodepro, file=txt_filename3, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

#关于新的degree参数做一些调整
library(dplyr)

# 调整每条边，确保节点顺序一致
adjusted_edges <- edge %>%
  rowwise() %>%
  mutate(from = min(c(from, to)),
         to = max(c(from, to))) %>%
  ungroup()

# 移除重复的边
unique_edges <- adjusted_edges %>%
  distinct(from, to) %>%
  nrow()

print(unique_edges)

nodepro_df <- data.frame(degree = nodepro)
total_degree_sum <- sum(nodepro_df$degree.igraph.degree)
head(total_degree_sum)
#计算phylum的degree
phylum_degrees <- aggregate(Degree ~ Class, data = nodes_dat, sum) # 按门分类汇总度数

# 您可以查看每个门的总连接数
print(phylum_degrees)
txt_filename4 <- sprintf("phylum_degrees%s.txt", month_pattern)
write.table(phylum_degrees, file=txt_filename4, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

### 南农rubustness评估与脆弱性评估
## 鲁棒性评估robustness simulation 
#input network matrix, percentage of randomly removed species, and ra of all species
#return the proportion of species remained

sp.ra2 <- sp.ra[colSums(abs(cor)) > 0]

rm.p.list = seq(0.05,0.2,by=0.05)

rand.remov.once<-function(netRaw, rm.percent, sp.ra, abundance.weighted=T){
  #-随机挑选出一定百分比的OTU
  id.rm<-sample(1:nrow(netRaw), round(nrow(netRaw)*rm.percent))
  net.Raw=netRaw
  #这些节点和其他节点连接全部去除
  net.Raw[id.rm,]=0;  net.Raw[,id.rm]=0;   ##remove all the links to these species
  if (abundance.weighted){
    #网络矩阵乘以物种平均丰度，改变相关性值的大小
    net.stength= net.Raw*sp.ra
  } else {
    net.stength= net.Raw
  }
  # 每一个节点的平均链接数
  sp.meanInteration<-colMeans(net.stength)
  
  id.rm2<- which(sp.meanInteration<=0)  ##remove species have negative interaction or no interaction with others
  remain.percent<-(nrow(netRaw)-length(id.rm2))/nrow(netRaw)
  #for simplicity, I only consider the immediate effects of removing the
  #'id.rm' species; not consider the sequential effects of extinction of
  # the 'id.rm2' species.
  
  #you can write out the network pruned
  #  net.Raw[id.rm2,]=0;  net.Raw[,id.rm2]=0;
  # write.csv( net.Raw,"network pruned.csv")
  remain.percent
}
library(tidyverse)

rmsimu<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se)
    names(result)<-c("remain.mean","remain.sd","remain.se")
    result
  }))
}

Weighted.simu<-rmsimu(netRaw=cor, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
Unweighted.simu<-rmsimu(netRaw=cor, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

dat1<-data.frame(Proportion.removed=rep(seq(0.05,1,by=0.05),2),rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),each=20),
                 year=rep(2014,40),treat=rep("Warmed",40))

currentdat<-dat1
txt_filename5 <- sprintf("random_removal_result%s.txt", month_pattern)
write.table(currentdat, file=txt_filename5, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
##plot
library(ggplot2)

currentdat$year

ggplot(currentdat[currentdat$weighted=="weighted",], aes(x=Proportion.removed, y=remain.mean, group=treat, color=treat)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("blue","red"))+
  xlab("Proportion of species removed")+
  ylab("Proportion of species remained")+
  theme_light()+
  facet_wrap(~year, ncol=3)


ggplot(currentdat[currentdat$weighted=="unweighted",], aes(x=Proportion.removed, y=remain.mean, group=treat, color=treat)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("blue","red"))+
  xlab("Proportion of species removed")+
  ylab("Proportion of species remained")+
  theme_light()
# OK 目前好像成功地做出来robustness的图，还需要对图像稍微修改一下
#脆弱性计算
g3 = graph_from_adjacency_matrix(as.matrix(cor), mode="undirected", weighted = NULL, diag = FALSE, add.colnames = NULL) 
iso_node_id = which(igraph::degree(g3)==0)
length(V(g3));length(E(g3))  
#check node number and links
node.vul<-info.centrality.vertex(g3)
igraph2 <- igraph 
E(igraph2)$weight <- abs(E(igraph2)$weight)
node.vul<-info.centrality.vertex(igraph2)
max(node.vul)
daaat_df <- as.data.frame(daaat)
# 计算 max(node.vul) 的值
max_vul_value <- max(node.vul)

# 添加至daaat_df中
daaat_df["max_vulnerability", ] <- max_vul_value

#cohension计算
library(vegan)
# 假设您有一个物种相对丰度矩阵 relative_abundance_matrix

taxa_shuffle_model <- function(abundance_matrix, num_iterations = 1000) {
  m <- ncol(abundance_matrix)  # 物种的数量
  n <- nrow(abundance_matrix)  # 样本的数量
  shuffled_cor <- matrix(0, m, m)  # 随机化后的相关性矩阵，尺寸应为 m×m
  
  for (i in 1:num_iterations) {
    # 随机打乱物种丰度矩阵的行（每一列代表一个物种）
    shuffled_matrix <- abundance_matrix[, sample(1:m, m)]
    
    # 计算随机化矩阵的相关性
    shuffled_cor <- shuffled_cor + cor(t(shuffled_matrix))  # 确保传递给cor函数的是转置后的矩阵
  }
  
  shuffled_cor <- shuffled_cor / num_iterations
  return(shuffled_cor)
}


# 计算Cohesion
connectedness_positive <- rowMeans(cor * (cor > 0), na.rm = TRUE)
connectedness_negative <- rowMeans(cor * (cor < 0), na.rm = TRUE)

positive_cohesion <- sum(sp.ra * connectedness_positive)
negative_cohesion <- sum(sp.ra * connectedness_negative)
daaat_df["positive_cohesion", ] <- positive_cohesion
daaat_df["negative_cohesion", ] <- negative_cohesion
txt_filename6 <- sprintf("network property%s.txt", month_pattern)
csv_filename6 <- sprintf("network property%s.csv", month_pattern)
write.table(daaat_df, file=csv_filename6, sep=",", quote=FALSE)
write.table(daaat_df, file=txt_filename6, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
}

final_plot <- plot_grid(plotlist = plots_list, ncol = 3)
ggsave("combined_plots.png", plot = final_plot, width = 30, height = 20)

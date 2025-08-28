
library("RColorBrewer")
library("ggplot2")
library("reshape2")
library("dplyr")
library("melt")
rm(list = ls())
windowsFonts(A=windowsFont("Times New Roman"),
             B=windowsFont("Arial"))
# 构造模拟数据：
data <- read.table("LINKS.txt", header = TRUE, sep = "\t" ,row.names = 1)

rownames(data) <- sub("^D_2__", "", rownames(data))

relative_abundance <- data

relative_abundance$TotalAbundance <- rowSums(relative_abundance[,-ncol(relative_abundance)])

relative_abundance$rowname <- row.names(relative_abundance) # 将行名添加为数据框的一列

# 排序并获取丰度前10的物种
top_species <- relative_abundance %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 15) %>%
  .$rowname

# 将不在前10的物种标记为'Others'
relative_abundance$Species <- ifelse(relative_abundance$rowname %in% top_species, relative_abundance$rowname, 'Others')

# 对'Others'进行汇总
relative_abundance_aggregated <- relative_abundance %>%
  group_by(Species) %>%
  summarise(across(-c(rowname, TotalAbundance), sum, .names = "sum_{.col}")) %>%
  ungroup()
col_names <- names(relative_abundance_aggregated)

# 使用正则表达式替换匹配的列名
new_col_names <- gsub("^sum_", "", col_names)

# 将新列名应用到数据框
names(relative_abundance_aggregated) <- new_col_names

relative_abundance <- relative_abundance_aggregated

data.melt <- melt(relative_abundance, id.vars = "Species")
data.melt$variable <- factor(data.melt$variable, levels = unique(data.melt$variable))

data.melt<-melt(relative_abundance,id.vars='Species')
data.melt$Species<-factor(data.melt$Species,levels = rev(relative_abundance$Species))
p1 <- ggplot(data=data.melt,aes(variable,value,fill=Species))+
  geom_bar(stat="identity", position="fill",color="black", width=0.6,size=0.4)+
  scale_fill_manual(values=c("#56B4E9",'gray', '#CCEBC5', '#BC80BD', '#FCCDE5', 
                             '#B3DE69', '#FDB462', '#80B1D3', '#FB8072',
                             '#BEBADA', '#FFFFB3', '#8DD3C7',"#F52F57","#094074","#FFEE93"))+
  theme(
    axis.title=element_text(size=15,face="plain",color="black"),
    axis.text = element_text(size=15,face="plain",color="black"),
    legend.title=element_text(size=15,face="plain",color="black"),
    legend.position = "right",
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black", size = 0.4))+theme_bw()+
  theme(text=element_text(family="A",size=20))+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )
p1

library(ggalluvial)

p2 <- ggplot(data=data.melt,aes(variable,value,fill=Species,stratum = Species, alluvium = Species)) +
  geom_stratum(color="black",width=0.6,size=0.5)+
  geom_flow(alpha = 0.5) +  #绘制同类别之间的连接线
  scale_fill_manual(values=c("#56B4E9",'gray', '#CCEBC5', '#BC80BD', '#FCCDE5', 
                             '#B3DE69', '#FDB462', '#80B1D3', '#FB8072',
                             '#BEBADA', '#FFFFB3', '#8DD3C7',"#F52F57","#094074","#FFEE93")) +
  theme(
    axis.title=element_text(size=15,face="plain",color="black"),
    axis.text = element_text(size=15,face="plain",color="black"),
    legend.title=element_text(size=15,face="plain",color="black"),
    legend.position = "right",
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black", size = 0.5))+theme_bw()+
  theme(text=element_text(family="A",size=20))+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )
p2

ggsave(filename = "interaction alage-pro 冲击图.png", p2, width=20, height=10)
ggsave(filename = "interaction alage-pro 冲击图.svg", p2, width=20, height=10)
rm(list = ls())
setwd("/Volumes/xiaowang/zhoushan/新分析/diversity/16S/α-diversity")

#load Rpackage
library(vegan)
library(picante)
library(ggplot2)
library(dplyr)
library(multcompView)

#Set_function     This function is coded by Rpackage "amplicon"
#---------------------------------------------------------------------------------------------------------------------
alpha_diversity <- function(x, tree = NULL) {
  observed_species <- estimateR(x)[1, ]
  Chao1 <- estimateR(x)[2, ]
  ACE <- estimateR(x)[4, ]
  shannon <- diversity(x, index = 'shannon',base = 2)
  simpson <- diversity(x, index = 'simpson')    #注意，这里是Gini-Simpson 指数
  goods_Coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  #保留四位小数
  shannon <- sprintf("%0.4f", shannon)
  simpson <- sprintf("%0.4f", simpson)
  goods_Coverage <- sprintf("%0.4f", goods_Coverage)
  
  
  result <- data.frame(observed_species, ACE,Chao1, shannon, simpson, goods_Coverage)
  
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
    
    result <- data.frame(observed_species, ACE,Chao1, Shannon, Simpson,
                         PD_whole_tree ,goods_Coverage)
  }
  result
}

alpha_boxplot1 <- function(alpha_div, metadata, index = "richness", groupID = "Group",facet_name = "richness",
                           outlier = TRUE
) {
  
  # 交叉筛选
  idx = rownames(metadata) %in% rownames(alpha_div)
  metadata = metadata[idx,,drop=F]
  alpha_div = alpha_div[rownames(metadata),]
  
  # 提取样品组信息,默认为group可指定
  sampFile = as.data.frame(metadata[, groupID],row.names = row.names(metadata))
  # colnames(sampFile)[1] = "group"
  
  # 合并alpha_div和metadata
  df = cbind(alpha_div[rownames(sampFile),index], sampFile)
  colnames(df) = c(index,"group")
  
  # 统计各种显著性
  model = aov(df[[index]] ~ group, data=df)
  # 计算Tukey显著性差异检验
  Tukey_HSD = TukeyHSD(model, ordered = TRUE, conf.level = 0.95)
  # 提取比较结果
  Tukey_HSD_table = as.data.frame(Tukey_HSD$group)
  
  # 保存统计结果
  # 保存一个制表符，解决存在行名时，列名无法对齐的问题
  write.table(paste(date(), "\nGroup\t", groupID, "\n\t", sep=""), file=paste("alpha_boxplot_TukeyHSD.txt",sep=""),append = T, quote = F, eol = "", row.names = F, col.names = F)
  # 保存统计结果，有warning正常
  suppressWarnings(write.table(Tukey_HSD_table, file=paste("alpha_boxplot_TukeyHSD.txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))
  
  # 函数：将Tukey检验结果P值转换为显著字母分组
  # 输入文件为图基检验结果和分组
  generate_label_df = function(TUKEY, variable){
    # library(multcompView)
    # 转换P值为字母分组
    ## 提取图基检验中分组子表的第4列P adjust值
    Tukey.levels = TUKEY[[variable]][,4]
    # 方法1.multcompLetters函数将两两p值转换为字母，data.frame并生成列名为Letters的数据框
    Tukey.labels = data.frame(multcompLetters(Tukey.levels)['Letters'])
    # 方法2. 解决字母顺序相反的问题
    # library(multcomp)
    # tuk <- cld(glht(model, alternative = 'two.sided', linfct = mcp(group = 'Tukey')), sig = p, decreasing = TRUE)
    # Tukey.labels <- data.frame(Letters=tuk$mcletters$Letters, stringsAsFactors = FALSE)
    
    # 按分组名字母顺序
    ## 提取字母分组行名为group组名
    Tukey.labels$group = rownames(Tukey.labels)
    # 按组名的字母顺序排列，默认的Levels
    Tukey.labels=Tukey.labels[order(Tukey.labels$group), ]
    return(Tukey.labels)
  }
  
  # 当只有两组时，用LSD标注字母
  if (length(unique(df$group)) == 2){
    # LSD检验，添加差异组字母
    library(agricolae)
    out = LSD.test(model, "group", p.adj="none")
    stat = out$groups
    # 分组结果添入Index
    df$stat=stat[as.character(df$group),]$groups
    # 当大于两组时，用multcompView标注字母
  }else{
    # library(multcompView)
    LABELS = generate_label_df(Tukey_HSD , "group")
    df$stat=LABELS[as.character(df$group),]$Letters
  }
  
  # 设置分组位置为各组y最大值+高的5%
  max=max(df[,c(index)])
  min=min(df[,index])
  x = df[,c("group",index)]
  #y = x %>% group_by(group) %>% summarise_(Max=paste('max(',index,')',sep=""))
  y = x %>% group_by(group) %>% summarise(Max = max(.data[[index]], na.rm = TRUE))
  y=as.data.frame(y)
  rownames(y)=y$group
  df$y=y[as.character(df$group),]$Max + (max-min)*0.05
  
  if (!is.null(facet_name )) {
    df$richness <- facet_name  # 确保为每个调用传递正确的facet名称
  }
  
  if (outlier) {
    # 绘图 plotting
    p = ggplot(df, aes(x=group, y=.data[[index]], color=group)) +
      geom_boxplot(alpha=1,
                   # outlier.shape = NA,
                   # outlier.size=0,
                   size=0.7,
                   width=0.5, fill="transparent") +
      labs(x="Groups", y=paste(index, "index"), color=groupID) + theme_bw() +
      facet_wrap(~richness, scales = "free",ncol = 1) +  
      geom_text(data=df, aes(x=group, y=y, color=group, label=stat)) +
      geom_jitter(position=position_jitter(0.17), size=1, alpha=0.7)+
      theme(text=element_text(family="sans", size=7),  # 通用文本设置
            strip.text.x = element_text(face="bold", size=10, family="sans"))
    p
  } else{
    # 绘图 plotting
    p = ggplot(df, aes(x=group, y=.data[[index]], color=group)) +
      geom_boxplot(alpha=1,
                   outlier.shape = NA,
                   outlier.size=0,
                   size=0.7,
                   width=0.5, fill="transparent") +
      labs(x="Groups", y=paste(index, "index"), color=groupID) + theme_bw() +
      facet_wrap(~richness, scales = "free",ncol = 1) +  # 使用facet_wrap来分割图表
      geom_text(data=df, aes(x=group, y=y, color=group, label=stat)) +
      geom_jitter(position=position_jitter(0.17), size=1, alpha=0.7)+
      theme(text=element_text(family="sans", size=7),  # 通用文本设置
            strip.text.x = element_text(face="bold", size=10, family="sans"))
    p
  }
  
  
}

#Input your own data
#-----------------------------------Alpha_diversity_calculation------------------------------------
metadata = read.table("sample-metadata.txt", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)
# 预览元数据前3行，注意分组列名
head(metadata, n = 3)
# 将Group列转换为因子类型
metadata$Group <- as.factor(metadata$Group)

OTUtable <- read.table("ASVtable.txt", header = TRUE, sep = "\t",fill = TRUE)

otu <- OTUtable %>%
  select(9:58)#discard taxonomy annotation
otu1 <- t(otu)
otu <- otu1
alpha_div <- alpha_diversity(otu)
write.table(alpha_div, "alpha_diversity.txt", sep = "\t", quote = FALSE, row.names = TRUE,col.names = TRUE)
#---------------------------------Test of normal distribution--------------------------------------
shannon_diversity <- alpha_div[,4]
simpson_diversity <- alpha_div[,5]

str(shannon_diversity)
shannon_diversity <- as.numeric(shannon_diversity)
simpson_diversity <- as.numeric(simpson_diversity)
sink("shapiro_test_results.txt")
print("Shannon diversity Shapiro-Wilk test:")
shapiro.test(shannon_diversity)
print("Simpson diversity Shapiro-Wilk test:")
shapiro.test(simpson_diversity)
sink()

#------------------------------Visualization-------------------------------------------------------
alpha_div = read.table("alpha_diversity.txt", header=T, row.names=1, sep="\t")
alpha_div$shannon <- as.numeric(alpha_div$shannon)
#alpha_div$simpson <- as.numeric(alpha_div$simpson)
#alpha_div$Chao1 <- as.numeric(alpha_div$Chao1)
metadata1 <- metadata
metadata$Group <- as.factor(metadata$Group)
metadata$Group = factor(metadata$Group, levels = c("2111","2201","2205","2207","2209","2211","2301","2303","2306"))
(p = alpha_boxplot1(alpha_div, metadata, "shannon","Group",facet_name = "shannon"))
metadata$Group = factor(metadata$Group, levels = c("2111","2201","2205","2207","2209","2211","2301","2303","2306"))
(p1 = alpha_boxplot1(alpha_div, metadata, "simpson", "Group",facet_name = "simpson"))
metadata$Group = factor(metadata$Group, levels = c("2111","2201","2205","2207","2209","2211","2301","2303","2306"))
(p2 = alpha_boxplot1(alpha_div, metadata, "Chao1", "Group",facet_name = "Chao1"))

ggsave(paste0("Gouqi shannon.svg"), p, width=140, height=90, units="mm")
ggsave(paste0("Gouqi simpson.svg"), p1, width=140, height=90, units="mm")

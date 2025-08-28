library(MetaCycle)
library(dplyr)
library(ggplot2)
library(tidyr)
library(rain)
#BiocManager::install("rain")

rm(list = ls())
setwd("/Volumes/xiaowang/zhoushan/新分析/diversity/16S18SV91_interaction/周期性检测")

pro16S <- read.table("sequences_rarefied_20240222.txt", header = TRUE, sep = "\t" , row.names = 1)
pro16S <-pro16S[,-c(1:7)]
phyto18S <- read.table("PR2_phyto_ASV.txt", header = TRUE, sep = "\t" , row.names = 1)
phyto18S <-phyto18S[,-c(1:9)]

total_abundance_16S <- rowSums(pro16S, na.rm = TRUE)
pro16S_relative <- pro16S / total_abundance_16S*100
pro16S_relative <- t(pro16S_relative)
#pro16S_relative$ASV_ID <-  rownames(pro16S_relative)
total_abundance_18S <- rowSums(phyto18S, na.rm = TRUE)
phyto18S_relative <- phyto18S / total_abundance_18S*100
phyto18S_relative <- t(phyto18S_relative)
#phyto18S_relative$ASV_ID <-  rownames(phyto18S_relative)
# 将 ASV_ID 列移动到第一列
#pro16S_relative <- pro16S_relative %>%
 # dplyr::select(ASV_ID, everything())
#phyto18S_relative <- phyto18S_relative %>%
 # dplyr::select(ASV_ID, everything())



#------------
?rain
phyto_year<-rain(phyto18S_relative,deltat=2,period=12,measure.sequence = c(4,4,0,6,6,6,6,6,6,6),adjp.method ="Bonferroni",verbose =TRUE )
#adjp.method ="Bonferroni",
pro_year<-rain(pro16S_relative,deltat=2,period=12,measure.sequence = c(4,4,0,6,6,6,6,6,6,6),adjp.method ="Bonferroni",verbose =TRUE)

write.csv(pro_year, file = "pro_year_results.csv", row.names = TRUE)
write.csv(phyto_year, file = "phyto_year_results.csv", row.names = TRUE)
#BH
#,method = "independent"
#period.delta=1,




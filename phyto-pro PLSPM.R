library(dplyr)
library(ggplot2)
library(vegan)
library(plspm)

rm(list = ls())

setwd("E:/Zhoushan eDNA/新分析/diversity/16S18SV91_interaction/PLSPM")

dat <- read.delim("Total.txt",header = T,row.names = 1)

dat_blocks <- list(
  Env = c('T',  'pH',"DO"), #
  Salinity = c("S"),
  Nurtrient = c('NO3.', 'NO2', 'PO4', 'SiO3', 'NH4'), 
  Biogeochemistry = c('DOC', 'POC', 'PON'), 
  Phyto =  c('phyto_Chao1','phyto_shannon','phyto_simpson',"phyto_MDS1"),
  Pro = c('pro_Chao1','pro_shannon','pro_simpson', "pro_MDS1")
)
#可以再加群落稳定性这些参数
dat_blocks
#通过 0-1 矩阵描述潜变量之间的关联，其中 0 代表变量间没有关联，1 代表有关联
Env <- c(0, 0, 0, 0, 0,0)
Salinity <- c(0, 0, 0, 0, 0,0)
Nurtrient <- c(0, 0, 0, 0, 0,0)
Biogeochemistry <- c(0, 0, 0, 0, 0,0)
Phyto <- c(1, 1,1, 1, 0, 0)
Pro <- c(1, 1, 1,1, 1,0 )



dat_path <- rbind(Env,Salinity ,Nurtrient, Biogeochemistry,  Phyto,Pro)
colnames(dat_path) <- rownames(dat_path)
dat_path
dat_modes <- rep('A', 6)
dat_modes
#?plspm
dat_pls <- plspm(dat, dat_path, dat_blocks, modes = dat_modes)
dat_pls
summary(dat_pls)

dat_pls$path_coefs
dat_pls$inner_model

innerplot(dat_pls, colpos = 'red', colneg = 'blue', show.values = TRUE, lcol = 'gray', box.lwd = 0)

dat_pls$inner_summary
dat_pls$gof





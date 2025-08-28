#弦图绘制＆需要写循环
library(reshape2)
library(circlize)
library(dplyr)
library(stringr)  # str_replace_all() 需要 stringr 包
library(Cairo)

rm(list = ls())
setwd("/Volumes/xiaowang/zhoushan/新分析/diversity/16S18SV91_interaction/network/PR2_0.6 0.05/弦图1/弦图")

dates <- c("2201")#"2111","2201","2205","2207","2209","2211","2301","2303"

for (date in dates) {
  pro_tax_filename <- paste0("pro_tax",date, ".txt")
  phyto_tax_filename <- paste0("phyto_tax",date, ".txt")
  edge_asv_filename <- paste0("edge_asv",date,".txt")

  png_filename <- paste0("Chord_Diagram_", date, ".svg")
  CairoSVG("Chord_Diagram.svg", width = 10, height = 10)
  svg(png_filename, width = 10, height = 10)


pro_tax <- read.table(pro_tax_filename, header = TRUE, sep = "\t")
pro_tax <- pro_tax %>% arrange(Class)
pro_group <- unique(pro_tax$Class)
pro_edge <- pro_tax$Pro


phyto_tax <- read.table(phyto_tax_filename, header = TRUE, sep = "\t")
phyto_tax <- phyto_tax %>% arrange(Phylum)
phyto_group <- unique(phyto_tax$Phylum)
phyto_edge <- phyto_tax$Phyto

#edge <- read.table("edge_asv.txt", header = TRUE, sep = "\t",row.names = 1)

edge <- read.table(edge_asv_filename, header = TRUE, sep = "\t",row.names = 1)

columns_to_select <- c("Dinophyceae", "Syndiniales", "Noctilucophyceae", 
                       "Mediophyceae", "OtherBacillariophyta", 
                       "OtherChlorophyta", "OtherCryptophyta", "OtherOchrophyta")

# 使用 dplyr 的 select() 函数选择存在的列
edge <- edge %>%
  select(all_of(intersect(columns_to_select, names(edge))))

# 初始化一个新的向量来存放最终的行名顺序
priority_rows <- c("SAR11clade","Rhodobacterales", "Puniceispirillales","OtherAlphaproteobacteria",  "Cellvibrionales","Betaproteobacteriales"
                   ,"Other Gammaproteobacteria","SAR86 clade","Flavobacteriales","Other Bacteroidia")
final_row_order <- vector()
# 检查每个优先行名是否存在于数据框中
for (row in priority_rows) {
  if (row %in% rownames(edge)) {
    final_row_order <- c(final_row_order, row)
  }
}
remaining_rows <- setdiff(rownames(edge), final_row_order)
final_row_order <- c(final_row_order, remaining_rows)
edge1 <- edge[final_row_order, ]
# 使用最终的行名顺序重新排序数据框
edge <- edge1


#kaishi
all_ID <- c(pro_edge, phyto_edge)

accum_pro <- rowSums(edge)

accum_phyto <- colSums(edge)

all_ID_xlim <- cbind(rep(0, length(all_ID)),data.frame(c(accum_pro, accum_phyto)))


#circlize 内圈连线数据

edge$otu_ID <- rownames(edge)

plot_data <- melt(edge, id = 'otu_ID') #此处使用了reshape2包中的melt()命令

colnames(plot_data)[2] <- 'Sample_ID'

plot_data <- plot_data[c(2, 1, 3, 3)]
print(plot_data)


#颜色设置


#color_otu <- c(
#  '#393b79', '#343b79', '#395b79', '#5654a3', '#5254a3', '#676ecf', '#686ecf', '#6b9ecf', '#9c8ede', '#9c9ede',
#  '#627939', '#637939', '#647939', '#8ca252', '#8ca252', '#b5cf6b', '#b5cf6b', '#b5cf6b'
#)
#'#de9ed6',
#'#8c6d31', '#8c6d31', '#8c6d31', '#bd9e39', '#bd9e39', '#e7ba52', '#e7ba52', '#e7cb94', '#e7cb94', '#e7cb94',
#'#843c39', '#843c39', '#ad424a', '#ad394a', '#ad494a', '#d6116b', '#d6616b', '#e7969c', '#e7919c', '#e7269c',
#'#7b4273', '#7b4173', '#a55194', '#a53194', '#a54194', '#ce61bd', '#ce6dbd', '#de1ed6', '#d4a017', '#d4b817',
#'#a8d017', '#b4e017', '#d4d917', '#e7c517', '#c59d17', '#ad7d17', '#b57617', '#d56117', '#e75117', '#cedb9c', '#cedb9c'

#49
#color_sample <- c(
#  '#3182bd', '#3182bd', '#6baed6', '#9ecae1', '#c6dbef', '#e6550d',
#  '#fd8d3c', '#fdae6b'
 # '#969696', '#bdbdbd', '#d9d9d9'
#)
#, '#fdd0a2', '#31a354', '#74c476', '#a1d99b',
#'#c7e9c0', '#756bb1', '#9e9ac8', '#bcbddc', '#dadaeb', '#636363'#,
#21
#color_phylum <-  c(
#  '#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a',
#  '#d62728', '#ff9896', '#9467bd'
  #'#17becf', '#c5b0d5', '#8c564b', '#c49c94',
#  '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d'#,
#)

#19
#color_group <- c("#fdd0a2", "#fd10a2", "#3182bd", "grey", '#d56117')#, '#e75117', '#bcbd22', '#dbdb8d',"grey", '#d56117'
#3

# 确保颜色向量长度匹配
#names(color_otu) <- pro_edge
#names(color_sample) <- phyto_edge


gap_size <- c(rep(3, length(pro_edge) - 1), 6, rep(3, length(phyto_edge) - 1), 6)

circos.par(cell.padding = c(0, 0, 0, 0), start.degree = 270, gap.degree = gap_size)

circos.initialize(factors = factor(all_ID, levels = all_ID), xlim = all_ID_xlim)


circos.trackPlotRegion(
  ylim = c(0, 1), track.height = 0.1, bg.border = NA, #track.height = 0.03
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data('sector.index')
    xlim = get.cell.meta.data('xlim')
    ylim = get.cell.meta.data('ylim')
  } )

#--------------设定颜色-------------------
predefined_colors <- c(
  Alphaproteobacteria = "#AEB8FE",
  Bacteroidia = "#87B6A7",
  Deltaproteobacteria = "#5654a3",
  Gammaproteobacteria = "#FFE066",
  Betaproteobacteriales = "#676ecf",
  Nitrososphaeria = "#686ecf",
  Nitrospinia = "#6b9ecf",
  Planctomycetacia = "#F4C3C2",
  Thermoplasmata = "#9c9ede",
  OhtersProkaryotes =  "#999999",
  Bacillariophyta = "#32908F",
  Chlorophyta =  "#88D498",
  Dinophyta =  "#9D79BC",
  OtherOchrophyta = "#999999",
  #----------then for the  circle[2]-----------
  Rhodobacterales = "#E6FDFF",
  SAR11clade = "#758BFD",
  Flavobacteriales =  "#FF8600",
  Nitrosopumilales = "#317B22" ,
  MarineGroupII = "#67E0A3" ,
  OtherAlphaproteobacteria = "#999999",
  OtherGammaproteobacteria = "#999999",
  OtherDeltaproteobacteria = "#999999",
  Betaproteobacteriales ="#3993DD" ,
  Mediophyceae = "#80CED7",
  OtherBacillariophyta = "#114B5F",
  OtherChlorophyta = "#999999",
  OtherCryptophyta = "#999999",
  Dinophyceae = "#D4C2FC",
  Noctilucophyceae = "#772D8B",
  Syndiniales = "#998FC7",
  OtherOchrophyta = "#999999",
  Default = "#999999"  # 默认颜色
)
#--------------设定颜色-------------------
for (i in 1:length(pro_group)) {
  
  group_name <- pro_group[i]
  
  group_color <- ifelse(!is.na(predefined_colors[group_name]), predefined_colors[group_name], predefined_colors["Default"])
  
  tax_OTU <- {subset(pro_tax, Class == pro_group[i])}$Pro
  
  highlight.sector(tax_OTU, track.index = 1, col = group_color, text = pro_group[i], cex = 1, text.col = 'black', niceFacing = FALSE)#text = pro_group[i]
  #col = color_phylum[i]#cex = 0.7
}



for (i in 1:length(phyto_group)) {
  
  group_name <- phyto_group[i]
  
  phyto_color <- ifelse(!is.na(predefined_colors[group_name]), predefined_colors[group_name], predefined_colors["Default"])
  
  group_sample <- {subset(phyto_tax, Phylum == phyto_group[i])}$Phyto
  
  highlight.sector(group_sample, track.index = 1, col = phyto_color, text =phyto_group[i] , cex = 1, text.col = 'black', niceFacing = FALSE)#text = phyto_group[i]
  #col = color_group[i]
}




#添加 OTU 百分比注释（第二圈）

#circos.trackPlotRegion(
#  
#  ylim = c(0, 1), track.height = 0.05, bg.border = NA,
  
#  panel.fun = function(x, y) {
    
#    sector.index = get.cell.meta.data('sector.index')
    
 #   xlim = get.cell.meta.data('xlim')
    
  #  ylim = get.cell.meta.data('ylim')
    
#  } )



#circos.track(
  
  #track.index = 2, bg.border = NA,
  
 # panel.fun = function(x, y) {
    
    #xlim = get.cell.meta.data('xlim')
    
   # ylim = get.cell.meta.data('ylim')
    
   # sector.name = get.cell.meta.data('sector.index')
    
   # xplot = get.cell.meta.data('xplot')
    
    
    
   # by = ifelse(abs(xplot[2] - xplot[1]) > 30, 0.25, 1)
    
 #   for (p in c(0, seq(by, 1, by = by))) circos.text(p*(xlim[2] - xlim[1]) + xlim[1], mean(ylim) + 0.4, paste0(p*100, '%'), cex = 0.4, adj = c(0.5, 0), niceFacing = FALSE)
    
    
    
 #   circos.lines(xlim, c(mean(ylim), mean(ylim)), lty = 3)
    
 # } )



#绘制 OTU、样本主区块（第三圈）

bg_colors <- sapply(all_ID, function(x) {
  if (x %in% names(predefined_colors)) {
    predefined_colors[x]
  } else {
    predefined_colors["Default"]  # 使用默认颜色，如果没有匹配
  }
})

circos.trackPlotRegion(
  
  ylim = c(0, 1), track.height = 0.05, bg.col = bg_colors, bg.border = NA, track.margin = c(0, 0.01),#bg.col = c(color_otu, color_sample)
  
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data('xlim')
    
    sector.name = get.cell.meta.data('sector.index')
    
    #circos.axis(h = 'top', labels.cex = 0.4, major.tick.percentage = 0.4, labels.niceFacing = FALSE)

  #  circos.text(mean(xlim), 0.2, sector.name, cex = 0.4, niceFacing = FALSE, adj = c(0.5, 0))
    
  } )



#绘制 OTU、样本副区块（第四圈）

#circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.03, track.margin = c(0, 0.01))

#绘制 OTU-样本关联连线（最内圈）
print(accum_pro)
print(accum_phyto)

#绘制 OTU-样本关联连线（最内圈）

add_alpha_to_color <- function(color, alpha = '70') {
  # 将十六进制颜色转换为透明颜色
  paste0(color, alpha)
}

for (i in seq_len(nrow(plot_data))) {
  
  otu_key <- plot_data[i, 1]  # 提取OTU键值
  otu_key_char <- as.character(otu_key)

  base_color <- ifelse(otu_key %in% names(predefined_colors),
                       predefined_colors[otu_key_char],
                       predefined_colors["Default"])
  final_color <- add_alpha_to_color(base_color)  # 添加默认透明度
  circos.link(
    
    plot_data[i,2], c(accum_pro[plot_data[i,2]], accum_pro[plot_data[i,2]] - plot_data[i,4]),
    
    plot_data[i,1], c(accum_phyto[plot_data[i,1]], accum_phyto[plot_data[i,1]] - plot_data[i,3]),
    
    col = final_color, border = NA)
  
#    col = paste0(color_otu[plot_data[i,2]], '70'), border = NA )
  
 # circos.rect(accum_pro[plot_data[i,2]], 0, accum_pro[plot_data[i,2]] - plot_data[i,4], 1, sector.index = plot_data[i,2], col = color_sample[plot_data[i,1]], border = NA)
  
  #circos.rect(accum_phyto[plot_data[i,1]], 0, accum_phyto[plot_data[i,1]] - plot_data[i,3], 1, sector.index = plot_data[i,1], col = color_otu[plot_data[i,2]], border = NA)
  
  
  
  accum_pro[plot_data[i,2]] = accum_pro[plot_data[i,2]] - plot_data[i,4]
  
  accum_phyto[plot_data[i,1]] = accum_phyto[plot_data[i,1]] - plot_data[i,3]
  
}

dev.off()  #
}

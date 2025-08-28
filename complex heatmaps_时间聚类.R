
rm(list = ls())
#setwd("/yourwordpath")

library(dplyr)
library(ggplot2)
library(stringr)
 
df <- read.table("ASVtable.txt", header = TRUE, sep = "\t",fill = TRUE)

#计算相对丰度
char_cols <- df %>% select(1:8)
relative_abundance <- df %>% 
  select(9:ncol(df)) %>%
  mutate(across(everything(), ~ . / sum(.) * 100))
result <- bind_cols(char_cols, relative_abundance)

# 对Order进行分组，并计算每个样品中的总数，先计算相对丰度，再合并
Order_sums <- result %>%
  group_by(Phylum,Order) %>%
  summarise(across(.cols = 7:56, .fns = sum), .groups = 'drop')

Order_sums <- Order_sums %>%
  mutate(Phylum = str_replace_all(Phylum, "D_1__", ""))

Order_sums <- Order_sums %>%
  mutate(Order = str_replace_all(Order, "D_3__", ""))

#这边就需要去除""空白内容了
result_cleaned <- Order_sums %>%
  filter(across(everything(), ~ . != ""))
#筛选>1%的Orders(value>0.01),if 一行有一个值大于0.01，则保留该行，否则删除该行
filtered_Orders <- result_cleaned %>%
  filter(if_any(3:52, ~ . > 1))
#这边站位需要排序 1 3 5 7 9 10？还是以9 7 5 1 3 10 
new_Order <- c("Phylum", "Order", "S2111_1", "S2111_5", "S2111_9", "S2111_10", "S2201_1", "S2201_5", "S2201_9", "S2201_10", 
               "S2205_1", "S2205_3","S2205_5", "S2205_7","S2205_9", "S2205_10", "S2207_1", "S2207_3","S2207_5", "S2207_7","S2207_9", "S2207_10", 
               "S2209_1", "S2209_3","S2209_5", "S2209_7","S2209_9", "S2209_10","S2211_1", "S2211_3","S2211_5", "S2211_7","S2211_9", "S2211_10",
               "S2301_1", "S2301_3","S2301_5", "S2301_7","S2301_9", "S2301_10","S2303_1", "S2303_3","S2303_5", "S2303_7","S2303_9", "S2303_10",
               "S2306_1", "S2306_3","S2306_5", "S2306_7","S2306_9", "S2306_10"
               )#……
filtered_Orders_reOrdered <- filtered_Orders %>%
  select(all_of(new_Order), everything())

filtered_Orders_reOrdered <- filtered_Orders_reOrdered %>%
  mutate(across(-c(Phylum, Order), as.numeric))

filtered_Orders_reOrdered <- filtered_Orders_reOrdered %>%
  mutate(across(-c(Phylum, Order), ~ ifelse(is.na(.), 0, .)))

#似乎就可以进行热图的制作了
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
row_names <- filtered_Orders_reOrdered$Order

# 准备矩阵数据，并移除 Class, Order, SampleType 列
data_matrix <- filtered_Orders_reOrdered %>%
  select(-Phylum, -Order) %>%
  as.matrix()

# 提取行名（Order）
rownames(data_matrix) <- row_names

row_groups <- filtered_Orders_reOrdered$Phylum

# 绘制热图


Phylum_colors <- c(
  "Phylum1" = "#1f77b4", "Phylum2" = "#ff7f0e", "Phylum3" = "#2ca02c", "Phylum4" = "#d62728", 
  "Phylum5" = "#9467bd", "Phylum6" = "#8c564b", "Phylum7" = "#e377c2", "Phylum8" = "#7f7f7f", 
  "Phylum9" = "#bcbd22", "Phylum10" = "#17becf", "Phylum11" = "#aec7e8", "Phylum12" = "#ffbb78",
  "Phylum13" = "#98df8a"
)#, "Phylum14" = "#ff9896", "Phylum15" = "#c5b0d5", "Phylum16" = "#c49c94"
names(Phylum_colors) <- unique(row_groups)

# 宏转录组采样信息通过metadata添加？季节，采样月份的颜色
metadata <- read.table("sample-metadata.txt", header = TRUE, sep = "\t",fill = TRUE)

top_annotation_data <- metadata %>%
  filter(sample.id %in% colnames(data_matrix))

# 创建颜色映射
#group_colors <- c("Spring" = "#2A9134", "Summer" = "#FF521B", "Autumn" = "#FED766", "Winter" = "#B5FFE9")
group_colors <- c("Nov 2021" = "#2A9134", "Jan 2022" = "#FEE440", "May 2022" = "#6247AA", "Jul 2022" = "#B5FFE9",
                  "Sep 2022" = "#F9A620", "Nov 2022" = "#D1B1C8",
                  "Jan 2023" = "#57C4E5","Mar 2023" = "#D1D646","Jun 2023" = "#F97068")
metatranscript_colors <- c("no" = "grey", "yes" = "green")
#top_annotation_data$Season <- factor(top_annotation_data$Season, levels = c("Spring", "Summer", "Autumn", "Winter"))
top_annotation_data$Group <- factor(top_annotation_data$Group, levels = c("Nov 2021", "Jan 2022", "May 2022", "Jul 2022",
                                                                          "Sep 2022" ,"Nov 2022","Jan 2023","Mar 2023","Jun 2023" ))
#注释内容
ha <- HeatmapAnnotation(
  #Season = top_annotation_data$Season,
  Month = top_annotation_data$Group,
  #Metatranscript = top_annotation_data$Metatranscrpt,
 # col = list(Season = group_colors, Metatranscript = metatranscript_colors),
 col = list(Month = group_colors),
  simple_anno_size = unit(3, "mm"),
  annotation_name_gp = gpar(fontsize = 6) 
  )
 #annotation_height =unit(c(1, 1), c("cm", "null")), height = unit(5, "cm") )


ht <- Heatmap(data_matrix,
        name = "Relative Abundance", # 图例标题
        row_title = "",
        column_title = "",
        show_column_names = FALSE, # 显示列名
        cluster_columns = FALSE, # 不进行列聚类
        row_dend_width = unit(15, "mm"),
        row_names_gp = gpar(fontsize = 6) ,#, fontfamily = "Times"
        width = unit(ncol(data_matrix) * 1, "mm"), # 热图宽度
        height = unit(nrow(data_matrix) * 2, "mm"), # 热图高度
        left_annotation = rowAnnotation(
          Phylum = row_groups,
          col = list(Phylum = Phylum_colors),
          annotation_name_gp = gpar(fontsize = 6),
          width = unit(6, "mm") # 设置left annotation的宽度
        ),
        top_annotation = ha ,
        col = colorRamp2(c(0, 1,70), # 10,20,30,40,50,60,
                         c("white", "yellow", "red")), # 颜色梯度# "lightblue", "lightgreen", "lightyellow", "orange", "lightcoral", "coral",
        heatmap_legend_param = list(title = "Relative Abundance", at = c(0, 1, 10, 20, 30, 40, 50, 60, 70), 
                                    labels = c("0", "1%", "10%", "20%", "30%", "40%", "50%", "60%", "70%"),
                                    title_position = "topcenter",  
                                    title_gp = gpar(fontsize = 5),  # 设置标题字体大小
                                    labels_gp = gpar(fontsize = 5), # 设置标签字体大小
                                    legend_height = unit(8, "cm"),  # 设置图例框的高度
                                    legend_width = unit(4, "cm")  ,  # 设置图例框的宽度
                                    direction = "horizontal",   
                                    gap = unit(1000, "mm")) # 图例设置
)




draw(ht, heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom")
#单独画图例


pdf("16S_heatmap.pdf", width = 5, height = 8)
draw(ht, heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom")
dev.off()


# 添加列的分组

annot_df <- data.frame(cyl = mtcars$cyl, am = mtcars$am,
                       mpg = mtcars$mpg)
# row.names(annot_df) = row.names(mtcars)
# Define colors for each levels of qualitative variables
# Define gradient color for continuous variable (mpg)
col = list(cyl = c("4" = "green", "6" = "gray", "8" = "darkred"),
           am = c("0" = "yellow", "1" = "orange"),
           mpg = circlize::colorRamp2(c(17, 25),
                                      c("lightblue", "purple")) )
# Create the heatmap annotation
ha <- HeatmapAnnotation(df = data.frame(cyl = mtcars$cyl, am = mtcars$am,
                                        mpg = mtcars$mpg), col = col)

# Combine the heatmap and the annotation
# df = t(df)
Heatmap(df, name = "mtcars",
        top_annotation = ha)


# Annotation data frame

annot_df <- data.frame(cyl = mtcars$cyl, am = mtcars$am,
                       mpg = mtcars$mpg)
# row.names(annot_df) = row.names(mtcars)
# Define colors for each levels of qualitative variables
# Define gradient color for continuous variable (mpg)
col = list(cyl = c("4" = "green", "6" = "gray", "8" = "darkred"),
           am = c("0" = "yellow", "1" = "orange"),
           mpg = circlize::colorRamp2(c(17, 25),
                                      c("lightblue", "purple")) )
# Create the heatmap annotation
ha <- HeatmapAnnotation(df = data.frame(cyl = mtcars$cyl, am = mtcars$am,
                                        mpg = mtcars$mpg), col = col)

# Combine the heatmap and the annotation
df = t(df)
Heatmap(df, name = "mtcars") +
  Heatmap(mtcars$mpg, name = "type", width = unit(5, "mm")) +
  Heatmap(mtcars$mpg, name = "type", width = unit(5, "mm"))

Heatmap(df, name = "mtcars") +
  Heatmap(annot_df , name = "type", width = unit(5, "mm"))







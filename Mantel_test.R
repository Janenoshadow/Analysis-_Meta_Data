rm(list = ls())
#setwd("")

library(dplyr)
library(ggplot2)
library(vegan)
library(linkET)
library(FD)

env <- read.table("eenv.txt",header = T,row.names = 1)#environment factors2205.txt
asv <- read.table("asvtable.txt",header = T,row.names = 1)

mantel <- mantel_test(asv, env,
                      spec_select = list(Alphaproteobacteria = 1:880,#your design object
                                         Bacteroidia = 881:1553,
                                         Gammaproteobacteria= 1554:2679,
                                         Nitrososphaeria= 2680:2768),#,                      Oxyphotobacteria=2768:2819
                      na_omit = TRUE) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
#> `mantel_test()` using 'bray' dist method for 'spec'.
#> `mantel_test()` using 'euclidean' dist method for 'env'.

my_plot <- qcorrplot(correlate(env), type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(aes(colour = pd, size = rd), 
              data = mantel, 
              curvature = nice_curvature()) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  scale_size_manual(values = c(0.5, 1, 3)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))
my_plot + theme(legend.text = element_text(face = "bold",size = 50))

my_plot
# 保存为PNG
ggsave("my_plot.png", plot = my_plot, width = 20, height = 14)

# 保存为PDF
ggsave("my_plot.svg", plot = my_plot, width = 20, height = 14)

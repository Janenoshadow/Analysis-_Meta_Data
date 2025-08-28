# contributor: Daliang Ning
rm(list = ls())
library(tidyr)
source("/ieggrtools.r")
#setwd("/yourworkpath")

comm1 <- read.table("sequences_rarefied_20240222.txt", header = TRUE, sep = "\t", fill = TRUE,row.names = 1)
comm_reduced <- comm1[,-(1:7)]
comm_transposed <- t(comm_reduced)
comm <- as.data.frame(comm_transposed)
treat <-read.table("sample-metadata_modified.txt", header = TRUE, sep = "\t", fill = TRUE,row.names = 1)#去除前8行以匹配数据

#####################

sum(is.na(comm)) # check NA
comm[is.na(comm)]=0# if have, should change to zero

sampc=match.name(rn.list = list(comm=comm,treat=treat))#整合为一个大数据狂
comm=sampc$comm
treat=sampc$treat

head(treat) # check which column is plot ID, which column is Group.
plot.lev=unique(treat$station)
Group.lev=sort(unique(treat$Group))

zeta.lev=2:length(Group.lev)

Group.windows=lapply(1:length(zeta.lev),
                  function(i)
                  {
                    zetai=zeta.lev[i]
                    lapply(1:(length(Group.lev)-zetai+1),function(j){Group.lev[j:(j+zetai-1)]})
                  })
names(Group.windows)=zeta.lev
Group.windows

# function of community stability among a group of samples
comstab<-function(subcom){((nrow(subcom)*sum(apply(subcom,2,min)))/sum(subcom))^0.5}

#######

stabl=lapply(1:length(Group.windows),
             function(i)
             {
               stabi=t(sapply(1:length(plot.lev),
                              function(j)
                              {
                                plotj=plot.lev[j]
                                sapply(1:length(Group.windows[[i]]),
                                       function(k)
                                       {
                                         Groupwdk=Group.windows[[i]][[k]]
                                         sampijk=rownames(treat)[which((treat$station==plotj) & (treat$Group %in% Groupwdk))]
                                         outijk=NA
                                         if(length(sampijk) < length(Groupwdk))
                                         {
                                           warning("plot ",plotj," has missing Group in Group window ",paste(Groupwdk,collapse = ","))
                                         }else if(length(sampijk) > length(Groupwdk)){
                                           warning("plot ",plotj," has duplicate samples in at least one Group of window ",paste(Groupwdk,collapse = ","))
                                         }else{
                                           comijk=comm[which(rownames(comm) %in% sampijk),,drop=FALSE]
                                           outijk=comstab(comijk)
                                         }
                                         outijk
                                       })
                              }))
               if(nrow(stabi)!=length(plot.lev) & nrow(stabi)==1){stabi=t(stabi)}
               rownames(stabi)=plot.lev
               colnames(stabi)=sapply(Group.windows[[i]],function(v){paste0("Zeta",zeta.lev[i],paste0(v,collapse = ""))})
               stabi
             }) #
stabm=Reduce(cbind,stabl)
#对于我自己，提取前八列 导出 进行图例的制作就可以了？
stabm_df <- as.data.frame(stabm)
stabm_df_reduced <- stabm_df[, 1:8]
stabm_df_reduced$station <- rownames(stabm_df_reduced)
write.csv(stabm_df_reduced, "stabm_reduced.csv", row.names = TRUE)
write.csv(stabm_df, "stabm_full.csv", row.names = TRUE)

library(mgcv)
library(ggplot2)


stabm_df_reduced <- read.table("stabm_reduced_gam.txt", header = TRUE, sep = "\t", fill = TRUE,row.names = 1)

stabm_long <- gather(stabm_df_reduced, key = "Time", value = "Stability", -station)

# 将时间点转换为连续变量，这里假设Zeta2M2111M2201表示一个时间点
# 例如，可以用正则表达式提取时间点，这里需要具体的时间点格式，以下仅为示例
stabm_long$Time <- gsub("X", "", stabm_long$Time)

library(lubridate)

library(dplyr)

# 假设stabm_long已经通过gather或pivot_longer创建，并且Time列包含原始的时间戳格式"XYYYY.MM.DD"

# 去除前缀'X'并转换时间格式为"YYYY-MM"
stabm_long <- stabm_long %>%
  mutate(
    TimeFormatted = gsub("^X", "", Time),  # 去除'X'
    YearMonth = format(as.Date(TimeFormatted, format="%Y.%m.%d"), "%Y-%m")  # 转换并格式化日期
  )
stabm_long$Date <- as.Date(paste0(stabm_long$YearMonth, "-01"))

# 提取年份和月份
stabm_long$Year <- as.numeric(format(stabm_long$Date, "%Y"))
stabm_long$Month <- as.numeric(format(stabm_long$Date, "%m"))

# 假设开始年份和月份
start_year <- 2022
start_month <- 1

# 计算从开始日期到每个观测日期的相对月份
stabm_long$MonthsSinceStart <- ((stabm_long$Year - start_year) * 12) + (stabm_long$Month - start_month)

gam_model <- gam(Stability ~ s(MonthsSinceStart, k = 5), data = stabm_long, method = "REML")

# 对每个群落类型进行GAM拟合
#gam_model <- gam(Stability ~ s(Time), data=stabm_long, method="REML") 我的数据不足以支持


summary(gam_model)

gam.check(gam_model)

#cv.gam(gam_model)




# 创建预测数据框，使用正确的时间变量
pred_data <- expand.grid(
  MonthsSinceStart = seq(min(stabm_long$MonthsSinceStart), max(stabm_long$MonthsSinceStart), length.out = 100)
)

# 进行预测，并添加置信区间
pred_data$Stability <- predict(gam_model, newdata=pred_data, type="response")
confidence_intervals <- predict(gam_model, newdata=pred_data, type="response", se=TRUE)
pred_data$UpperCI <- pred_data$Stability + 1.96 * confidence_intervals$se.fit
pred_data$LowerCI <- pred_data$Stability - 1.96 * confidence_intervals$se.fit

# 绘图
#p <- ggplot(stabm_long, aes(x = MonthsSinceStart, y = Stability)) +
#  geom_ribbon(data=pred_data, aes(x=MonthsSinceStart, ymin=LowerCI, ymax=UpperCI), alpha=0.2, fill="blue") +
#  geom_line(data=pred_data, aes(x=MonthsSinceStart, y=Stability), color="blue") +
#  geom_point(aes(x=MonthsSinceStart, y=Stability)) +
#  theme_minimal() +
#  labs(x = "Months Since Start", y = "Community Stability", title = "GAM fit for Community Stability Over Time")
#p
monthly_stats <- stabm_long %>%
  group_by(MonthsSinceStart) %>%
  summarize(mean_stability = mean(Stability, na.rm = TRUE),
            se = sd(Stability, na.rm = TRUE) / sqrt(n()))
library(ggthemes)


smooth_term_p_value <- summary(gam_model)$p.values[1] # Replace with actual index or extraction method

# Create the annotation text

explained_variance <- "45.4%"  # Placeholder for the explained variance
significance_level <- "p < 0.001"  

# 绘制每个月份的点和误差线


p2 <- ggplot() +
  geom_point(data=monthly_stats, aes(x=MonthsSinceStart, y=mean_stability)) +
  geom_errorbar(data=monthly_stats, aes(x=MonthsSinceStart, ymin=mean_stability-se, ymax=mean_stability+se), width=0.5) +
  geom_line(data=pred_data, aes(x=MonthsSinceStart, y=Stability), color="#2B2D42", size = 1) +
  geom_ribbon(data=pred_data, aes(x=MonthsSinceStart, ymin=LowerCI, ymax=UpperCI), alpha=0.2, fill="#92DCE5") +
  theme_base() +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1, color = "black"),  # x 轴字体大小和样式
    axis.text.y = element_text(size = 6, color = "black"),  # y 轴字体大小和样式
    axis.title.x = element_text(size = 6, face = "bold"),  # x 轴标题字体大小
    axis.title.y = element_text(size = 6, face = "bold")   # y 轴标题字体大小
  ) +
  labs(x = "", y = "Stability", title = NULL) +
  scale_y_continuous(
    limits = c(0.15, 0.83),  # 设置 y 轴的范围，比如 0 到 1
    breaks = seq(0, 1, by = 0.2)
  )+
  scale_x_continuous(breaks = c(0, 4, 6, 8, 10,12,14,17),  # 指定所需的x轴刻度
                     labels = c("Nov-Jan", "Jan-May", "May-Jul", "Jul-Sep", "Sep-Nov", "Nov-Jan", "Jan-Mar", "Mar-Jun")) +  # 对应的标签
  annotate("text", x = Inf, y = Inf, label = paste("Explained variance:", explained_variance, "\nSignificance:", significance_level),
           hjust = 1, vjust = 1, size = 4, color = "#2B2D42", position = position_nudge(x = -1.5, y = -1.5))
p2
ggsave("ok.svg",plot = p2,width = 4,height = 3.5)
ggsave("ok.png",plot = p2,width = 4,height = 3.5)





# 导出预测数据及置信区间
write.csv(pred_data, "predicted_stability_with_CIs.csv", row.names = FALSE)
# 导出处理后的原始数据
write.csv(stabm_long, "processed_original_data.csv", row.names = FALSE)
# 捕获模型摘要的输出并写入文件
capture.output(summary(gam_model), file = "gam_model_summary.txt")
# 捕获gam.check的输出并写入文件
capture.output(gam.check(gam_model), file = "gam_model_check.txt")











treat.rm09=treat[which(treat$Group!="Y09"),,drop=FALSE] # just to annotate the treatment information of each plot, remove the special Group
output=data.frame(treat.rm09[match(rownames(stabm),treat.rm09$station),1:5,drop=FALSE],stabm,stringsAsFactors = FALSE)
output=output[order(output$Warming,output$Clipping,output$Precip,output$station),,drop=FALSE]
head(output)
#save.file(output,filename = "MultiOrderStability.csv")

# plot Figure S7abcde - multiorder compositional stability. Adapted from the concept of Zeta diveristy in Hui and McGeoch American Naturalist 2014.
cs.mo = output

cs.long <- gather(cs.mo, Group, cs, Zeta2Y09Y10:Zeta6Y09Y10Y11Y12Y13Y14, factor_key=TRUE)
cs.long <- cs.long[-which(is.na(cs.long$cs)),]
cs.long$EndGroup = as.numeric(gsub(".*Y", "", cs.long$Group))

get_text <- function(y, x){
  lm = lm(y~x)
  lm_p = round(summary(lm)$coefficients[2,4],3)
  lm_r = round(summary(lm)$adj.r.squared,3)
  lm_s = round(summary(lm)$coefficients[2,1],3)
  txt = paste(lm_s, " (", lm_r, ", ", lm_p, ")", sep="")
  return(txt)
}

plot_cs_full <- function(xc, yc, xw, yw){
  #  yl=c(min(yc, yw), max(yc, yw))

  plot(yw~xw, ylim=c(0.2,1), xlim=c(10,14), ylab="Compositional stability", xlab="Group", pch=19, cex.axis=0.7, cex.lab=0.7, tck=-0.03, col="#e7211f")
  abline(lm(yw~xw), col="#e7211f")

  points(yc~xc, col="#214da0")
  abline(lm(yc~xc), col="#214da0")

  txt1 = get_text(yw,xw)
  txt2 = get_text(yc,xc)
  mtext(txt1, side=3, line=0.8, cex=0.5, col="#e7211f")
  mtext(txt2, side=3, line=0.1, cex=0.5, col="#214da0")
}

plot_cs_zeta6 <- function(xc, yc, xw, yw){
  #  yl=c(min(yc, yw), max(yc, yw))

  plot(yw~xw, ylim=c(0.2,1), xlim=c(10,14), ylab="Compositional stability", xlab="Group", pch=19, cex.axis=0.7, cex.lab=0.7, tck=-0.03, col="#e7211f")
  points(yc~xc, col="#214da0")

  wt = wilcox.test(yw, yc)
  mtext(paste(wt$statistic, wt$p.value, sep=","), cex=0.5)
}

par(mfrow=c(2,3))
plot_cs_full(xc=cs.long$EndGroup[intersect(which(cs.long$Warming=="N"), grep("Zeta2", cs.long$Group))],
             yc=cs.long$cs[intersect(which(cs.long$Warming=="N"), grep("Zeta2", cs.long$Group))],
             xw=cs.long$EndGroup[intersect(which(cs.long$Warming=="W"), grep("Zeta2", cs.long$Group))],
             yw=cs.long$cs[intersect(which(cs.long$Warming=="W"), grep("Zeta2", cs.long$Group))])

plot_cs_full(xc=cs.long$EndGroup[intersect(which(cs.long$Warming=="N"), grep("Zeta3", cs.long$Group))],
             yc=cs.long$cs[intersect(which(cs.long$Warming=="N"), grep("Zeta3", cs.long$Group))],
             xw=cs.long$EndGroup[intersect(which(cs.long$Warming=="W"), grep("Zeta3", cs.long$Group))],
             yw=cs.long$cs[intersect(which(cs.long$Warming=="W"), grep("Zeta3", cs.long$Group))])

plot_cs_full(xc=cs.long$EndGroup[intersect(which(cs.long$Warming=="N"), grep("Zeta4", cs.long$Group))],
             yc=cs.long$cs[intersect(which(cs.long$Warming=="N"), grep("Zeta4", cs.long$Group))],
             xw=cs.long$EndGroup[intersect(which(cs.long$Warming=="W"), grep("Zeta4", cs.long$Group))],
             yw=cs.long$cs[intersect(which(cs.long$Warming=="W"), grep("Zeta4", cs.long$Group))])

plot_cs_full(xc=cs.long$EndGroup[intersect(which(cs.long$Warming=="N"), grep("Zeta5", cs.long$Group))],
             yc=cs.long$cs[intersect(which(cs.long$Warming=="N"), grep("Zeta5", cs.long$Group))],
             xw=cs.long$EndGroup[intersect(which(cs.long$Warming=="W"), grep("Zeta5", cs.long$Group))],
             yw=cs.long$cs[intersect(which(cs.long$Warming=="W"), grep("Zeta5", cs.long$Group))])

plot_cs_zeta6(xc=cs.long$EndGroup[intersect(which(cs.long$Warming=="N"), grep("Zeta6", cs.long$Group))],
              yc=cs.long$cs[intersect(which(cs.long$Warming=="N"), grep("Zeta6", cs.long$Group))],
              xw=cs.long$EndGroup[intersect(which(cs.long$Warming=="W"), grep("Zeta6", cs.long$Group))],
              yw=cs.long$cs[intersect(which(cs.long$Warming=="W"), grep("Zeta6", cs.long$Group))])


# plot Figure 3d. Compositional stability of adjacent two Groups

colors = c(rep("#214da0",5), rep("#e7211f", 5))
syms = c(1,1,1,1,1, 19,19,19,19,19)

plot_cs <- function(yw, yc, ebw, ebc){
  x=c(1,2,3,4,5)
  yl=c(min(c(yw-ebw,yc-ebc)), max(c(yw+ebw,yc+ebc)))

  plot(yw~x, ylim=yl, ylab="Compositional stability", xlab="Group of warming", pch=19, cex.axis=0.7, cex.lab=0.7, tck=-0.03, col="red")
  abline(lm(yw~x), col="red")
  arrows(x, yw-ebw, x, yw+ebw, length=0.05, angle=90, code=3, col="red")

  points(yc~x, col="blue")
  abline(lm(yc~x), lty=2, col="blue")
  arrows(x, yc-ebc, x, yc+ebc, length=0.05, angle=90, code=3, col="blue")

  txt1 = get_text(yw,x)
  txt2 = get_text(yc,x)
  mtext(txt1, side=3, line=0.8, cex=0.5, col="red")
  mtext(txt2, side=3, line=0.1, cex=0.5, col="blue")
}

zeta2 = cs.long[grep("Zeta2", cs.long$Group),]
zeta2_means = aggregate(zeta2$cs, list(zeta2$Warming, zeta2$EndGroup), FUN=mean)
zeta2_sds = aggregate(zeta2$cs, list(zeta2$Warming, zeta2$EndGroup), FUN=sd)

plot_cs(yw=zeta2_means$x[which(zeta2_means$Group.1=="W")],
        yc=zeta2_means$x[which(zeta2_means$Group.1=="N")],
        ebw=zeta2_sds$x[which(zeta2_sds$Group.1=="W")],
        ebc=zeta2_sds$x[which(zeta2_sds$Group.1=="N")])

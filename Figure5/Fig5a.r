library("ggplot2")
library("dplyr")
setwd("/home/yiliao/Documents/Pepper_2021/Final/6_Figure_loop/FINAL")
#####Loop size
mydata1 <-read.table("All.bed.bed.bed",head=T)
head(mydata1)
mydata1$Resolution<- as.factor(mydata1$Resolution)
mydata1$Tissues<- as.factor(mydata1$Tissues)
ggplot(mydata1, aes(x=Resolution, y=Size, fill=factor(Tissues))) + geom_boxplot(outlier.shape = NA) + stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red") + ylim(0,8000000) +labs(fill = "Tissues")

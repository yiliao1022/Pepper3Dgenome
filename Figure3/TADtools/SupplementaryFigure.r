library("ggplot2")
library("gplots")
library("amap")
library("RColorBrewer")

### For tadtool 
## For num
setwd("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_1/TadTool")
data <- read.table("Tadtool.num.mat", sep= "\t", header=FALSE)
data
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
Dis <- data.matrix(data[1:8])
Distance <-1-Dis
Distance
Dist<-as.dist(Distance)
Dist
hc <- hclust(Dist)                # apply hirarchical clustering 
plot(hc)
Dis
heatmap.2(Dis, Rowv=as.dendrogram(hc), symm=T, trace="none",
          col=hmcol, margins=c(4,4), main="")
dev.copy2pdf(file="/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_1/TadTool/TadTool.num.pdf")


## For boundaries
setwd("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_1/TadTool")
data <- read.table("Tadtool.boundaries.mat", sep= "\t", header=FALSE)
data
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
Dis <- data.matrix(data[1:8])
Distance <-1-Dis
Distance
Dist<-as.dist(Distance)
Dist
hc <- hclust(Dist)                # apply hirarchical clustering 
plot(hc)
Dis
heatmap.2(Dis, Rowv=as.dendrogram(hc), symm=T, trace="none",
          col=hmcol, margins=c(4,4), main="")
dev.copy2pdf(file="/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_1/TadTool/TadTool.boundaries.pdf")








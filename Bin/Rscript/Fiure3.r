library("ggplot2")
library("gplots")
library("amap")
library("RColorBrewer")

###For Jaccard distance genome cov
setwd("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_1")
data <- read.table("Cov.mat", sep= "\t", header=FALSE)
data
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
d <- as.matrix(data[1:8])   # find distance matrix 
hc <- hclust(d)                # apply hirarchical clustering 
plot(hc)
heatmap.2(d, Rowv=as.dendrogram(hc), symm=T, trace="none",
          col=hmcol, margins=c(4,4), main="Pearson correlation of each
          sample")


### For extended Figures

## For num
setwd("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_1/extended")
data <- read.table("pair.num.mat", sep= "\t", header=FALSE)
data
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
d <- as.matrix(data)   # find distance matrix 
hc <- hclust(d)                # apply hirarchical clustering 
plot(hc)
heatmap.2(d, Rowv=as.dendrogram(hc), symm=T, trace="none",
          col=hmcol, margins=c(4,4), main="")

## For boundaries
setwd("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_1/extended")
data <- read.table("pari.boundaries.mat", sep= "\t", header=FALSE)
data
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
d <- as.matrix(data)   # find distance matrix 
hc <- hclust(d)                # apply hirarchical clustering 
plot(hc)
heatmap.2(d, Rowv=as.dendrogram(hc), symm=T, trace="none",
          col=hmcol, margins=c(4,4), main="")

## For cov
setwd("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_1/extended")
data <- read.table("pair.cov.mat", sep= "\t", header=FALSE)
data
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
d <- as.matrix(data)   # find distance matrix 
hc <- hclust(d)                # apply hirarchical clustering 
plot(hc)
heatmap.2(d, Rowv=as.dendrogram(hc), symm=T, trace="none",
          col=hmcol, margins=c(4,4), main="")

########
install.packages("plyr")
library(ggplot2)
library(plyr)
setwd("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_1")
fre<-read.table("TAD.topdom.size.bed.bed",header = T)
fre
mu <- ddply(fre, "Tissues", summarise, grp.mean=mean(Length))

ggplot(fre, aes(x=Length,color=Tissues)) + 
  geom_density(stat="bin",bins=100) + ylab("Count") +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Tissues),
             linetype="dashed") +
  xlim(10000, 8000000)



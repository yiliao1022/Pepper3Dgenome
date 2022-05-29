library("ggplot2")
library("gplots")
library("amap")
library("RColorBrewer")
library(dplyr)
setwd("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/Figure5_b/2022_revision")
data<-read.table(file = "matrixYP2.all.domain.raw.Cluster.new.bed.bed",
                 header = TRUE, sep="\t")
Means <-colMeans(data[sapply(data, is.numeric)],na.rm=TRUE)

MeansVec <- c(Means['Gene'], Means['LTR'], Means['Expression'],Means['Methy1'],Means['Methy2'],Means['H3K4_1'],Means['H3K4_2'],Means['H3K9_1'],Means['H3K9_2'],Means['H3K27_1'],Means['H3K27_2'])
NnormEO <- sweep(data[,2:12], 2, MeansVec, `/`)
logDatEO <-log2(NnormEO+0.01)
logDatEO['TAD']<-data['TAD']

#logDatEO

min_max_norm <- function(x) {
  2*((x - min(x,na.rm=TRUE)) /(max(x,na.rm=TRUE) - min(x,na.rm=TRUE)))-1
}

LogDatEOsub <- logDatEO[c("LTR","Gene","Methy1","H3K4_1","H3K9_1","H3K27_1")]
iris_norm <- as.data.frame(lapply(LogDatEOsub[1:6], min_max_norm ))
otter.matrix <- as.matrix((iris_norm[, 1:6]))


dimnames(otter.matrix)[1] <- list(data$TAD)





dend1 <- as.dendrogram(hclust(dist(otter.matrix,method="euclidean"),method="complete"))

clust<-hclust(dist(otter.matrix, method="euclidean"),method="complete")
groups<-cutree(clust, k=3)
write.table(groups,file="YP2.group.bed")
groups


c_group <- 3
dend1 <- color_branches(dend1, k = c_group, col = rainbow_hcl)
dend1 <- color_labels(dend1, k = c_group, col = rainbow_hcl)  
rMeans <- rowMeans(otter.matrix, na.rm = T)
dend1 <- reorder(dend1, rowMeans(otter.matrix, na.rm = T), agglo.FUN = mean)
col_labels <- get_leaves_branches_col(dend1)
col_labels <- col_labels[order(order.dendrogram(dend1))]
dend1 <- set(dend1, "labels_cex", 0.8)
par(mar = c(1,1,1,14))
plot_horiz.dendrogram(dend1, side = F) 

##########


par(cex.main=0.8)                   # adjust font size of titles
heatmap.2(otter.matrix, #main = '',
          # reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean),
          # order by branch mean so the deepest color is at the top
          dendrogram = "row",        # no dendrogram for columns
          Rowv = dend1,              # * use self-made dendrogram
          Colv = "NA",               # make sure the columns follow data's order
          col = rev(rainbow(20*10, start = 0/6, end = 4/6)),         # color pattern of the heatmap
          
          trace="none",              # hide trace
          density.info="none",       # hide histogram
          
          margins = c(5,18),         # margin on top(bottom) and left(right) side.
          cexRow=1, cexCol = 0.8,      # size of row / column labels
          #xlab = "#",
          srtCol=0, adjCol = c(0.5,1), # adjust the direction of row label to be horizontal
          # margin for the color key
          # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
          key.par=list(mar=c(5,1,3,1)),
          RowSideColors = col_labels, # to add nice colored strips        
          colRow = col_labels         # add color to label
)




###################






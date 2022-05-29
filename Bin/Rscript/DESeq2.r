library("DESeq2")
library("RColorBrewer")
library("gplots")
library("amap")
library("ggplot2")
library("BiocParallel")

setwd("/home/yiliao/Documents/Pepper_2021/CD/500kb")
data <- read.table("CA59.all.lst.out.out.out", header=T, row.names=1, com='', quote='',
                   check.names=F, sep="\t")


sample <- read.table("sampleFile", header=T, row.names=1, com='',
                     quote='', check.names=F, sep="\t", colClasses="factor")

sample <- sample[match(colnames(data), rownames(sample)),, drop=F]
sample_rowname <- rownames(sample)
sample <- data.frame(lapply(sample, function(x) factor(x, levels=unique(x))))
rownames(sample) <- sample_rowname
sample



ddsFullCountTable <- DESeqDataSetFromMatrix(countData = data,
                                            colData = sample,  design= ~ batch + conditions)
dds <- DESeq(ddsFullCountTable)
dds
normalized_counts <- counts(dds, normalized=TRUE)
head(normalized_counts)
normalized_counts 
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
rd <-round(normalized_counts)

# log转换后的结果
rld <- rlog(rd, blind=FALSE)
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]
rMat <- limma::removeBatchEffect(rlogMat,c(sample$batch))

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

# 计算相关性pearson correlation
pearson_cor <- as.matrix(cor(rd, method="pearson"))

# 层级聚类
hc<-hcluster(t(rMat), method="pearson")
hc
#pdf("ehbio_trans.Count_matrix.xls.DESeq2.normalized.rlog.pearson.pdf", pointsize=10)

heatmap.2(pearson_cor, Rowv=as.dendrogram(hc), symm=T, trace="none",
          col=hmcol, margins=c(11,11), main="The pearson correlation of each
sample")

pcadata<-plotPCA(rld, intgroup=c("condition2"), returnData=T, ntop=5000)


percentVar <- round(100*attr(pcadata,"percentVar"),1)
percentVar
ggplot(pcadata, aes(PC1, PC2, color=Group, shape=Batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text_repel(aes(PC1, PC2,color=Group,label=colnames(vsd)),size=3) +
  theme_bw()

plot(pcaData[,1:2],pch=19,col=c("red","red","blue","blue","orange","orange","black","black"))
text(pcaData[,1],pcaData[,2]+0.2,row.names(pcaData),cex=0.5)
pcaData





















LimmaNor <-normalizeCyclicLoess(data, weights = NULL, span=0.7, iterations = 3, method = "pairs")

head(LimmaNor)

LimmaNor[,5]

LJGR1 <- density(LimmaNor[,1])
plot(LJGR1, xlim=c(1,500),main="LJGR1")
polygon(LJGR1, col="red", border="blue")

LJGR2 <- density(LimmaNor[,2])
plot(LJGR2, xlim=c(1,500), main="LJGR2")
polygon(LJGR2, col="red", border="blue")


LJHL1 <- density(LimmaNor[,3])
plot(LJHL1, main="LJHL1")
polygon(LJHL1, col="red", border="blue")

LJHL2 <- density(LimmaNor[,4])
plot(LJHL2, main="LJHL2")
polygon(LJHL2, col="red", border="blue")


LJHL1 <- density(LimmaNor[5,])
plot(LJHL1, main="LJHL1")
polygon(LJHL1, col="red", border="blue")

LJHL2 <- density(LimmaNor[6,])
plot(LJHL2, main="LJHL2")
polygon(LJHL2, col="red", border="blue")

LJHL1 <- density(LimmaNor[7,])
plot(LJHL1, main="LJHL1")
polygon(LJHL1, col="red", border="blue")

LJHL2 <- density(LimmaNor[8,])
plot(LJHL2, main="LJHL2")
polygon(LJHL2, col="red", border="blue")



rMat <- limma::removeBatchEffect(LimmaNor,c(sample$batch))


hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

# 计算相关性pearson correlation
pearson_cor <- as.matrix(cor(rMat, method="pearson"))

# 层级聚类
hc<-hcluster(t(rMat), method="pearson")

heatmap.2(pearson_cor, Rowv=as.dendrogram(hc), symm=T, trace="none",
          col=hmcol, margins=c(11,11), main="The pearson correlation of each
sample")














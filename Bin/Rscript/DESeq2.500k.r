
library("DESeq2")
library("RColorBrewer")
library("gplots")
library("amap")
library("ggplot2")
library("limma")


setwd("/home/yiliao/Documents/Pepper_2021/CD/500kb")
data <- read.table("CA59.all.lst.out.out.out", header=T, row.names=1, com='', quote='', check.names=F, sep="\t")
data <- data[rowSums(data)>5,]
head(data)
dim(data)

sample <- read.table("sampleFile", header=T, row.names=1, com='', quote='', check.names=F, sep="\t", colClasses="factor")
sample <- sample[match(colnames(data), rownames(sample)),, drop=F]
sample_rowname <- rownames(sample)
sample <- data.frame(lapply(sample, function(x) factor(x, levels=unique(x))))
rownames(sample) <- sample_rowname
sample


ddsFullCountTable <- DESeqDataSetFromMatrix(countData = data, colData = sample,  design= ~ batch + conditions)
dds <- DESeq(ddsFullCountTable)


res =results(dds,contrast=c("conditions","LJGR","LJTZ"))
res = res[order(res$pvalue),]
summary(res)

res1 =results(dds,contrast=c("conditions","LJYP","LJTZ"))
res1 = res1[order(res1$pvalue),]
summary(res1)

res2 =results(dds,contrast=c("conditions","LJYP","LJGR"))
res2 = res2[order(res2$pvalue),]
summary(res2)

res3 =results(dds,contrast=c("conditions","LJYP","LJHL"))
res3 = res3[order(res3$pvalue),]
summary(res3)

data_nor<-normalizeCyclicLoess(data, weights = NULL, span=0.7, iterations = 3, method = "fast")

rMat <- limma::removeBatchEffect(dds,c(sample$batch))

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

# 计算相关性pearson correlation
pearson_cor <- as.matrix(cor(rMat, method="pearson"))

# 层级聚类
hc<-hcluster(t(rMat), method="pearson")
hc
#pdf("ehbio_trans.Count_matrix.xls.DESeq2.normalized.rlog.pearson.pdf", pointsize=10)

heatmap.2(pearson_cor, Rowv=as.dendrogram(hc), symm=T, trace="none",
          col=hmcol, margins=c(11,11), main="The pearson correlation of each sample")

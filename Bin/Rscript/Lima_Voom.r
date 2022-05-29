if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")

library(edgeR)


setwd("/home/yiliao/Documents/Pepper_2021/geneExp")


counts <- read.delim("Ca59.gene.tsv", row.names = 1)
head(counts)
dim(counts)
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
cutoff <- 8
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
dim(d)
snames <- colnames(counts)
snames
tissues <- substr(snames,1,4)
tissues
plotMDS(d, col = as.numeric(tissues))

mm <- model.matrix(~0 + tissues)
y <- voom(d, mm, plot = T)

fit <- lmFit(y, mm)
head(coef(fit))

contr <- makeContrasts(tissuesLJGR - tissuesLJHL, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
top.table
head(top.table, 5977)
write.table(top.table,"LJGRvsLJHL.txt")

length(which(top.table$adj.P.Val < 0.05))




##############################################################
setwd("/home/yiliao/Documents/Pepper_2021/geneExp")


counts <- read.delim("ind.lst.out.out", row.names = 1)
head(counts)
dim(counts)
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
dim(d)
snames <- colnames(counts)
snames
tissues <- substr(snames,1,5)
tissues
plotMDS(d, col = as.numeric(tissues))

mm <- model.matrix(~0 + tissues)
y <- voom(d, mm, plot = T)

fit <- lmFit(y, mm)
head(coef(fit))

contr <- makeContrasts(tissuesLJGR - tissuesLJHL, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 200)
length(which(top.table$adj.P.Val < 0.05))








########################







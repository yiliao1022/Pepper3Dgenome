#install.packages("DESeq2")
library("edgeR")
library("dplyr")
####library("DESeq2")
        
setwd("./")
counts <- read.delim("featureCounts.rawcount.txt.addsite.bed.bed", row.names = 1)
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
CPM <-fpkm(d0)
cutoff <- 0.1
drop <- which(apply(fpkm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
CPM<-fpkm(d)
write.csv(CPM,file="CA59.RNA.40k.cpm.csv")
snames <- colnames(counts) # Sample names
tissues <- substr(snames,1,4)
#batch <- factor(x=c("A", "B", "A", "A", "A", "B", "A","B"),levels = c("A", "B"))
#mm <- model.matrix(~0 + tissues + batch)
mm <- model.matrix(~0 + tissues)
y <- voom(d, mm, plot = F)
fit <- lmFit(y, mm)
#head(coef(fit))
contr1 <- makeContrasts(tissuesLJYP - tissuesLJGR, levels = colnames(coef(fit)))
tmp1 <- contrasts.fit(fit, contr1)
tmp1 <- eBayes(tmp1)
top.table1 <- topTable(tmp1, sort.by = "P", n = Inf)
#head(top.table, 200)
#length(which(top.table$adj.P.Val < 0.05))
write.csv(top.table1,file="LJYPvsLJGR.RNA.all.csv")
rm(contr1)
rm(tmp1)
rm(top.table1)

###
contr2 <- makeContrasts(tissuesLJYP - tissuesLJHL, levels = colnames(coef(fit)))
tmp2 <- contrasts.fit(fit, contr2)
tmp2 <- eBayes(tmp2)
top.table2 <- topTable(tmp2, sort.by = "P", n = Inf)
write.csv(top.table2,file="LJYPvsLJHL.RNA.all.csv")
rm(contr2)
rm(tmp2)
rm(top.table2)

###
contr3 <- makeContrasts(tissuesLJYP - tissuesLJTZ, levels = colnames(coef(fit)))
tmp3 <- contrasts.fit(fit, contr3)
tmp3 <- eBayes(tmp3)
top.table3 <- topTable(tmp3, sort.by = "P", n = Inf)
write.csv(top.table3,file="LJYPvsLJTZ.RNA.all.csv")
rm(contr3)
rm(tmp3)
rm(top.table3)

###
contr4 <- makeContrasts(tissuesLJHL - tissuesLJTZ, levels = colnames(coef(fit)))
tmp4 <- contrasts.fit(fit, contr4)
tmp4 <- eBayes(tmp4)
top.table4 <- topTable(tmp4, sort.by = "P", n = Inf)
write.csv(top.table4,file="LJHLvsLJTZ.RNA.all.csv")
rm(contr4)
rm(tmp4)
rm(top.table4)

###
contr5 <- makeContrasts(tissuesLJHL - tissuesLJGR, levels = colnames(coef(fit)))
tmp5 <- contrasts.fit(fit, contr5)
tmp5 <- eBayes(tmp5)
top.table5 <- topTable(tmp5, sort.by = "P", n = Inf)
write.csv(top.table5,file="LJHLvsLJGR.RNA.all.csv")
rm(contr5)
rm(tmp5)
rm(top.table5)

###
contr6 <- makeContrasts(tissuesLJTZ - tissuesLJGR, levels = colnames(coef(fit)))
tmp6 <- contrasts.fit(fit, contr6)
tmp6 <- eBayes(tmp6)
top.table6 <- topTable(tmp6, sort.by = "P", n = Inf)
write.csv(top.table6,file="LJTZvsLJGR.RNA.all.csv")
rm(contr6)
rm(tmp6)
rm(top.table6)






























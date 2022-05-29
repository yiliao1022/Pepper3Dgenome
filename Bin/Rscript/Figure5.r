#################Figure5
library("edgeR")
library("ggplot2")
library("dplyr")
setwd("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD")
rawdata <- read.table("RNA.40k.5.bed.bed.bed.bed", header=T, row.names=1, com='', quote='',
                      check.names=F, sep="\t")
cpmdata <- edgeR::cpm(rawdata)
write.table(cpmdata,"RNA.40K.CPM.5.bed",sep="\t")

#######################Figure5
install.packages("PairedData")
install.packages("ggpubr")
setwd("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD")
library("ggpubr")
my_data<-read.table("Genecomparison.bed.bed",header=T)
ggboxplot(my_data, x = "Types", y = "Density", 
          color = "Types", palette = c("#00AFBB", "#E7B800"),
          order = c("interTAD", "withinTAD"),
          ylab = "Density", xlab = "Types")
res <- wilcox.test(Density ~ Types, data = my_data, paired = F)

my_data1<-read.table("TEcomparison.bed",header=T)
ggboxplot(my_data1, x = "Types", y = "Density", 
          color = "Types", palette = c("#00AFBB", "#E7B800"),
          order = c("interTAD", "withinTAD"),
          ylab = "Density", xlab = "Types")
res <- wilcox.test(Density ~ Types, data = my_data1, paired = F)

my_data2<-read.table("RNA.tau.comparsion.bed",header=T)
ggboxplot(my_data2, x = "Types", y = "Density", 
          color = "Types", palette = c("#00AFBB", "#E7B800"),
          order = c("interTAD", "withinTAD"),
          ylab = "tau", xlab = "Types")
res <- wilcox.test(Density ~ Types, data = my_data2, paired = F)
res

my_data3<-read.table("SpecificUniqinterComparison.bed",header=T)
ggboxplot(my_data3, x = "Types", y = "Density", 
          color = "Types", palette = c("#00AFBB", "#E7B800","orange"),
          order = c("Specific", "Shared","Inter"),
          ylab = "tau", xlab = "Types")
res <- wilcox.test(Density ~ Types, data = my_data3, paired = F)
res

my_data4<-read.table("GR1boundariesCompasion.bed",header=T)
ggboxplot(my_data4, x = "Types", y = "Density", 
          color = "Types", palette = c("#00AFBB", "#E7B800"),
          order = c("Specific", "Shared"),
          ylab = "tau", xlab = "Types")
res <- wilcox.test(Density ~ Types, data = my_data4, paired = F)
res

my_data4<-read.table("YP1.tau.Comparison.bed",header=T)
ggboxplot(my_data4, x = "Types", y = "Density", 
          color = "Types", palette = c("#00AFBB", "#E7B800"),
          order = c("Specific", "Shared"),
          ylab = "tau", xlab = "Types")
res <- wilcox.test(Density ~ Types, data = my_data4, paired = F)
res


















my_data5<-read.table("GRRNAfoldComparison.bed.bed",header=T)
ggboxplot(my_data5, x = "Types", y = "GRvsHL", 
          color = "Types", palette = c("#00AFBB", "#E7B800","orange"),
          order = c("Specific", "Shared","interTAD"),
          ylab = "RNA fold change", xlab = "Types")

my_data6 <- my_data5 %>% filter(Types=="Specific" | Types=="Shared")
my_data6
res <- wilcox.test(TZvsYP ~ Types, data = my_data6, paired = F)
res

my_data4<-read.table("RNAboundariesComparison.bed.bed",header=T)
ggboxplot(my_data4, x = "Types", y = "GRvsTZ", 
          color = "Types", palette = c("#00AFBB", "#E7B800"),
          order = c("Specific", "Shared"),
          ylab = "tau", xlab = "Types")
res <- wilcox.test(GRvsYP ~ Types, data = my_data4, paired = F)
res

my_data4<-read.table("GRvsHL.boundaries.RNA.bed.bed",header=T)
ggboxplot(my_data4, x = "Types", y = "GRvsHL", 
          color = "Types", palette = c("#00AFBB", "#E7B800"),
          order = c("Overlap", "Specific"),
          ylab = "tau", xlab = "Types")
res <- wilcox.test(GRvsHL ~ Types, data = my_data4, paired = F)
res

my_data4<-read.table("PairbounvsRNA.bed.bed",header=T)
ggboxplot(my_data4, x = "Types", y = "GRvsTZ", 
          color = "Types", palette = c("#00AFBB", "#E7B800"),
          order = c("Shared", "Specific"),
          ylab = "tau", xlab = "Types")
res <- wilcox.test(GRvsTZ ~ Types, data = my_data4, paired = F)
res

###############Figure5c


##TAD number
install.packages("VennDiagram")
library("VennDiagram")
venn.plot <- draw.triple.venn(
  area1 = 1680,
  area2 = 4663,
  area3 = 2641,
  n12 = 669,
  n23 = 1241,
  n13 = 683,
  n123 = 341,
  category = c("Arrowhead", "HiCExplorer", "TopDom"),
  fill = c("blue", "red", "green"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("blue", "red", "green")
);

###TAD genome coverage
install.packages("VennDiagram")
library("VennDiagram")
venn.plot <- draw.triple.venn(
  area1 = 1677,
  area2 = 3038,
  area3 = 3043,
  n12 = 570,
  n23 = 905,
  n13 = 776,
  n123 = 294,
  category = c("Arrowhead", "HiCExplorer", "TopDom"),
  fill = c("blue", "red", "green"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("blue", "red", "green")
);


### TAD boundaries 
install.packages("VennDiagram")
library("VennDiagram")
venn.plot <- draw.triple.venn(
  area1 = 3239,
  area2 = 4675,
  area3 = 2699,
  n12 = 708,
  n23 = 1901,
  n13 = 1164,
  n123 = 438,
  category = c("Arrowhead", "HiCExplorer", "TopDom"),
  fill = c("blue", "red", "green"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("blue", "red", "green")
);


#################Figure5d
library(ggplot2)
library(plyr)
setwd("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/figure5d")
dataSize<-read.table("All.3.size.bed",head=T)
dataSize1 <- dataSize %>% filter (Tissues=='YP1')

mu <- ddply(dataSize1, "Methods", summarise, grp.mean=mean(Sizes))

ggplot(dataSize1, aes(x=Sizes,color=Methods)) + 
  geom_density() +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Methods),
             linetype="dashed") +
  xlim(1000, 6000000)

###################Figure5e
jfunc <- function(x) { s <- sd(x); m <- mean(x); return(c(m+c(-1, 0, 1)*s)) }
setwd("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/figure5e")
dat1 <- matrix(as.numeric(scan(file = 'Pepper_gene_40k.density.background.bed')), nc = 51)
rs1 = rowSums(dat1)
tmp1 <- apply(X=dat1, MARGIN = 2, function(x) x) 
ans1 <- apply(X = tmp1, MARGIN = 2, FUN = jfunc)
HQ <- as.numeric(scan(file = 'Pepper_gene_40k.density.peak.cov.bed.bed'))
plot(xlab="Distance to TAD boundaries (kp)",ylab="Gene content per 40-kb bin [%]",-1,-1,ylim = c(0,0.08), xlim = c(0,51)) ##
polygon(y = t(ans1[c(1,3),]), x = c(1:51, 51:1), col = 'gray', border = NA)
lines(ans1[2,],col='maroon',lty="dashed")
lines(HQ, col = 'black', lwd=1)
abline(v=26,lty="dashed")
abline(v=27,lty="dashed")



dat1 <- matrix(as.numeric(scan(file = 'Pepper_LTR_40k.density.bed.background.bed')), nc = 51)
dat1
rs1 = rowSums(dat1)
tmp1 <- apply(X=dat1, MARGIN = 2, function(x) x) 
ans1 <- apply(X = tmp1, MARGIN = 2, FUN = jfunc)
tmp1
ans1
HQ <- as.numeric(scan(file = 'Pepper_LTR_40k.density.bed.out'))
plot(xlab="Distance to TAD boundaries (kp)",ylab="LTR content per 40-kb bin [%]",-1,-1,ylim = c(0.65,0.8), xlim = c(0,51)) ##
polygon(y = t(ans1[c(1,3),]), x = c(1:51, 51:1), col = 'gray', border = NA)
lines(ans1[2,],col='maroon',lty="dashed")
lines(HQ, col = 'black', lwd=1)
abline(v=26,lty="dashed")
abline(v=27,lty="dashed")




################Figure5f
install.packages("vegan")
install.packages("pheatmap")
library("vegan")
library("pheatmap")

jaccard_mat <- read.table ("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/figure5f/Jaccard.matrixnew")
jaccard_mat<-as.matrix(jaccard_mat)
hc<-hcluster(jaccard_mat, method="pearson")
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

heatmap.2(jaccard_mat, Rowv=as.dendrogram(hc), symm=T, trace="none",
          col=hmcol, margins=c(4,4), main="Pearson correlation of each
          sample")
################Figure5h
mydata1 <-read.table("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/figure5h/LTR.density.Comparison.bed.bed",head=T)
head(mydata1)
mydata1$Tissues <- as.factor(mydata1$Tissues )
ggplot(mydata1, aes(x=Tissues, y=Density, fill=factor(Types))) + geom_boxplot(outlier.shape = NA) + ylim(0,1)+labs(fill = "Types")

mydata1_filter <-mydata1 %>% filter(Tissues=='GR') %>% filter(Types=='InterTAD' | Types=='Shared')
wilcox.test(Density~Types, data = mydata1_filter, paired = FALSE)

mydata1_filter <-mydata1 %>% filter(Tissues=='HL') %>% filter(Types=='InterTAD' | Types=='Shared')
wilcox.test(Density~Types, data = mydata1_filter, paired = FALSE)

mydata1_filter <-mydata1 %>% filter(Tissues=='TZ') %>% filter(Types=='InterTAD' | Types=='Shared')
wilcox.test(Density~Types, data = mydata1_filter, paired = FALSE)

mydata1_filter <-mydata1 %>% filter(Tissues=='YP') %>% filter(Types=='InterTAD' | Types=='Shared')
wilcox.test(Density~Types, data = mydata1_filter, paired = FALSE)


mydata2 <-read.table("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/figure5h/Gene.density.Comparison.bed.bed.bed",head=T)
head(mydata2)
mydata2$Tissues <- as.factor(mydata2$Tissues )
ggplot(mydata2, aes(x=Tissues, y=Density, fill=factor(Types))) + geom_boxplot(outlier.shape = NA) + ylim(0,0.6)+labs(fill = "Types")

mydata2_filter <-mydata2 %>% filter(Tissues=='GR') %>% filter(Types=='InterTAD' | Types=='Shared')
wilcox.test(Density~Types, data = mydata2_filter, paired = FALSE)

mydata2_filter <-mydata2 %>% filter(Tissues=='HL') %>% filter(Types=='InterTAD' | Types=='Shared')
wilcox.test(Density~Types, data = mydata2_filter, paired = FALSE)

mydata2_filter <-mydata2 %>% filter(Tissues=='TZ') %>% filter(Types=='InterTAD' | Types=='Shared')
wilcox.test(Density~Types, data = mydata2_filter, paired = FALSE)

mydata2_filter <-mydata2 %>% filter(Tissues=='YP') %>% filter(Types=='InterTAD' | Types=='Shared')
wilcox.test(Density~Types, data = mydata2_filter, paired = FALSE)
########################


###############Figure5i

##Left
mydata3 <-read.table("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/figure5i/Tau.density.Comparison.bed.bed",head=T)
head(mydata3)
mydata3$Tissues <- as.factor(mydata3$Tissues )
ggplot(mydata3, aes(x=Tissues, y=Density, fill=factor(Types))) + geom_boxplot(outlier.shape = NA) + ylim(0,1)+labs(fill = "Types")

mydata3_filter <-mydata3 %>% filter(Tissues=='GR') %>% filter(Types=='InterTAD' | Types=='Specific')
wilcox.test(Density~Types, data = mydata3_filter, paired = FALSE)
mydata3_filter <-mydata3 %>% filter(Tissues=='GR') %>% filter(Types=='Shared' | Types=='Specific')
wilcox.test(Density~Types, data = mydata3_filter, paired = FALSE)


mydata3_filter <-mydata3 %>% filter(Tissues=='HL') %>% filter(Types=='InterTAD' | Types=='Specific')
wilcox.test(Density~Types, data = mydata3_filter, paired = FALSE)
mydata3_filter <-mydata3 %>% filter(Tissues=='HL') %>% filter(Types=='Shared' | Types=='Specific')
wilcox.test(Density~Types, data = mydata3_filter, paired = FALSE)

mydata3_filter <-mydata3 %>% filter(Tissues=='TZ') %>% filter(Types=='InterTAD' | Types=='Specific')
wilcox.test(Density~Types, data = mydata3_filter, paired = FALSE)
mydata3_filter <-mydata3 %>% filter(Tissues=='TZ') %>% filter(Types=='Shared' | Types=='Specific')
wilcox.test(Density~Types, data = mydata3_filter, paired = FALSE)

mydata3_filter <-mydata3 %>% filter(Tissues=='YP') %>% filter(Types=='InterTAD' | Types=='Specific')
wilcox.test(Density~Types, data = mydata3_filter, paired = FALSE)
mydata3_filter <-mydata3 %>% filter(Tissues=='YP') %>% filter(Types=='Shared' | Types=='Specific')
wilcox.test(Density~Types, data = mydata3_filter, paired = FALSE)


mydata3 <-read.table("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/figure5i/Boundaries.tau.topdom.out.out",head=T)
head(mydata3)
mydata3$Tissues <- as.factor(mydata3$Tissues )
ggplot(mydata3, aes(x=Tissues, y=Density, fill=factor(Types))) + geom_boxplot(outlier.shape = NA) + ylim(-0.1,1.1)+labs(fill = "Types")

mydata3_filter <-mydata3 %>% filter(Tissues=='GR')  %>% filter(Types=='conserved' | Types=='specific')
wilcox.test(Density~Types, data = mydata3_filter, paired = FALSE) 
head(mydata3_filter)
t.test(Density~Types, data = mydata3_filter, var.equal = FALSE)





mydata3 <-read.table("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/figure5i/Boundaries.fold.new.arrowhead.out.out",head=T)
head(mydata3)
mydata3$Tissues <- as.factor(mydata3$Tissues )
ggplot(mydata3, aes(x=Tissues, y=Density, fill=factor(Types))) + geom_boxplot(outlier.shape = NA) + ylim(0,5)+labs(fill = "Types")


mydata3_filter <-mydata3 %>% filter(Tissues=='YPvsHL') %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter, alternative = "less", paired = FALSE)
head(mydata3_filter)
t.test(Density~Types, data = mydata3_filter, var.equal = FALSE)


###tmp
mydata3 <-read.table("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/figure5i/YP.tad.tmp.bed.bed",head=T)
head(mydata3)
mydata3$Tissues <- as.factor(mydata3$Tissues )
ggplot(mydata3, aes(x=Tissues, y=Density, fill=factor(Types))) + geom_boxplot(outlier.shape = NA) + ylim(0,1)+labs(fill = "Types")
mydata3_filter <-mydata3 %>% filter(Tissues=='YP') %>% filter(Types=='specific' | Types=='shared')
wilcox.test(Density~Types, data = mydata3_filter,  alternative = "less", paired = FALSE)
t.test(Density~Types, data = mydata3_filter,  alternative = "less", var.equal = FALSE)

###Right

mydata4 <-read.table("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/figure5i/Tau.boundaries.density.Comparison.bed.bed",head=T)
head(mydata4)
mydata4$Tissues <- as.factor(mydata4$Tissues )
ggplot(mydata4, aes(x=Tissues, y=Density, fill=factor(Types))) + geom_boxplot(outlier.shape = NA) + ylim(0,1)+labs(fill = "Types")

mydata4_filter <-mydata4 %>% filter(Tissues=='GR') %>% filter(Types=='InterTAD' | Types=='Specific')
wilcox.test(Density~Types, data = mydata4_filter, paired = FALSE)
mydata4_filter <-mydata4 %>% filter(Tissues=='GR') %>% filter(Types=='InterTAD' | Types=='Shared')
wilcox.test(Density~Types, data = mydata4_filter, paired = FALSE)
mydata4_filter <-mydata4 %>% filter(Tissues=='GR') %>% filter(Types=='Shared' | Types=='Specific')
wilcox.test(Density~Types, data = mydata4_filter, paired = FALSE)


mydata4_filter <-mydata4 %>% filter(Tissues=='HL') %>% filter(Types=='InterTAD' | Types=='Specific')
wilcox.test(Density~Types, data = mydata4_filter, paired = FALSE)
mydata4_filter <-mydata4 %>% filter(Tissues=='HL') %>% filter(Types=='InterTAD' | Types=='Shared')
wilcox.test(Density~Types, data = mydata4_filter, paired = FALSE)
mydata4_filter <-mydata4 %>% filter(Tissues=='HL') %>% filter(Types=='Shared' | Types=='Specific')
wilcox.test(Density~Types, data = mydata4_filter, paired = FALSE)

mydata4_filter <-mydata4 %>% filter(Tissues=='TZ') %>% filter(Types=='InterTAD' | Types=='Specific')
wilcox.test(Density~Types, data = mydata4_filter, paired = FALSE)
mydata4_filter <-mydata4 %>% filter(Tissues=='TZ') %>% filter(Types=='InterTAD' | Types=='Shared')
wilcox.test(Density~Types, data = mydata4_filter, paired = FALSE)
mydata4_filter <-mydata4 %>% filter(Tissues=='TZ') %>% filter(Types=='Shared' | Types=='Specific')
wilcox.test(Density~Types, data = mydata4_filter, paired = FALSE)

mydata4_filter <-mydata4 %>% filter(Tissues=='YP') %>% filter(Types=='InterTAD' | Types=='Specific')
wilcox.test(Density~Types, data = mydata4_filter, paired = FALSE)
mydata4_filter <-mydata4 %>% filter(Tissues=='YP') %>% filter(Types=='InterTAD' | Types=='Shared')
wilcox.test(Density~Types, data = mydata4_filter, paired = FALSE)
mydata4_filter <-mydata4 %>% filter(Tissues=='YP') %>% filter(Types=='Shared' | Types=='Specific')
wilcox.test(Density~Types, data = mydata4_filter, paired = FALSE)


#################




##################Figure4 b
library("ggplot2")
library("gplots")
library("amap")
library("RColorBrewer")

setwd("/home/yiliao/Documents/Pepper_2021/bnbc")
fre<-read.table("All.chr01.txt.txt.2m.out")
pearson_cor <- as.matrix(cor(fre, method="pearson"))

hc<-hcluster(t(fre), method="pearson")
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(pearson_cor, Rowv=as.dendrogram(hc), symm=T, trace="none",
          col=hmcol, margins=c(4,4), main="Pearson correlation of each
          sample")
##########################################
##Figure4c Excel
####################
##Figure4d
library(ggplot2)
library(plyr)
setwd("/home/yiliao/Documents/Pepper_2021/ABcompartments/draw/FINAL/All")
fre<-read.table("All.domain.length.bed.bed",header = T)
fre
mu <- ddply(fre, "Tissues", summarise, grp.mean=mean(Length))

ggplot(fre, aes(x=Length,color=Tissues)) + 
  geom_density() +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Tissues),
                             linetype="dashed") +
  xlim(10000, 1000000)


library(ggplot2)
library(plyr)
setwd("/home/yiliao/Documents/Pepper_2021/ABcompartments/draw/2022_revision/4papper")
fre<-read.table("All.length.bed.bed",header = T)
fre
mu <- ddply(fre, "Tissues", summarise, grp.mean=mean(Length))

ggplot(fre, aes(x=Length,color=Tissues)) + 
  geom_density() +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Tissues),
             linetype="dashed") +
  xlim(9000, 1000000)


################################################
### Figure4e
setwd("/home/yiliao/Documents/Pepper_2021/ABcompartments/draw/FINAL/All")

library("edgeR")
library("ggplot2")
library("dplyr")

mydata1 <-read.table("ABvsRNA.bed.bed.bed.1.bed.bed",head=T)
head(mydata1)
mydata1$Rank <- as.factor(mydata1$Rank )
ggplot(mydata1, aes(x=Rank, y=log2(RNA), fill=factor(Tissues))) + geom_boxplot(outlier.shape = NA) + ylim(0,12)+labs(fill = "Tissues") 

###################################################
###Figure4f perl

###Figure4g perl
##################################################
###Figure4h

#####
#mydata2 <-read.table("/home/yiliao/Documents/Pepper_2021/ABcompartments/draw/FINAL/All/20211015/GeneBased/Tworeplicate/ABvsRNA.one.gene.bed.bed",head=T)
#head(mydata2)
#mydata2$Tissues <- as.factor(mydata2$Tissues )
#ggplot(mydata2, aes(x=Tissues, y=RNAFold, fill=factor(Rank))) + geom_boxplot(outlier.shape = NA) + ylim(-10,10)+labs(fill = "Rank")

####
#mydata2 <-read.table("~/Documents/Pepper_2021/ABcompartments/draw/FINAL/All/Final/ABvsRNA.bed.bed.dif.bed",head=T)
#head(mydata2)
#mydata2$Tissues <- as.factor(mydata2$Tissues )
#ggplot(mydata2, aes(x=Tissues, y=RNAFOLD, fill=factor(Rank))) + geom_boxplot(outlier.shape = NA) + ylim(-10,10)+labs(fill = "Rank")

#mydata3 <-mydata2 %>% filter(Tissues=='GRvsHL') %>% filter(Rank==-0.25 | Rank==0)
#wilcox.test(RNAFOLD~Rank, data = mydata3, paired = FALSE,
#            alternative = "less")

#mydata4 <-mydata2 %>% filter(Tissues=='GRvsHL') %>% filter(Rank==0 | Rank==0.25)
#wilcox.test(RNAFOLD~Rank, data = mydata4, paired = FALSE,
#            alternative = "less")

#mydata5 <-mydata2 %>% filter(Tissues=='GRvsYP') %>% filter(Rank==-0.25 | Rank==0)
#wilcox.test(RNAFOLD~Rank, data = mydata5, paired = FALSE,
#            alternative = "less")

#mydata6 <-mydata2 %>% filter(Tissues=='GRvsYP') %>% filter(Rank==0 | Rank==0.25)
#wilcox.test(RNAFOLD~Rank, data = mydata6, paired = FALSE,
#            alternative = "less")

#mydata7 <-mydata2 %>% filter(Tissues=='HLvsYP') %>% filter(Rank==-0.25 | Rank==0)
#wilcox.test(RNAFOLD~Rank, data = mydata7, paired = FALSE,
#            alternative = "less")

#mydata8 <-mydata2 %>% filter(Tissues=='HLvsYP') %>% filter(Rank==0 | Rank==0.25)
#wilcox.test(RNAFOLD~Rank, data = mydata8, paired = FALSE,
 #           alternative = "less")

#mydata9 <-mydata2 %>% filter(Tissues=='TZvsGR') %>% filter(Rank==-0.25 | Rank==0)
#wilcox.test(RNAFOLD~Rank, data = mydata9, paired = FALSE,
 #           alternative = "less")

#mydata10 <-mydata2 %>% filter(Tissues=='TZvsGR') %>% filter(Rank==0 | Rank==0.25)
#wilcox.test(RNAFOLD~Rank, data = mydata10, paired = FALSE,
 #           alternative = "less")

#mydata11 <-mydata2 %>% filter(Tissues=='TZvsHL') %>% filter(Rank==-0.25 | Rank==0)
#wilcox.test(RNAFOLD~Rank, data = mydata11, paired = FALSE,
#            alternative = "less")

#mydata12 <-mydata2 %>% filter(Tissues=='TZvsHL') %>% filter(Rank==0 | Rank==0.25)
#wilcox.test(RNAFOLD~Rank, data = mydata12, paired = FALSE,
#            alternative = "less")

#mydata13 <-mydata2 %>% filter(Tissues=='TZvsYP') %>% filter(Rank==-0.25 | Rank==0)
#wilcox.test(RNAFOLD~Rank, data = mydata13, paired = FALSE,
#            alternative = "less")

#mydata14 <-mydata2 %>% filter(Tissues=='TZvsYP') %>% filter(Rank==0 | Rank==0.25)
#wilcox.test(RNAFOLD~Rank, data = mydata14, paired = FALSE,
#            alternative = "less")

#
mydata2 <-read.table("/home/yiliao/Documents/Pepper_2021/ABcompartments/draw/FINAL/All/20211015/GeneBased/Tworeplicate/Pip3/ABvsRNA.two.gene.bed.bed", head=T)

mydata2 <-read.table("/home/yiliao/Documents/Pepper_2021/ABcompartments/draw/FINAL/All/20211015/GeneBased/Tworeplicate/3DvsGeneExp/ABvsRNA.two.gene.bed.bed", head=T)

mydata2 <-read.table("/home/yiliao/Documents/Pepper_2021/ABcompartments/draw/FINAL/All/20211015/GeneBased/Tworeplicate/3DvsGeneExp/overall/ABvsRNA.two.gene.bed.bed", head=T)

mydata2 <-read.table("/home/yiliao/Documents/Pepper_2021/ABcompartments/draw/FINAL/All/20211015/BinBased/Tworeplicate/001/ABvsRNA.two.gene.bed.bed", head=T)

mydata2 <-read.table("/home/yiliao/Documents/Pepper_2021/ABcompartments/draw/FINAL/All/20211015/BinBased/Onereplicate/Highcoverage/ABvsRNA.one.gene.bed.bed", head=T)

mydata2 <-read.table("/home/yiliao/Documents/Pepper_2021/ABcompartments/draw/FINAL/All/20211015/GeneBased/Onereplicate/Highcoverage/ABvsRNA.one.gene.bed.bed", head=T)

mydata2 <-read.table("/home/yiliao/Documents/Pepper_2021/ABcompartments/draw/FINAL/All/20211015/GeneBased/Tworeplicate/005/ABvsRNA.two.gene.bed.bed", head=T)


mydata2 <-read.table("/home/yiliao/Documents/Pepper_2021/ABcompartments/draw/FINAL/All/20211015/GeneBased/NEW/One/ABvsRNA.one.gene.bed.bed", head=T)
head(mydata2)
mydata2$Tissues <- as.factor(mydata2$Tissues )
mydata9 <-mydata2 %>% filter(Pvalue < 0.01) 
ggplot(mydata9, aes(x=Tissues, y=RNAFold, fill=factor(Rank))) + geom_boxplot(outlier.shape = NA) + ylim(-15,15)+labs(fill = "Rank")

mydata3 <-mydata2 %>% filter(Tissues=='HLvsGR') %>% filter(Rank=="down" | Rank=="stable")
wilcox.test(RNAFold~Rank, data = mydata3, paired = FALSE,
            alternative = "less")
mydata4 <-mydata2 %>% filter(Tissues=='HLvsGR') %>% filter(Rank=="stable" | Rank=="up")
wilcox.test(RNAFold~Rank, data = mydata4, paired = FALSE,
            alternative = "less")
mydata4 <-mydata2 %>% filter(Tissues=='HLvsGR') %>% filter(Rank=="down" | Rank=="up")
wilcox.test(RNAFold~Rank, data = mydata4, paired = FALSE,
            alternative = "less")

mydata3 <-mydata2 %>% filter(Tissues=='HLvsTZ') %>% filter(Rank=="down" | Rank=="stable")
wilcox.test(RNAFold~Rank, data = mydata3, paired = FALSE,
            alternative = "less")
mydata4 <-mydata2 %>% filter(Tissues=='HLvsTZ') %>% filter(Rank=="stable" | Rank=="up")
wilcox.test(RNAFold~Rank, data = mydata4, paired = FALSE,
            alternative = "less")
mydata4 <-mydata2 %>% filter(Tissues=='HLvsTZ') %>% filter(Rank=="down" | Rank=="up")
wilcox.test(RNAFold~Rank, data = mydata4, paired = FALSE,
            alternative = "less")

mydata3 <-mydata2 %>% filter(Tissues=='TZvsGR') %>% filter(Rank=="down" | Rank=="stable")
wilcox.test(RNAFold~Rank, data = mydata3, paired = FALSE,
            alternative = "less")
mydata4 <-mydata2 %>% filter(Tissues=='TZvsGR') %>% filter(Rank=="stable" | Rank=="up")
wilcox.test(RNAFold~Rank, data = mydata4, paired = FALSE,
            alternative = "less")
mydata4 <-mydata2 %>% filter(Tissues=='TZvsGR') %>% filter(Rank=="down" | Rank=="up")
wilcox.test(RNAFold~Rank, data = mydata4, paired = FALSE,
            alternative = "less")

mydata3 <-mydata2 %>% filter(Tissues=='YPvsGR') %>% filter(Rank=="down" | Rank=="stable")
wilcox.test(RNAFold~Rank, data = mydata3, paired = FALSE,
            alternative = "less")
mydata4 <-mydata2 %>% filter(Tissues=='YPvsGR') %>% filter(Rank=="stable" | Rank=="up")
wilcox.test(RNAFold~Rank, data = mydata4, paired = FALSE,
            alternative = "less")
mydata4 <-mydata2 %>% filter(Tissues=='YPvsGR') %>% filter(Rank=="down" | Rank=="up")
wilcox.test(RNAFold~Rank, data = mydata4, paired = FALSE,
            alternative = "less")

mydata3 <-mydata2 %>% filter(Tissues=='YPvsHL') %>% filter(Rank=="down" | Rank=="stable")
wilcox.test(RNAFold~Rank, data = mydata3, paired = FALSE,
            alternative = "less")
mydata4 <-mydata2 %>% filter(Tissues=='YPvsHL') %>% filter(Rank=="stable" | Rank=="up")
wilcox.test(RNAFold~Rank, data = mydata4, paired = FALSE,
            alternative = "less")
mydata4 <-mydata2 %>% filter(Tissues=='YPvsHL') %>% filter(Rank=="down" | Rank=="up")
wilcox.test(RNAFold~Rank, data = mydata4, paired = FALSE,
            alternative = "less")

mydata3 <-mydata2 %>% filter(Tissues=='YPvsTZ') %>% filter(Rank=="down" | Rank=="stable")
wilcox.test(RNAFold~Rank, data = mydata3, paired = FALSE,
            alternative = "less")
mydata4 <-mydata2 %>% filter(Tissues=='YPvsTZ') %>% filter(Rank=="stable" | Rank=="up")
wilcox.test(RNAFold~Rank, data = mydata4, paired = FALSE,
            alternative = "less")
mydata4 <-mydata2 %>% filter(Tissues=='YPvsTZ') %>% filter(Rank=="down" | Rank=="up")
wilcox.test(RNAFold~Rank, data = mydata4, paired = FALSE,
            alternative = "less")


####Figure4i
library(ggplot2)
library(plyr)
library(lattice)
mydata4 <-read.table("/home/yiliao/Documents/Pepper_2021/ABcompartments/draw/FINAL/All/Final/ABvsRNA.bed.bed.RNA.dif.bed",head=T)
mydata4 
mydata4$RNAFOLDRank <- as.factor(mydata4$RNAFOLDRank)
mydata4$RankChange <- as.numeric(mydata4$RankChange)
mydata4$RankChange[mydata4$RankChange>0.2]<-0.25
mydata4$RankChange[mydata4$RankChange<(-0.2)]<-(-0.25)

mydata5 <-mydata4 %>% filter(RNAFOLDRank=='stable')
mydata6 <-mydata4 %>% filter(RNAFOLDRank=='up')
mydata7 <-mydata4 %>% filter(RNAFOLDRank=='down')
dim(mydata5)
dim(mydata6)
dim(mydata7)
mydata4
 mydata8<-mydata4 %>% filter(Tissues=='GRvsYP' | Tissues=='TZvsYP' | Tissues =='HLvsYP')
 mydata8
ggplot(mydata8, aes(x= RankChange, group=RNAFOLDRank)) + 
  geom_bar(aes(y = ..prop.., fill =factor(..x..)), stat="count") +
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.5) +
  labs(y = "Percent", fill="RankChange") +
  facet_grid(~RNAFOLDRank) +
  scale_y_continuous(labels = scales::percent)

wilcox.test(RankChange~RNAFOLDRank, data = mydata5, paired = FALSE,
            alternative = "less")


########################


####################################
library(ggplot2)
library(plyr)
library(lattice)
mydata4 <-read.table("/home/yiliao/Documents/Pepper_2021/ABcompartments/draw/FINAL/All/Final/ABvsRNA.bed.bed.RNA.dif.bed",head=T)
mydata4 
mydata4$RNAFOLDRank <- as.factor(mydata4$RNAFOLDRank)
mydata4$RankChange <- as.numeric(mydata4$RankChange)
mydata4$RankChange[mydata4$RankChange>0.2]<-0.25
mydata4$RankChange[mydata4$RankChange<(-0.2)]<-(-0.25)

mydata5 <-mydata4 %>% filter(RNAFOLDRank=='stable')
mydata6 <-mydata4 %>% filter(RNAFOLDRank=='up')
mydata7 <-mydata4 %>% filter(RNAFOLDRank=='down')
dim(mydata5)
dim(mydata6)
dim(mydata7)
mydata4
mydata8<-mydata4 %>% filter(Tissues=='GRvsHL')
mydata5 <-mydata8 %>% filter(RNAFOLDRank=='stable')
mydata6 <-mydata8 %>% filter(RNAFOLDRank=='up')
mydata7 <-mydata8 %>% filter(RNAFOLDRank=='down')
dim(mydata5)
dim(mydata6)
dim(mydata7)

ggplot(mydata8, aes(x= RankChange, group=RNAFOLDRank)) + 
  geom_bar(aes(y = ..prop.., fill =factor(..x..)), stat="count") +
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.5) +
  labs(y = "Percent", fill="RankChange") +
  facet_grid(~RNAFOLDRank) +
  scale_y_continuous(labels = scales::percent)

wilcox.test(RankChange~RNAFOLDRank, data = mydata5, paired = FALSE,
            alternative = "less")







####Figure 7a
mydata2 <-read.table("/home/yiliao/Documents/Pepper_2021/ABcompartments/draw/FINAL/All/20211015/GeneBased/NEW/One/ABvsRNA.one.gene.bed.bed", head=T)

mydata2 <-read.table("/home/yiliao/Documents/Pepper_2021/ABcompartments/draw/FINAL/All/20211015/GeneBased/NEW/two/ABvsRNA.one.gene.bed.bed", head=T)

mydata2 <-read.table("/home/yiliao/Documents/Pepper_2021/ABcompartments/draw/FINAL/All/20211015/BinBased/New/one/ABvsRNA.one.gene.bed.bed", head=T)

mydata2 <-read.table("/home/yiliao/Documents/Pepper_2021/ABcompartments/draw/FINAL/All/20211015/BinBased/New/two/ABvsRNA.one.gene.bed.bed", head=T)


head(mydata2)
mydata2$Tissues <- as.factor(mydata2$Tissues )
mydata9 <-mydata2 %>% filter(Pvalue < 0.01) 
ggplot(mydata9, aes(x=Tissues, y=RNAFold, fill=factor(Rank))) + geom_boxplot(outlier.shape = NA) + ylim(-15,15)+labs(fill = "Rank") +ylab("log2(fold change)") + xlab("Comparisons")

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



####Figure 7b
library(ggplot2)
library(plyr)
library(lattice)
mydata <-read.table("/home/yiliao/Documents/Pepper_2021/ABcompartments/draw/FINAL/All/Final/20211109/ABvsRNA.bed.bed.RNA.dif.bed",head=T)

mydata$RNAFOLDRank <- as.factor(mydata$RNAFOLDRank)
mydata$RankChange <- as.numeric(mydata$RankChange)

mydata$RankChange[mydata$RankChange>0.2]<-0.25
mydata$RankChange[mydata$RankChange<(-0.2)]<-(-0.25)

#mydata1 <-mydata %>% filter(RNAFOLDRank=='stable')
#mydata2 <-mydata %>% filter(RNAFOLDRank=='up')
#mydata3 <-mydata %>% filter(RNAFOLDRank=='down')
#mydata8<-mydata4 %>% filter(Tissues=='GRvsYP' | Tissues=='TZvsYP' | Tissues =='HLvsYP')
mydata4<-mydata %>% filter(Tissues=='TZvsGR')

mydata1 <-mydata4 %>% filter(RNAFOLDRank=='stable')
mydata2 <-mydata4 %>% filter(RNAFOLDRank=='up')
mydata3 <-mydata4 %>% filter(RNAFOLDRank=='down')

dim(mydata1)
dim(mydata2)
dim(mydata3)


ggplot(mydata4, aes(x= RankChange, group=RNAFOLDRank)) + 
  geom_bar(aes(y = ..prop.., fill =factor(..x..)), stat="count") +
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.5) +
  labs(y = "Percent", fill="RankChange") +
  facet_grid(~RNAFOLDRank) +
  scale_y_continuous(labels = scales::percent)




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

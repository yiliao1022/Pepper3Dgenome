install.packages("ggplot2")
install.packages("gplots")
install.packages("amap")
install.packages("RColorBrewer")
install.packages("dplyr")
install.packages("ggpubr")

library("ggplot2")
library("gplots")
library("amap")
library("RColorBrewer")
library("dplyr")
library("ggpubr")


#setwd("/home/yiliao/Documents/Pepper_2021/ABcompartments/draw/2022_revision/4papper/Rscript")
#Dat<-read.table("sample.lst.bed.bed.sort.txt.bed",header = TRUE)
#Dat<-read.table("ABvsRNA.40k.bed.bed.bed.txt.bed",header = TRUE)

### For 10-kb
setwd("/home/yiliao/igv/pepper/10k")
Dat<-read.table("sample.lst.bed.bed.txt.bed",header = TRUE)

### For 40-kb
setwd("/home/yiliao/igv/pepper/40k/Final")
Dat<-read.table("ABvsRNA.40k.bed.bed.final.bed.txt.bed",header = TRUE)


head(Dat)

#boxplot(CHH~YP, data = Dat, outline=FALSE)
#regline <- lm(CHH ~ as.numeric(YP), data = Dat)

#X = sort(unique(Dat$YP))
#lines(x=1:length(X),
 #     y=predict(regline,data.frame(YP=X)),col="blue",lty=8)

#cor.test(x=1:length(X),
 #        y=predict(regline,data.frame(YP=X)), method = "spearman")
#cor.test(x=Dat$YP,
 #        y=Dat$CHH, method = "spearman")


Means <-colMeans(Dat[sapply(Dat, is.numeric)],na.rm=TRUE) 


###AB
mydataA <-Dat %>% filter(AB=='A') 
mydataB <-Dat %>% filter(AB=='B')
MeansA <-colMeans(mydataA[sapply(mydataA, is.numeric)],na.rm=TRUE)  
MeansB <-colMeans(mydataB[sapply(mydataB, is.numeric)],na.rm=TRUE)  

GeneAB <- c(MeansA[['Gene']]/Means[['Gene']], 
           MeansB[['Gene']]/Means[['Gene']])
log2GeneAB <- sapply(GeneAB, function(x) log2(GeneAB))
log2GeneAB[,1]

LTRAB <- c(MeansA[['LTR']]/Means[['LTR']], 
          MeansB[['LTR']]/Means[['LTR']])
log2LTRAB <- sapply(LTRAB, function(x) log2(LTRAB))
log2LTRAB[,1]

ExpAB <- c(MeansA[['Expression']]/Means[['Expression']], 
          MeansB[['Expression']]/Means[['Expression']])
log2ExpAB <- sapply(ExpAB, function(x) log2(ExpAB))
log2ExpAB[,1]

Methy1AB <- c(MeansA[['Methy1']]/Means[['Methy1']], 
             MeansB[['Methy1']]/Means[['Methy1']])
log2Methy1AB <- sapply(Methy1AB, function(x) log2(Methy1AB))
log2Methy1AB[,1]


CG1AB <- c(MeansA[['CG']]/Means[['CG']], 
              MeansB[['CG']]/Means[['CG']])
log2CG1AB <- sapply(CG1AB, function(x) log2(CG1AB))
log2CG1AB[,1]

CHG1AB <- c(MeansA[['CHG']]/Means[['CHG']], 
              MeansB[['CHG']]/Means[['CHG']])
log2CHG1AB <- sapply(CHG1AB, function(x) log2(CHG1AB))
log2CHG1AB[,1]

CHH1AB <- c(MeansA[['CHH']]/Means[['CHH']], 
              MeansB[['CHH']]/Means[['CHH']])
log2CHH1AB <- sapply(CHH1AB, function(x) log2(CHH1AB))
log2CHH1AB[,1]


Methy2AB <- c(MeansA[['Methy2']]/Means[['Methy2']], 
             MeansB[['Methy2']]/Means[['Methy2']])
log2Methy2AB <- sapply(Methy2AB, function(x) log2(Methy2AB))
log2Methy2AB[,1]

H3K4_1AB <-c (MeansA[['H3K4_1']]/Means[['H3K4_1']],
             MeansB[['H3K4_1']]/Means[['H3K4_1']])
log2H3K4_1AB <- sapply(H3K4_1AB, function(x) log2(H3K4_1AB))
log2H3K4_1AB[,1]

H3K4_2AB <-c(MeansA[['H3K4_2']]/Means[['H3K4_2']],
            MeansB[['H3K4_2']]/Means[['H3K4_2']])
log2H3K4_2AB <- sapply(H3K4_2AB, function(x) log2(H3K4_2AB))
log2H3K4_2AB[,1]

H3K9_1AB <- c(MeansA[['H3K9_1']]/Means[['H3K9_1']],
             MeansB[['H3K9_1']]/Means[['H3K9_1']])
log2H3K9_1AB <- sapply(H3K9_1AB, function(x) log2(H3K9_1AB))
log2H3K9_1AB[,1]

H3K9_2AB <- c(MeansA[['H3K9_2']]/Means[['H3K9_2']],
             MeansB[['H3K9_2']]/Means[['H3K9_2']])
log2H3K9_2AB <- sapply(H3K9_2AB, function(x) log2(H3K9_2AB))
log2H3K9_2AB[,1]

H3K27_1AB <- c(MeansA[['H3K27_1']]/Means[['H3K27_1']],
              MeansB[['H3K27_1']]/Means[['H3K27_1']])
log2H3K27_1AB <- sapply(H3K27_1AB, function(x) log2(H3K27_1AB))
log2H3K27_1AB[,1]

H3K27_2AB <- c(MeansA[['H3K27_2']]/Means[['H3K27_2']],
              MeansB[['H3K27_2']]/Means[['H3K27_2']])
log2H3K27_2AB <- sapply(H3K27_2AB, function(x) log2(H3K27_2AB))
log2H3K27_2AB[,1]

##AABB
mydataA1 <-Dat %>% filter(AABB=='A1') 
mydataA2 <-Dat %>% filter(AABB=='A2')
mydataB1 <-Dat %>% filter(AABB=='B1') 
mydataB2 <-Dat %>% filter(AABB=='B2')

MeansA1 <-colMeans(mydataA1[sapply(mydataA1, is.numeric)],na.rm=TRUE)  
MeansA2 <-colMeans(mydataA2[sapply(mydataA2, is.numeric)],na.rm=TRUE) 
MeansB1 <-colMeans(mydataB1[sapply(mydataB1, is.numeric)],na.rm=TRUE)  
MeansB2 <-colMeans(mydataB2[sapply(mydataB2, is.numeric)],na.rm=TRUE)

Gene4 <- c(MeansA1[['Gene']]/Means[['Gene']], 
          MeansA2[['Gene']]/Means[['Gene']], 
          MeansB1[['Gene']]/Means[['Gene']], 
          MeansB2[['Gene']]/Means[['Gene']])
log2Gene4 <- sapply(Gene4, function(x) log2(Gene4))
log2Gene4[,1]

LTR4 <- c(MeansA1[['LTR']]/Means[['LTR']], 
           MeansA2[['LTR']]/Means[['LTR']], 
           MeansB1[['LTR']]/Means[['LTR']], 
           MeansB2[['LTR']]/Means[['LTR']])
log2LTR4 <- sapply(LTR4, function(x) log2(LTR4))
log2LTR4[,1]

Exp4 <- c(MeansA1[['Expression']]/Means[['Expression']], 
          MeansA2[['Expression']]/Means[['Expression']], 
          MeansB1[['Expression']]/Means[['Expression']], 
          MeansB2[['Expression']]/Means[['Expression']])
log2Exp4 <- sapply(Exp4, function(x) log2(Exp4))
log2Exp4[,1]

CG41 <- c(MeansA1[['CG']]/Means[['CG']], 
            MeansA2[['CG']]/Means[['CG']], 
            MeansB1[['CG']]/Means[['CG']], 
            MeansB2[['CG']]/Means[['CG']])
log2CG41 <- sapply(CG41, function(x) log2(CG41))
log2CG41[,1]


CHG41 <- c(MeansA1[['CHG']]/Means[['CHG']], 
          MeansA2[['CHG']]/Means[['CHG']], 
          MeansB1[['CHG']]/Means[['CHG']], 
          MeansB2[['CHG']]/Means[['CHG']])
log2CHG41 <- sapply(CHG41, function(x) log2(CHG41))
log2CHG41[,1]


CHH41 <- c(MeansA1[['CHH']]/Means[['CHH']], 
           MeansA2[['CHH']]/Means[['CHH']], 
           MeansB1[['CHH']]/Means[['CHH']], 
           MeansB2[['CHH']]/Means[['CHH']])
log2CHH41 <- sapply(CHH41, function(x) log2(CHH41))
log2CHH41[,1]

Methy41 <- c(MeansA1[['Methy1']]/Means[['Methy1']], 
             MeansA2[['Methy1']]/Means[['Methy1']], 
             MeansB1[['Methy1']]/Means[['Methy1']], 
             MeansB2[['Methy1']]/Means[['Methy1']])
log2Methy41 <- sapply(Methy41, function(x) log2(Methy41))
log2Methy41[,1]


Methy42 <- c(MeansA1[['Methy2']]/Means[['Methy2']], 
            MeansA2[['Methy2']]/Means[['Methy2']], 
            MeansB1[['Methy2']]/Means[['Methy2']], 
            MeansB2[['Methy2']]/Means[['Methy2']])
log2Methy42 <- sapply(Methy42, function(x) log2(Methy42))
log2Methy42[,1]

H3K4_41 <-c (MeansA1[['H3K4_1']]/Means[['H3K4_1']],
            MeansA2[['H3K4_1']]/Means[['H3K4_1']],
            MeansB1[['H3K4_1']]/Means[['H3K4_1']],
            MeansB2[['H3K4_1']]/Means[['H3K4_1']])
log2H3K4_41 <- sapply(H3K4_41, function(x) log2(H3K4_41))
log2H3K4_41[,1]

H3K4_42 <-c(MeansA1[['H3K4_2']]/Means[['H3K4_2']],
           MeansA2[['H3K4_2']]/Means[['H3K4_2']],
           MeansB1[['H3K4_2']]/Means[['H3K4_2']],
           MeansB2[['H3K4_2']]/Means[['H3K4_2']])
log2H3K4_42 <- sapply(H3K4_42, function(x) log2(H3K4_42))
log2H3K4_42[,1]

H3K9_41 <- c(MeansA1[['H3K9_1']]/Means[['H3K9_1']],
            MeansA2[['H3K9_1']]/Means[['H3K9_1']],
            MeansB1[['H3K9_1']]/Means[['H3K9_1']],
            MeansB2[['H3K9_1']]/Means[['H3K9_1']])
log2H3K9_41 <- sapply(H3K9_41, function(x) log2(H3K9_41))
log2H3K9_41[,1]

H3K9_42 <- c(MeansA1[['H3K9_2']]/Means[['H3K9_2']],
            MeansA2[['H3K9_2']]/Means[['H3K9_2']],
            MeansB1[['H3K9_2']]/Means[['H3K9_2']],
            MeansB2[['H3K9_2']]/Means[['H3K9_2']])
log2H3K9_42 <- sapply(H3K9_42, function(x) log2(H3K9_42))
log2H3K9_42[,1]

H3K27_41 <- c(MeansA1[['H3K27_1']]/Means[['H3K27_1']],
             MeansA2[['H3K27_1']]/Means[['H3K27_1']],
             MeansB1[['H3K27_1']]/Means[['H3K27_1']],
             MeansB2[['H3K27_1']]/Means[['H3K27_1']])
log2H3K27_41 <- sapply(H3K27_41, function(x) log2(H3K27_41))
log2H3K27_41[,1]

H3K27_42 <- c(MeansA1[['H3K27_2']]/Means[['H3K27_2']],
             MeansA2[['H3K27_2']]/Means[['H3K27_2']],
             MeansB1[['H3K27_2']]/Means[['H3K27_2']],
             MeansB2[['H3K27_2']]/Means[['H3K27_2']])
log2H3K27_42 <- sapply(H3K27_42, function(x) log2(H3K27_42))
log2H3K27_42[,1]


#AAAABBBB
mydataA1.1 <-Dat %>% filter(AAAABBBB=='A1.1') 
mydataA1.2 <-Dat %>% filter(AAAABBBB=='A1.2')
mydataA2.1 <-Dat %>% filter(AAAABBBB=='A2.1') 
mydataA2.2 <-Dat %>% filter(AAAABBBB=='A2.2')
mydataB1.1 <-Dat %>% filter(AAAABBBB=='B1.1') 
mydataB1.2 <-Dat %>% filter(AAAABBBB=='B1.2')
mydataB2.1 <-Dat %>% filter(AAAABBBB=='B2.1') 
mydataB2.2 <-Dat %>% filter(AAAABBBB=='B2.2')

MeansA1.1 <-colMeans(mydataA1.1[sapply(mydataA1.1, is.numeric)],na.rm=TRUE)  
MeansA1.2 <-colMeans(mydataA1.2[sapply(mydataA1.2, is.numeric)],na.rm=TRUE) 
MeansA2.1 <-colMeans(mydataA2.1[sapply(mydataA2.1, is.numeric)],na.rm=TRUE)  
MeansA2.2 <-colMeans(mydataA2.2[sapply(mydataA2.2, is.numeric)],na.rm=TRUE) 
MeansB1.1 <-colMeans(mydataB1.1[sapply(mydataB1.1, is.numeric)],na.rm=TRUE)  
MeansB1.2 <-colMeans(mydataB1.2[sapply(mydataB1.2, is.numeric)],na.rm=TRUE) 
MeansB2.1 <-colMeans(mydataB2.1[sapply(mydataB2.1, is.numeric)],na.rm=TRUE)  
MeansB2.2 <-colMeans(mydataB2.2[sapply(mydataB2.2, is.numeric)],na.rm=TRUE) 

Gene <- c(MeansA1.1[['Gene']]/Means[['Gene']], 
         MeansA1.2[['Gene']]/Means[['Gene']], 
         MeansA2.1[['Gene']]/Means[['Gene']], 
         MeansA2.2[['Gene']]/Means[['Gene']], 
         MeansB1.1[['Gene']]/Means[['Gene']], 
         MeansB1.2[['Gene']]/Means[['Gene']], 
         MeansB2.1[['Gene']]/Means[['Gene']], 
         MeansB2.2[['Gene']]/Means[['Gene']])

log2Gene <- sapply(Gene, function(x) log2(Gene))
log2Gene[,1]

LTR <- c(MeansA1.1[['LTR']]/Means[['LTR']], 
          MeansA1.2[['LTR']]/Means[['LTR']], 
          MeansA2.1[['LTR']]/Means[['LTR']], 
          MeansA2.2[['LTR']]/Means[['LTR']], 
          MeansB1.1[['LTR']]/Means[['LTR']], 
          MeansB1.2[['LTR']]/Means[['LTR']], 
          MeansB2.1[['LTR']]/Means[['LTR']], 
          MeansB2.2[['LTR']]/Means[['LTR']])
log2LTR <- sapply(LTR, function(x) log2(LTR))
log2LTR[,1]


Methy1 <- c(MeansA1.1[['Methy1']]/Means[['Methy1']], 
         MeansA1.2[['Methy1']]/Means[['Methy1']], 
         MeansA2.1[['Methy1']]/Means[['Methy1']], 
         MeansA2.2[['Methy1']]/Means[['Methy1']], 
         MeansB1.1[['Methy1']]/Means[['Methy1']], 
         MeansB1.2[['Methy1']]/Means[['Methy1']], 
         MeansB2.1[['Methy1']]/Means[['Methy1']], 
         MeansB2.2[['Methy1']]/Means[['Methy1']])
log2Methy1 <- sapply(Methy1, function(x) log2(Methy1))
log2Methy1[,1]

Methy2 <- c(MeansA1.1[['Methy2']]/Means[['Methy2']], 
            MeansA1.2[['Methy2']]/Means[['Methy2']], 
            MeansA2.1[['Methy2']]/Means[['Methy2']], 
            MeansA2.2[['Methy2']]/Means[['Methy2']], 
            MeansB1.1[['Methy2']]/Means[['Methy2']], 
            MeansB1.2[['Methy2']]/Means[['Methy2']], 
            MeansB2.1[['Methy2']]/Means[['Methy2']], 
            MeansB2.2[['Methy2']]/Means[['Methy2']])
log2Methy2 <- sapply(Methy2, function(x) log2(Methy2))
log2Methy2[,1]

CG <- c(MeansA1.1[['CG']]/Means[['CG']], 
            MeansA1.2[['CG']]/Means[['CG']], 
            MeansA2.1[['CG']]/Means[['CG']], 
            MeansA2.2[['CG']]/Means[['CG']], 
            MeansB1.1[['CG']]/Means[['CG']], 
            MeansB1.2[['CG']]/Means[['CG']], 
            MeansB2.1[['CG']]/Means[['CG']], 
            MeansB2.2[['CG']]/Means[['CG']])
log2CG <- sapply(CG, function(x) log2(CG))
log2CG[,1]

CHG <- c(MeansA1.1[['CHG']]/Means[['CHG']], 
        MeansA1.2[['CHG']]/Means[['CHG']], 
        MeansA2.1[['CHG']]/Means[['CHG']], 
        MeansA2.2[['CHG']]/Means[['CHG']], 
        MeansB1.1[['CHG']]/Means[['CHG']], 
        MeansB1.2[['CHG']]/Means[['CHG']], 
        MeansB2.1[['CHG']]/Means[['CHG']], 
        MeansB2.2[['CHG']]/Means[['CHG']])
log2CHG <- sapply(CG, function(x) log2(CHG))
log2CHG[,1]

CHH <- c(MeansA1.1[['CHH']]/Means[['CHH']], 
        MeansA1.2[['CHH']]/Means[['CHH']], 
        MeansA2.1[['CHH']]/Means[['CHH']], 
        MeansA2.2[['CHH']]/Means[['CHH']], 
        MeansB1.1[['CHH']]/Means[['CHH']], 
        MeansB1.2[['CHH']]/Means[['CHH']], 
        MeansB2.1[['CHH']]/Means[['CHH']], 
        MeansB2.2[['CHH']]/Means[['CHH']])
log2CHH <- sapply(CHH, function(x) log2(CHH))
log2CHH[,1]



Exp <- c(MeansA1.1[['Expression']]/Means[['Expression']], 
         MeansA1.2[['Expression']]/Means[['Expression']], 
         MeansA2.1[['Expression']]/Means[['Expression']], 
         MeansA2.2[['Expression']]/Means[['Expression']], 
         MeansB1.1[['Expression']]/Means[['Expression']], 
         MeansB1.2[['Expression']]/Means[['Expression']], 
         MeansB2.1[['Expression']]/Means[['Expression']], 
         MeansB2.2[['Expression']]/Means[['Expression']])
log2Exp <- sapply(Exp, function(x) log2(Exp))
log2Exp[,1]

H3K4_1 <-c (MeansA1.1[['H3K4_1']]/Means[['H3K4_1']],
           MeansA1.2[['H3K4_1']]/Means[['H3K4_1']],
           MeansA2.1[['H3K4_1']]/Means[['H3K4_1']],
           MeansA2.2[['H3K4_1']]/Means[['H3K4_1']],
           MeansB1.1[['H3K4_1']]/Means[['H3K4_1']],
           MeansB1.2[['H3K4_1']]/Means[['H3K4_1']],
           MeansB2.1[['H3K4_1']]/Means[['H3K4_1']],
           MeansB2.2[['H3K4_1']]/Means[['H3K4_1']])
log2H3K4_1 <- sapply(H3K4_1, function(x) log2(H3K4_1))
log2H3K4_1[,1]

H3K4_2 <-c(MeansA1.1[['H3K4_2']]/Means[['H3K4_2']],
          MeansA1.2[['H3K4_2']]/Means[['H3K4_2']],
          MeansA2.1[['H3K4_2']]/Means[['H3K4_2']],
          MeansA2.2[['H3K4_2']]/Means[['H3K4_2']],
          MeansB1.1[['H3K4_2']]/Means[['H3K4_2']],
          MeansB1.2[['H3K4_2']]/Means[['H3K4_2']],
          MeansB2.1[['H3K4_2']]/Means[['H3K4_2']],
          MeansB2.2[['H3K4_2']]/Means[['H3K4_2']])
log2H3K4_2 <- sapply(H3K4_2, function(x) log2(H3K4_2))
log2H3K4_2[,1]

H3K9_1 <- c(MeansA1.1[['H3K9_1']]/Means[['H3K9_1']],
           MeansA1.2[['H3K9_1']]/Means[['H3K9_1']],
           MeansA2.1[['H3K9_1']]/Means[['H3K9_1']],
           MeansA2.2[['H3K9_1']]/Means[['H3K9_1']],
           MeansB1.1[['H3K9_1']]/Means[['H3K9_1']],
           MeansB1.2[['H3K9_1']]/Means[['H3K9_1']],
           MeansB2.1[['H3K9_1']]/Means[['H3K9_1']],
           MeansB2.2[['H3K9_1']]/Means[['H3K9_1']])
log2H3K9_1 <- sapply(H3K9_1, function(x) log2(H3K9_1))
log2H3K9_1[,1]

H3K9_2 <- c(MeansA1.1[['H3K9_2']]/Means[['H3K9_2']],
           MeansA1.2[['H3K9_2']]/Means[['H3K9_2']],
           MeansA2.1[['H3K9_2']]/Means[['H3K9_2']],
           MeansA2.2[['H3K9_2']]/Means[['H3K9_2']],
           MeansB1.1[['H3K9_2']]/Means[['H3K9_2']],
           MeansB1.2[['H3K9_2']]/Means[['H3K9_2']],
           MeansB2.1[['H3K9_2']]/Means[['H3K9_2']],
           MeansB2.2[['H3K9_2']]/Means[['H3K9_2']])
log2H3K9_2 <- sapply(H3K9_2, function(x) log2(H3K9_2))
log2H3K9_2[,1]

H3K27_1 <- c(MeansA1.1[['H3K27_1']]/Means[['H3K27_1']],
            MeansA1.2[['H3K27_1']]/Means[['H3K27_1']],
            MeansA2.1[['H3K27_1']]/Means[['H3K27_1']],
            MeansA2.2[['H3K27_1']]/Means[['H3K27_1']],
            MeansB1.1[['H3K27_1']]/Means[['H3K27_1']],
            MeansB1.2[['H3K27_1']]/Means[['H3K27_1']],
            MeansB2.1[['H3K27_1']]/Means[['H3K27_1']],
            MeansB2.2[['H3K27_1']]/Means[['H3K27_1']])
log2H3K27_1 <- sapply(H3K27_1, function(x) log2(H3K27_1))
log2H3K27_1[,1]

H3K27_2 <- c(MeansA1.1[['H3K27_2']]/Means[['H3K27_2']],
            MeansA1.2[['H3K27_2']]/Means[['H3K27_2']],
            MeansA2.1[['H3K27_2']]/Means[['H3K27_2']],
            MeansA2.2[['H3K27_2']]/Means[['H3K27_2']],
            MeansB1.1[['H3K27_2']]/Means[['H3K27_2']],
            MeansB1.2[['H3K27_2']]/Means[['H3K27_2']],
            MeansB2.1[['H3K27_2']]/Means[['H3K27_2']],
            MeansB2.2[['H3K27_2']]/Means[['H3K27_2']])
log2H3K27_2 <- sapply(H3K27_2, function(x) log2(H3K27_2))
log2H3K27_2[,1]

#######################AB
Name <- c("A","B")
DatMatAB <- data.frame (Name,log2GeneAB[,1],log2LTRAB[,1],log2Methy1AB[,1],log2Methy2AB[,1],log2CG1AB[,1],log2CHG1AB[,1],log2CHH1AB[,1],log2ExpAB[,1],log2H3K4_1AB[,1],log2H3K4_2AB[,1],log2H3K27_1AB[,1],log2H3K27_2AB[,1],log2H3K9_1AB[,1],log2H3K9_2AB[,1])
rnames <- DatMatAB[,1]
rnames
mat_dataAB <- data.matrix(DatMatAB[,2:ncol(DatMatAB)])
mat_dataAB
rownames(mat_dataAB) <- rnames
matdata_TAB <- t(mat_dataAB)
matdata_TAB

#######################AABB
Name <- c("A1", "A2", "B1","B2")
DatMat4 <- data.frame (Name,log2Gene4[,1],log2LTR4[,1],log2Methy41[,1],log2Methy42[,1],log2CG41[,1],log2CHG41[,1],log2CHH41[,1],log2Exp4[,1],log2H3K4_41[,1],log2H3K4_42[,1],log2H3K27_41[,1],log2H3K27_42[,1],log2H3K9_41[,1],log2H3K9_42[,1])
rnames <- DatMat4[,1]
rnames
mat_data4 <- data.matrix(DatMat4[,2:ncol(DatMat4)])
mat_data4
rownames(mat_data4) <- rnames
matdata_T4 <- t(mat_data4)
matdata_T4

##########################AAAABBBB
Name <- c("A1.1", "A1.2", "A2.1", "A2.2", "B1.1","B1.2","B2.1","B2.2")
DatMat <- data.frame (Name,log2Gene[,1],log2LTR[,1],log2Methy1[,1],log2Methy2[,1],log2CG[,1],log2CHG[,1],log2CHH[,1],log2Exp[,1],log2H3K4_1[,1],log2H3K4_2[,1],log2H3K27_1[,1],log2H3K27_2[,1],log2H3K9_1[,1],log2H3K9_2[,1])
rnames <- DatMat[,1]
rnames
mat_data <- data.matrix(DatMat[,2:ncol(DatMat)])
mat_data
rownames(mat_data) <- rnames
matdata_T <- t(mat_data)
matdata_T

###########################################
my_palette <- colorRampPalette(c("blue", "yellow", "red"))(n = 299)
col_breaks = c(seq(-2,-0.51,length=100), # for red
               seq(-0.5,0.5,length=100),  # for yellow
               seq(0.51,2,length=100)) # for green
########################################################
#dev.off()
#par(mfrow=c(1,3))
#par(mfrow= c(3,1) )

heatmap.2(matdata_TAB,
          cellnote = matdata_TAB,  # same data set for cell labels
          main = "Correlation", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(10,10),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="column",     # only draw a row dendrogram
          #Colv="NA",
          key.par=list(mar=c(3.5,0,3,0))
) 


heatmap.2(matdata_T,
          cellnote = matdata_T,  # same data set for cell labels
          main = "Correlation", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(10,10),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="column",     # only draw a row dendrogram
          #Colv="NA",
          key.par=list(mar=c(3.5,0,3,0))
          ) 

heatmap.2(matdata_T4,
          cellnote = matdata_T4,  # same data set for cell labels
          main = "Correlation", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(10,10),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="column",     # only draw a row dendrogram
          #Colv="NA",
          key.par=list(mar=c(3.5,0,3,0))
) 

##########################
dev.copy2pdf(file="/home/yiliao/igv/pepper/40k/Final/AB2file.pdf")





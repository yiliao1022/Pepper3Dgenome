library(ggplot2)
install.packages("ggpubr")
library(ggpubr)

setwd("~/Documents/Pepper_2021/Final/5_Figure5_TAD/Figure5_b/2022_revision")

############### Features within TADs
Sdata<-read.csv(file = "matrixYP2.all.domain.raw.Cluster.new.bed.bed.bed.bed",header=TRUE,sep="\t")
head(Sdata)
Sdata["Group"][Sdata["Group"] == 1] <- "Active"
Sdata["Group"][Sdata["Group"] == 2] <- "Inactive"
Sdata["Group"][Sdata["Group"] == 3] <- "HP1"

Sdata$group <-as.factor(Sdata$Group)
Sdata$Size <-as.numeric(Sdata$Size)


p <- ggboxplot(Sdata, x = "Group", y = "Methy1", ylim=c(0.25,0.75),
               color = "Group", palette =c("#00AFBB", "#E7B800", "#FC4E07"), notch = TRUE,
               shape = "Group", xlab = "TAD group", ylab = "Total methylation level", title = "DNA methylation",  names= c("Active","Inactive","HP1") )

my_comparisons <- list( c("Active", "Inactive"), c("Active", "HP1"), c("Inactive", "HP1") )
p + stat_compare_means(comparisons = my_comparisons, method =  "wilcox.test") + # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50) 



p<-boxplot (Methy1~Group, data=Sdata, xlab="Domain group", ylab="Total methylation level", main="DNA Methylation", notch= TRUE, varwidth= TRUE, col=c("green","yellow","purple"), names=c("Active","Inactive","HP1"))




methy <-ggplot(Sdata, aes(x=Group,y=Methy1,fill=group)  ) + ylim (0.25,0.75) + geom_boxplot(outlier.colour="black", outlier.shape=16,
                                                                                          outlier.size=1)  +
  scale_fill_manual(values=c("#008080","#E69F00","#8B4513")) + labs(size= "",
                                                                    x = "",
                                                                    y = "Intensity",
                                                                    title = "DNA Methylation")

H3K9<-ggplot(Sdata, aes(x=Group,y=H3K9_1,fill=group)) + ylim (0,10) + geom_boxplot(outlier.colour="black", outlier.shape=16,
                                                                                   outlier.size=1)   +
  scale_fill_manual(values=c("#008080","#E69F00","#8B4513")) + labs(size= "",
                                                                  x = "",
                                                                  y = "Intensity",
                                                                  title = "H3K9me2")

H3K27<-ggplot(Sdata, aes(x=Group,y=H3K27_1,fill=group)) + ylim (0,10) + geom_boxplot(outlier.colour="black", outlier.shape=16,
                                                                                     outlier.size=1)   +
  scale_fill_manual(values=c("#008080","#E69F00","#8B4513")) + labs(size= "",
                                                                    x = "",
                                                                    y = "Intensity",
                                                                    title = "H3K27me3")

H3K4<-ggplot(Sdata, aes(x=Group,y=H3K4_1,fill=group)) + ylim (0,10) + geom_boxplot(outlier.colour="black", outlier.shape=16,
                                                                                   outlier.size=1)   +
  scale_fill_manual(values=c("#008080","#E69F00","#8B4513")) + labs(size= "",
                                                                      x = "",
                                                                      y = "Intensity",
                                                                      title = "H3K4me3")

LTR<-ggplot(Sdata, aes(x=Group,y=LTR,fill=group)) + ylim (0,1) + geom_boxplot(outlier.colour="black", outlier.shape=16,
                                                                              outlier.size=1)   +
  scale_fill_manual(values=c("#008080","#E69F00","#8B4513")) +  labs(size= "",
                                                                      x = "",
                                                                      y = "Intensity",
                                                                      title = "LTR")

Gene<-ggplot(Sdata, aes(x=Group,y=Gene,fill=group)) + ylim (0,1) + geom_boxplot(outlier.colour="black", outlier.shape=16,
                                                                                outlier.size=1)   +
  scale_fill_manual(values=c("#008080","#E69F00","#8B4513")) + labs(size= "",
                                                                    x = "",
                                                                    y = "Intensity",
                                                                    title = "Gene")

Size <-ggplot(Sdata, aes(x=Group,y=Size,fill=group)) + ylim (100000,4000000) + geom_boxplot(outlier.colour="black", outlier.shape=16,
                                                                                outlier.size=1)   +
  scale_fill_manual(values=c("#008080","#E69F00","#8B4513")) + labs (
                                                                    x = "",
                                                                    y = "Intensity",
                                                                    title = "Size")

ggarrange(Gene, LTR, methy, H3K9, H3K4, H3K27, Size,
          labels = c("A", "B","C","D","E","F","G"),
          ncol = 4, nrow = 2)

############### Features around TAD boundaries

setwd("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/Figure5_b/2022_revision/Boundary")

## For gene
deletion1 <-read.table("A2A.bed.Gene.bedgraph.peak.cov.bed.bed.bed")
deletion2 <-read.table("A2I.bed.Gene.bedgraph.peak.cov.bed.bed.bed")
deletion3 <-read.table("A2H.bed.Gene.bedgraph.peak.cov.bed.bed.bed")
deletion4 <-read.table("I2I.bed.Gene.bedgraph.peak.cov.bed.bed.bed")
deletion5 <-read.table("I2H.bed.Gene.bedgraph.peak.cov.bed.bed.bed")
deletion6 <-read.table("H2H.bed.Gene.bedgraph.peak.cov.bed.bed.bed")

colnames(deletion1)<-c("label","counts")
deletion1 <- cbind(deletion1, type="active_active")
deletion1$counts <- log2(deletion1$counts/mean(deletion1$counts))

colnames(deletion2) <- c("label","counts")
deletion2 <- cbind(deletion2 , type = "active_inactive")
deletion2$counts <- log2(deletion2$counts/mean(deletion2$counts))

colnames(deletion3)<-c("label","counts")
deletion3 <- cbind(deletion3, type="active_medium")
deletion3$counts <- log2(deletion3$counts/mean(deletion3$counts))

colnames(deletion4)<-c("label","counts")
deletion4 <- cbind(deletion4, type="inactive_inactive")
deletion4$counts <- log2(deletion4$counts/mean(deletion4$counts))

colnames(deletion5)<-c("label","counts")
deletion5 <- cbind(deletion5, type="medium_inactive")
deletion5$counts <- log2(deletion5$counts/mean(deletion5$counts))

colnames(deletion6)<-c("label","counts")
deletion6 <- cbind(deletion6, type="medium_medium")
deletion6$counts <- log2(deletion6$counts/mean(deletion6$counts))


mydf1 <- rbind(deletion1,deletion2,deletion3,deletion4,deletion5,deletion6)
mydf1
Gene <- ggplot(data=mydf1,aes(x=label, y = counts, color=type))  + ggtitle("Gene") +
  geom_smooth(method="loess",span=0.2,se=T) + 
  scale_color_manual(values=c("#CD5C5C","#FFA07A","#FFD700","black","#008B8B","#4169E1")) + 
  labs(color="Types",x="Distance to TAD boundary",y="log2(Observed/average)") #+
  #theme( legend.position=c(0.8,0.1),legend.background = element_blank(), legend.title= element_text(size=11, face="bold", color='blue'))

## For H3K9
deletion1 <-read.table("A2A.bed.H3K9.bedgraph.peak.cov.bed.bed.bed")
deletion2 <-read.table("A2I.bed.H3K9.bedgraph.peak.cov.bed.bed.bed")
deletion3 <-read.table("A2H.bed.H3K9.bedgraph.peak.cov.bed.bed.bed")
deletion4 <-read.table("I2I.bed.H3K9.bedgraph.peak.cov.bed.bed.bed")
deletion5 <-read.table("I2H.bed.H3K9.bedgraph.peak.cov.bed.bed.bed")
deletion6 <-read.table("H2H.bed.H3K9.bedgraph.peak.cov.bed.bed.bed")

colnames(deletion1)<-c("label","counts")
deletion1 <- cbind(deletion1, type="active_active")
deletion1$counts <- log2(deletion1$counts/mean(deletion1$counts))

colnames(deletion2) <- c("label","counts")
deletion2 <- cbind(deletion2 , type = "active_inactive")
deletion2$counts <- log2(deletion2$counts/mean(deletion2$counts))

colnames(deletion3)<-c("label","counts")
deletion3 <- cbind(deletion3, type="active_medium")
deletion3$counts <- log2(deletion3$counts/mean(deletion3$counts))

colnames(deletion4)<-c("label","counts")
deletion4 <- cbind(deletion4, type="inactive_inactive")
deletion4$counts <- log2(deletion4$counts/mean(deletion4$counts))

colnames(deletion5)<-c("label","counts")
deletion5 <- cbind(deletion5, type="medium_inactive")
deletion5$counts <- log2(deletion5$counts/mean(deletion5$counts))

colnames(deletion6)<-c("label","counts")
deletion6 <- cbind(deletion6, type="medium_medium")
deletion6$counts <- log2(deletion6$counts/mean(deletion6$counts))


mydf1 <- rbind(deletion1,deletion2,deletion3,deletion4,deletion5,deletion6)
mydf1
H3K9<- ggplot(data=mydf1,aes(x=label, y = counts, color=type))  + ggtitle("H3K9me2") +
  geom_smooth(method="loess",span=0.2,se=T) + 
  scale_color_manual(values=c("#CD5C5C","#FFA07A","#FFD700","black","#008B8B","#4169E1")) + 
  labs(color="Types",x="Distance to TAD boundary",y="log2(Observed/average)") #+
 # theme( legend.position=c(0.8,0.1),legend.background = element_blank(), legend.title= element_text(size=11, face="bold", color='blue'))

## For H3K27
deletion1 <-read.table("A2A.bed.H3K27.bedgraph.peak.cov.bed.bed.bed")
deletion2 <-read.table("A2I.bed.H3K27.bedgraph.peak.cov.bed.bed.bed")
deletion3 <-read.table("A2H.bed.H3K27.bedgraph.peak.cov.bed.bed.bed")
deletion4 <-read.table("I2I.bed.H3K27.bedgraph.peak.cov.bed.bed.bed")
deletion5 <-read.table("I2H.bed.H3K27.bedgraph.peak.cov.bed.bed.bed")
deletion6 <-read.table("H2H.bed.H3K27.bedgraph.peak.cov.bed.bed.bed")

colnames(deletion1)<-c("label","counts")
deletion1 <- cbind(deletion1, type="active_active")
deletion1$counts <- log2(deletion1$counts/mean(deletion1$counts))

colnames(deletion2) <- c("label","counts")
deletion2 <- cbind(deletion2 , type = "active_inactive")
deletion2$counts <- log2(deletion2$counts/mean(deletion2$counts))

colnames(deletion3)<-c("label","counts")
deletion3 <- cbind(deletion3, type="active_medium")
deletion3$counts <- log2(deletion3$counts/mean(deletion3$counts))

colnames(deletion4)<-c("label","counts")
deletion4 <- cbind(deletion4, type="inactive_inactive")
deletion4$counts <- log2(deletion4$counts/mean(deletion4$counts))

colnames(deletion5)<-c("label","counts")
deletion5 <- cbind(deletion5, type="medium_inactive")
deletion5$counts <- log2(deletion5$counts/mean(deletion5$counts))

colnames(deletion6)<-c("label","counts")
deletion6 <- cbind(deletion6, type="medium_medium")
deletion6$counts <- log2(deletion6$counts/mean(deletion6$counts))


mydf1 <- rbind(deletion1,deletion2,deletion3,deletion4,deletion5,deletion6)
mydf1
H3K27<- ggplot(data=mydf1,aes(x=label, y = counts, color=type))   + ggtitle("H3K27me3") +
  geom_smooth(method="loess",span=0.2,se=T) + 
  scale_color_manual(values=c("#CD5C5C","#FFA07A","#FFD700","black","#008B8B","#4169E1")) + 
  labs(color="Types",x="Distance to TAD boundary",y="log2(Observed/average)") #+
  #theme( legend.position=c(0.8,0.1),legend.background = element_blank(), legend.title= element_text(size=11, face="bold", color='blue'))


## For H3K4
deletion1 <-read.table("A2A.bed.H3K4.bedgraph.peak.cov.bed.bed.bed")
deletion2 <-read.table("A2I.bed.H3K4.bedgraph.peak.cov.bed.bed.bed")
deletion3 <-read.table("A2H.bed.H3K4.bedgraph.peak.cov.bed.bed.bed")
deletion4 <-read.table("I2I.bed.H3K4.bedgraph.peak.cov.bed.bed.bed")
deletion5 <-read.table("I2H.bed.H3K4.bedgraph.peak.cov.bed.bed.bed")
deletion6 <-read.table("H2H.bed.H3K4.bedgraph.peak.cov.bed.bed.bed")

colnames(deletion1)<-c("label","counts")
deletion1 <- cbind(deletion1, type="active_active")
deletion1$counts <- log2(deletion1$counts/mean(deletion1$counts))

colnames(deletion2) <- c("label","counts")
deletion2 <- cbind(deletion2 , type = "active_inactive")
deletion2$counts <- log2(deletion2$counts/mean(deletion2$counts))

colnames(deletion3)<-c("label","counts")
deletion3 <- cbind(deletion3, type="active_medium")
deletion3$counts <- log2(deletion3$counts/mean(deletion3$counts))

colnames(deletion4)<-c("label","counts")
deletion4 <- cbind(deletion4, type="inactive_inactive")
deletion4$counts <- log2(deletion4$counts/mean(deletion4$counts))

colnames(deletion5)<-c("label","counts")
deletion5 <- cbind(deletion5, type="medium_inactive")
deletion5$counts <- log2(deletion5$counts/mean(deletion5$counts))

colnames(deletion6)<-c("label","counts")
deletion6 <- cbind(deletion6, type="medium_medium")
deletion6$counts <- log2(deletion6$counts/mean(deletion6$counts))


mydf1 <- rbind(deletion1,deletion2,deletion3,deletion4,deletion5,deletion6)
mydf1
H3K4<- ggplot(data=mydf1,aes(x=label, y = counts, color=type))  +  ggtitle("H3K4me3") +
  geom_smooth(method="loess",span=0.2,se=T) + 
  scale_color_manual(values=c("#CD5C5C","#FFA07A","#FFD700","#556B2F","#008B8B","#4169E1")) + 
  labs(color="Types",x="Distance to TAD boundary",y="log2(Observed/average)") #+
  #theme( legend.position=c(0.8,0.1),legend.background = element_blank(), legend.title= element_text(size=11, face="bold", color='blue'))

## For methylation
deletion1 <-read.table("A2A.bed.Methy.bedgraph.peak.cov.bed.bed.bed")
deletion2 <-read.table("A2I.bed.Methy.bedgraph.peak.cov.bed.bed.bed")
deletion3 <-read.table("A2H.bed.Methy.bedgraph.peak.cov.bed.bed.bed")
deletion4 <-read.table("I2I.bed.Methy.bedgraph.peak.cov.bed.bed.bed")
deletion5 <-read.table("I2H.bed.Methy.bedgraph.peak.cov.bed.bed.bed")
deletion6 <-read.table("H2H.bed.Methy.bedgraph.peak.cov.bed.bed.bed")

colnames(deletion1)<-c("label","counts")
deletion1 <- cbind(deletion1, type="active_active")
deletion1$counts <- log2(deletion1$counts/mean(deletion1$counts))

colnames(deletion2) <- c("label","counts")
deletion2 <- cbind(deletion2 , type = "active_inactive")
deletion2$counts <- log2(deletion2$counts/mean(deletion2$counts))

colnames(deletion3)<-c("label","counts")
deletion3 <- cbind(deletion3, type="active_medium")
deletion3$counts <- log2(deletion3$counts/mean(deletion3$counts))

colnames(deletion4)<-c("label","counts")
deletion4 <- cbind(deletion4, type="inactive_inactive")
deletion4$counts <- log2(deletion4$counts/mean(deletion4$counts))

colnames(deletion5)<-c("label","counts")
deletion5 <- cbind(deletion5, type="medium_inactive")
deletion5$counts <- log2(deletion5$counts/mean(deletion5$counts))

colnames(deletion6)<-c("label","counts")
deletion6 <- cbind(deletion6, type="medium_medium")
deletion6$counts <- log2(deletion6$counts/mean(deletion6$counts))


mydf1 <- rbind(deletion1,deletion2,deletion3,deletion4,deletion5,deletion6)
mydf1
Methy<- ggplot(data=mydf1,aes(x=label, y = counts, color=type))  + ggtitle("DNA methylation") +
  geom_smooth(method="loess",span=0.2,se=T) + 
  scale_color_manual(values=c("#CD5C5C","#FFA07A","#FFD700","#556B2F","#008B8B","#4169E1")) + 
  labs(color="Types",x="Distance to TAD boundary",y="log2(Observed/average)") #+
 # theme( legend.position=c(0.8,0.1),legend.background = element_blank(), legend.title= element_text(size=11, face="bold", color='blue'))

## For LTR
deletion1 <-read.table("A2A.bed.LTR.bedgraph.peak.cov.bed.bed.bed")
deletion2 <-read.table("A2I.bed.LTR.bedgraph.peak.cov.bed.bed.bed")
deletion3 <-read.table("A2H.bed.LTR.bedgraph.peak.cov.bed.bed.bed")
deletion4 <-read.table("I2I.bed.LTR.bedgraph.peak.cov.bed.bed.bed")
deletion5 <-read.table("I2H.bed.LTR.bedgraph.peak.cov.bed.bed.bed")
deletion6 <-read.table("H2H.bed.LTR.bedgraph.peak.cov.bed.bed.bed")

colnames(deletion1)<-c("label","counts")
deletion1 <- cbind(deletion1, type="active_active")
deletion1$counts <- log2(deletion1$counts/mean(deletion1$counts))

colnames(deletion2) <- c("label","counts")
deletion2 <- cbind(deletion2 , type = "active_inactive")
deletion2$counts <- log2(deletion2$counts/mean(deletion2$counts))

colnames(deletion3)<-c("label","counts")
deletion3 <- cbind(deletion3, type="active_medium")
deletion3$counts <- log2(deletion3$counts/mean(deletion3$counts))

colnames(deletion4)<-c("label","counts")
deletion4 <- cbind(deletion4, type="inactive_inactive")
deletion4$counts <- log2(deletion4$counts/mean(deletion4$counts))

colnames(deletion5)<-c("label","counts")
deletion5 <- cbind(deletion5, type="medium_inactive")
deletion5$counts <- log2(deletion5$counts/mean(deletion5$counts))

colnames(deletion6)<-c("label","counts")
deletion6 <- cbind(deletion6, type="medium_medium")
deletion6$counts <- log2(deletion6$counts/mean(deletion6$counts))


mydf1 <- rbind(deletion1,deletion2,deletion3,deletion4,deletion5,deletion6)
mydf1
LTR<- ggplot(data=mydf1,aes(x=label, y = counts, color=type))  +  ggtitle("LTR") +
  geom_smooth(method="loess",span=0.2,se=T) + 
  scale_color_manual(values=c("#CD5C5C","#FFA07A","#FFD700","#556B2F","#008B8B","#4169E1")) + 
  labs(color="Types",x="Distance to TAD boundary",y="log2(Observed/average)") #+
 # theme( legend.position=c(0.8,0.1),legend.background = element_blank(), legend.title= element_text(size=11, face="bold", color='blue'))

ggarrange(Gene, LTR, Methy, H3K9, H3K4, H3K27,
          labels = c("A", "B","C","D","E","F"),
          ncol = 3, nrow = 2)


library(ggplot2)
setwd("~/Documents/Pepper_2021/Figure6/Peak")

### For close coverage 

deletion1 <-read.table("CA59.Bacca.sing.sorted.bed.peak.cov.bed.bed")
deletion2 <-read.table("CA59.Chine.sing.sorted.bed.peak.cov.bed.bed")
deletion3 <-read.table("CA59.CM334.sing.sorted.bed.peak.cov.bed.bed")
deletion4 <-read.table("CA59.Glabr.sing.sorted.bed.peak.cov.bed.bed")
deletion5 <-read.table("CA59.Zunla.sing.sorted.bed.peak.cov.bed.bed")

colnames(deletion1)<-c("label","counts")
deletion1 <- cbind(deletion1, type="Bacca")

colnames(deletion2)<-c("label","counts")
deletion2 <- cbind(deletion2, type="Chine")

colnames(deletion3)<-c("label","counts")
deletion3 <- cbind(deletion3, type="CM334")

colnames(deletion4)<-c("label","counts")
deletion4 <- cbind(deletion4, type="Glabr")

colnames(deletion5)<-c("label","counts")
deletion5 <- cbind(deletion5, type="Zunla")

mydf1 <- rbind(deletion1,deletion2,deletion3,deletion4,deletion5)
ggplot(data=mydf1,aes(x=label, y = counts, color=type))  + 
  geom_smooth(method="loess",span=0.05,se=F) + 
  scale_color_manual(values=c("blue","Orange","purple","green","black")) + 
  labs(color="Types",x="Distance to TAD boundary",y="Coverage") +
  theme( legend.position=c(0.8,0.1),legend.background = element_blank(), legend.title= element_text(size=11, face="bold", color='blue'))


### For distantly coverage 

deletion1 <-read.table("CA59.HQ.sing.sorted.bed.peak.cov.bed.bed")
deletion2 <-read.table("CA59.RH89A.sing.sorted.bed.peak.cov.bed.bed")
deletion3 <-read.table("CA59.RH89B.sing.sorted.bed.peak.cov.bed.bed")
deletion4 <-read.table("CA59.SL4.sing.sorted.bed.peak.cov.bed.bed")

colnames(deletion1)<-c("label","counts")
deletion1 <- cbind(deletion1, type="Eggplant")

colnames(deletion2)<-c("label","counts")
deletion2 <- cbind(deletion2, type="PotatoA")

colnames(deletion3)<-c("label","counts")
deletion3 <- cbind(deletion3, type="PotatoB")

colnames(deletion4)<-c("label","counts")
deletion4 <- cbind(deletion4, type="Tomato")

mydf1 <- rbind(deletion1,deletion2,deletion3,deletion4)
ggplot(data=mydf1,aes(x=label, y = counts, color=type))  + 
  geom_smooth(method="loess",span=0.05,se=F) + 
  scale_color_manual(values=c("blue","Orange","purple","green")) + 
  labs(color="Types",x="Distance to TAD boundary",y="Coverage") +
  theme( legend.position=c(0.8,0.1),legend.background = element_blank(), legend.title= element_text(size=11, face="bold", color='blue'))


library(ggplot2)
setwd("/home/yiliao/Documents/Pepper_2021/Figure6/Peak/tomato/Pepper_SNP_DEL/HiCExplorer")

### For close coverage 
deletion1 <-read.table("Zunla.bed.bed.bed")
deletion2 <-read.table("CM334.bed.bed.bed")
deletion3 <-read.table("Glabr.bed.bed.bed")
deletion4 <-read.table("Chine.bed.bed.bed")
deletion5 <-read.table("Bacca.bed.bed.bed")

colnames(deletion1)<-c("label","counts","missing","SNP","DEL")
deletion1 <- cbind(deletion1, type="Zunla")
deletion1$SNP<-deletion1$SNP*201/sum(deletion1$SNP)
deletion1$DEL<-deletion1$DEL*201/sum(deletion1$DEL)

colnames(deletion2)<-c("label","counts","missing","SNP","DEL")
deletion2 <- cbind(deletion2, type="CM334")
deletion2$SNP<-deletion2$SNP*201/sum(deletion2$SNP)
deletion2$DEL<-deletion2$DEL*201/sum(deletion2$DEL)

colnames(deletion3)<-c("label","counts","missing","SNP","DEL")
deletion3 <- cbind(deletion3, type="Glabr")
deletion3$SNP<-deletion3$SNP*201/sum(deletion3$SNP)
deletion3$DEL<-deletion3$DEL*201/sum(deletion3$DEL)

colnames(deletion4)<-c("label","counts","missing","SNP","DEL")
deletion4 <- cbind(deletion4, type="Chine")
deletion4$SNP<-deletion4$SNP*201/sum(deletion4$SNP)
deletion4$DEL<-deletion4$DEL*201/sum(deletion4$DEL)

colnames(deletion5)<-c("label","counts","missing","SNP","DEL")
deletion5 <- cbind(deletion5, type="Bacca")
deletion5$SNP<-deletion5$SNP*201/sum(deletion5$SNP)
deletion5$DEL<-deletion5$DEL*201/sum(deletion5$DEL)

mydf1 <- rbind(deletion1,deletion2,deletion3,deletion4,deletion5)
ggplot(data=mydf1,aes(x=label, y = DEL, color=type))  + 
  geom_smooth(method="loess",span=0.05,se=F) + 
  scale_color_manual(values=c("#E69F00","#56B4E9","#009E73","#CC79A7","#0072B2")) + 
  labs(color="Types",x="Distance to TAD boundary",y="Obs/Exp deletion coverage") +
  theme( legend.position=c(0.8,0.1),legend.background = element_blank(), legend.title= element_text(size=11, face="bold", color='blue'))

ggplot(data=mydf1,aes(x=label, y = SNP, color=type))  + 
  geom_smooth(method="loess",span=0.05,se=F) + 
  scale_color_manual(values=c("#E69F00","#56B4E9","#009E73","#CC79A7","#0072B2")) + 
  labs(color="Types",x="Distance to TAD boundary",y="Obs/Exp SNP density") +
  theme( legend.position=c(0.8,0.1),legend.background = element_blank(), legend.title= element_text(size=11, face="bold", color='blue'))

###########################################################################################

library(ggplot2)
setwd("/home/yiliao/Documents/Pepper_2021/Figure6/Peak/tomato/Pepper_SNP_DEL/TopDom/YP2")
### For close coverage 
deletion1 <-read.table("Zunla.bed.bed.bed")
deletion2 <-read.table("CM334.bed.bed.bed")
deletion3 <-read.table("Glabr.bed.bed.bed")
deletion4 <-read.table("Chine.bed.bed.bed")
deletion5 <-read.table("Bacca.bed.bed.bed")

colnames(deletion1)<-c("label","counts","missing","SNP","DEL")
deletion1 <- cbind(deletion1, type="Zunla")
deletion1$SNP<-deletion1$SNP*201/sum(deletion1$SNP)
deletion1$DEL<-deletion1$DEL*201/sum(deletion1$DEL)

colnames(deletion2)<-c("label","counts","missing","SNP","DEL")
deletion2 <- cbind(deletion2, type="CM334")
deletion2$SNP<-deletion2$SNP*201/sum(deletion2$SNP)
deletion2$DEL<-deletion2$DEL*201/sum(deletion2$DEL)

colnames(deletion3)<-c("label","counts","missing","SNP","DEL")
deletion3 <- cbind(deletion3, type="Glabr")
deletion3$SNP<-deletion3$SNP*201/sum(deletion3$SNP)
deletion3$DEL<-deletion3$DEL*201/sum(deletion3$DEL)

colnames(deletion4)<-c("label","counts","missing","SNP","DEL")
deletion4 <- cbind(deletion4, type="Chine")
deletion4$SNP<-deletion4$SNP*201/sum(deletion4$SNP)
deletion4$DEL<-deletion4$DEL*201/sum(deletion4$DEL)

colnames(deletion5)<-c("label","counts","missing","SNP","DEL")
deletion5 <- cbind(deletion5, type="Bacca")
deletion5$SNP<-deletion5$SNP*201/sum(deletion5$SNP)
deletion5$DEL<-deletion5$DEL*201/sum(deletion5$DEL)

mydf1 <- rbind(deletion1,deletion2,deletion3,deletion4,deletion5)
ggplot(data=mydf1,aes(x=label, y = DEL, color=type))  + 
  geom_smooth(method="loess",span=0.05,se=F) + 
  scale_color_manual(values=c("#E69F00","#56B4E9","#009E73","#CC79A7","#0072B2")) + 
  labs(color="Types",x="Distance to TAD boundary",y="Obs/Exp deletion coverage") +
  theme( legend.position=c(0.8,0.1),legend.background = element_blank(), legend.title= element_text(size=11, face="bold", color='blue'))

ggplot(data=mydf1,aes(x=label, y = SNP, color=type))  + 
  geom_smooth(method="loess",span=0.05,se=F) + 
  scale_color_manual(values=c("#E69F00","#56B4E9","#009E73","#CC79A7","#0072B2")) + 
  labs(color="Types",x="Distance to TAD boundary",y="Obs/Exp SNP density") +
  theme( legend.position=c(0.8,0.1),legend.background = element_blank(), legend.title= element_text(size=11, face="bold", color='blue'))

#############################################################################################################################################

library(ggplot2)
setwd("/home/yiliao/Documents/Pepper_2021/Figure6/Peak/tomato/Pepper_SNP_DEL/TopDom")

### For close coverage 
deletion1 <-read.table("Zunla.bed.bed.bed")
deletion2 <-read.table("CM334.bed.bed.bed")
deletion3 <-read.table("Glabr.bed.bed.bed")
deletion4 <-read.table("Chine.bed.bed.bed")
deletion5 <-read.table("Bacca.bed.bed.bed")

colnames(deletion1)<-c("label","counts","missing","SNP","DEL")
deletion1 <- cbind(deletion1, type="Zunla")
deletion1$SNP<-deletion1$SNP*201/sum(deletion1$SNP)
deletion1$DEL<-deletion1$DEL*201/sum(deletion1$DEL)

colnames(deletion2)<-c("label","counts","missing","SNP","DEL")
deletion2 <- cbind(deletion2, type="CM334")
deletion2$SNP<-deletion2$SNP*201/sum(deletion2$SNP)
deletion2$DEL<-deletion2$DEL*201/sum(deletion2$DEL)

colnames(deletion3)<-c("label","counts","missing","SNP","DEL")
deletion3 <- cbind(deletion3, type="Glabr")
deletion3$SNP<-deletion3$SNP*201/sum(deletion3$SNP)
deletion3$DEL<-deletion3$DEL*201/sum(deletion3$DEL)

colnames(deletion4)<-c("label","counts","missing","SNP","DEL")
deletion4 <- cbind(deletion4, type="Chine")
deletion4$SNP<-deletion4$SNP*201/sum(deletion4$SNP)
deletion4$DEL<-deletion4$DEL*201/sum(deletion4$DEL)

colnames(deletion5)<-c("label","counts","missing","SNP","DEL")
deletion5 <- cbind(deletion5, type="Bacca")
deletion5$SNP<-deletion5$SNP*201/sum(deletion5$SNP)
deletion5$DEL<-deletion5$DEL*201/sum(deletion5$DEL)

mydf1 <- rbind(deletion1,deletion2,deletion3,deletion4,deletion5)
ggplot(data=mydf1,aes(x=label, y = DEL, color=type))  + 
  geom_smooth(method="loess",span=0.1,se=F) + 
  scale_color_manual(values=c("#E69F00","#56B4E9","#009E73","#CC79A7","#0072B2")) + 
  labs(color="Types",x="Distance to TAD boundary",y="Obs/Exp deletion coverage") +
  theme( legend.position=c(0.8,0.1),legend.background = element_blank(), legend.title= element_text(size=11, face="bold", color='blue'))

ggplot(data=mydf1,aes(x=label, y = SNP, color=type))  + 
  geom_smooth(method="loess",span=0.1,se=F) + 
  scale_color_manual(values=c("#E69F00","#56B4E9","#009E73","#CC79A7","#0072B2")) + 
  labs(color="Types",x="Distance to TAD boundary",y="Obs/Exp SNP density") +
  theme( legend.position=c(0.8,0.1),legend.background = element_blank(), legend.title= element_text(size=11, face="bold", color='blue'))


#################################################################################################
library(ggplot2)
setwd("/home/yiliao/Documents/Pepper_2021/Figure6/Peak/tomato/Pepper_SNP_DEL/Arrowhead")

### For close coverage 
deletion1 <-read.table("Zunla.bed.bed.bed.bed.bed")
deletion2 <-read.table("CM334.bed.bed.bed.bed.bed")
deletion3 <-read.table("Glabr.bed.bed.bed.bed.bed")
deletion4 <-read.table("Chine.bed.bed.bed.bed.bed")
deletion5 <-read.table("Bacca.bed.bed.bed.bed.bed")

colnames(deletion1)<-c("label","counts","missing","SNP","DEL")
deletion1 <- cbind(deletion1, type="Zunla")
deletion1$SNP<-deletion1$SNP*201/sum(deletion1$SNP)
deletion1$DEL<-deletion1$DEL*201/sum(deletion1$DEL)

colnames(deletion2)<-c("label","counts","missing","SNP","DEL")
deletion2 <- cbind(deletion2, type="CM334")
deletion2$SNP<-deletion2$SNP*201/sum(deletion2$SNP)
deletion2$DEL<-deletion2$DEL*201/sum(deletion2$DEL)

colnames(deletion3)<-c("label","counts","missing","SNP","DEL")
deletion3 <- cbind(deletion3, type="Glabr")
deletion3$SNP<-deletion3$SNP*201/sum(deletion3$SNP)
deletion3$DEL<-deletion3$DEL*201/sum(deletion3$DEL)

colnames(deletion4)<-c("label","counts","missing","SNP","DEL")
deletion4 <- cbind(deletion4, type="Chine")
deletion4$SNP<-deletion4$SNP*201/sum(deletion4$SNP)
deletion4$DEL<-deletion4$DEL*201/sum(deletion4$DEL)

colnames(deletion5)<-c("label","counts","missing","SNP","DEL")
deletion5 <- cbind(deletion5, type="Bacca")
deletion5$SNP<-deletion5$SNP*201/sum(deletion5$SNP)
deletion5$DEL<-deletion5$DEL*201/sum(deletion5$DEL)

mydf1 <- rbind(deletion1,deletion2,deletion3,deletion4,deletion5)
ggplot(data=mydf1,aes(x=label, y = DEL, color=type))  + 
  geom_smooth(method="loess",span=0.05,se=F) + 
  scale_color_manual(values=c("#E69F00","#56B4E9","#009E73","#CC79A7","#0072B2")) + 
  labs(color="Types",x="Distance to TAD boundary",y="Obs/Exp deletion coverage") +
  theme( legend.position=c(0.8,0.1),legend.background = element_blank(), legend.title= element_text(size=11, face="bold", color='blue'))

ggplot(data=mydf1,aes(x=label, y = SNP, color=type))  + 
  geom_smooth(method="loess",span=0.2,se=F) + 
  scale_color_manual(values=c("#E69F00","#56B4E9","#009E73","#CC79A7","#0072B2")) + 
  labs(color="Types",x="Distance to TAD boundary",y="Obs/Exp SNP density") +
  theme( legend.position=c(0.8,0.1),legend.background = element_blank(), legend.title= element_text(size=11, face="bold", color='blue'))

library(ggplot2)
setwd("~/Documents/Pepper_2021/Figure6/Peak/tomato")


### For tomato SNPs and InDels
snp<-read.table("all.snp.bed.peak.combined.bed")
colnames(snp)<-c("label","counts")
snp <- cbind(snp, type="SNP")
snp$counts<-snp$counts*201/sum(snp$counts)

del<-read.table("All.raw.deletions.1_49bp.bed.peak.combined.bed")
colnames(del)<-c("label","counts")
del <- cbind(del, type="1-49bp deletions")
del$counts<-del$counts*201/sum(del$counts)

ins<-read.table("All.raw.ins.bed.peak.combined.bed")
colnames(ins)<-c("label","counts")
ins <- cbind(ins, type="1-49bp insertions")
ins$counts<-ins$counts*201/sum(ins$counts)

mydata <- rbind(snp,del,ins)
ggplot(data=mydata,aes(x=label, y=counts, color=type))  + 
  geom_smooth(method="loess",span=0.1,se=F) + 
  scale_color_manual(values=c("#E69F00","#0072B2","#CC79A7")) + 
  labs(color="Types",x="Distance to TAD boundary",y="Obs/Exp SNP and InDel density") +
  theme( legend.position=c(0.8,0.1),legend.background = element_blank(), legend.title= element_text(size=11, face="bold", color='blue'))

### For SVs
SVDel<-read.table("SVDel.breaks.bed")
colnames(SVDel)<-c("label","counts")
SVDel <- cbind(SVDel, type="50-20kbp deletion density")
SVDel$counts<-SVDel$counts*201/sum(SVDel$counts)

SVDelCov<-read.table("SVDel.covs.bed")
colnames(SVDelCov)<-c("label","counts")
SVDelCov <- cbind(SVDelCov, type="50-20kbp deletion coverage")
SVDelCov$counts<-SVDelCov$counts*201/sum(SVDelCov$counts)

SVIns<-read.table("SVIns.breaks.bed")
colnames(SVIns)<-c("label","counts")
SVIns <- cbind(SVIns, type="50-20kbp insertion density")
SVIns$counts<-SVIns$counts*201/sum(SVIns$counts)

SVInsCov<-read.table("SVIns.covs.bed")
colnames(SVInsCov)<-c("label","counts")
SVInsCov <- cbind(SVInsCov, type="50-20kbp insertion coverage")
SVInsCov$counts<-SVInsCov$counts*201/sum(SVInsCov$counts)


mydata <- rbind(SVDel,SVDelCov,SVIns,SVInsCov)
ggplot(data=mydata,aes(x=label, y=counts, color=type))  + 
  geom_smooth(method="loess",span=0.1,se=F) + 
  scale_color_manual(values=c("#E69F00","#0072B2","#CC79A7","#56B4E9")) + 
  labs(color="Types",x="Distance to TAD boundary",y="Obs/Exp SNP and InDel density") +
  theme( legend.position=c(0.8,0.1),legend.background = element_blank(), legend.title= element_text(size=11, face="bold", color='blue'))

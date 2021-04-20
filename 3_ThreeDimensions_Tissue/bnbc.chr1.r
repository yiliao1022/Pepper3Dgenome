library(bnbc)

setwd("/home/yiliao/OLDISK/genome_assembly/hic_explorer/7_hic2cool/40k/cooltools_raw/ajust/Final/Chr01")

df1<-read.table("/home/yiliao/OLDISK/genome_assembly/hic_explorer/7_hic2cool/40k/cooltools_raw/ajust/Final/LJGR1/LJGR1.Ca_59Chr01.tsv.sym.mat")
df2<-read.table("/home/yiliao/OLDISK/genome_assembly/hic_explorer/7_hic2cool/40k/cooltools_raw/ajust/Final/LJGR2/LJGR2.Ca_59Chr01.tsv.sym.mat")
df3<-read.table("/home/yiliao/OLDISK/genome_assembly/hic_explorer/7_hic2cool/40k/cooltools_raw/ajust/Final/LJHL1/LJHL1.Ca_59Chr01.tsv.sym.mat")
df4<-read.table("/home/yiliao/OLDISK/genome_assembly/hic_explorer/7_hic2cool/40k/cooltools_raw/ajust/Final/LJHL2/LJHL2.Ca_59Chr01.tsv.sym.mat")
df5<-read.table("/home/yiliao/OLDISK/genome_assembly/hic_explorer/7_hic2cool/40k/cooltools_raw/ajust/Final/LJTZ1/LJTZ1.Ca_59Chr01.tsv.sym.mat")
df6<-read.table("/home/yiliao/OLDISK/genome_assembly/hic_explorer/7_hic2cool/40k/cooltools_raw/ajust/Final/LJTZ2/LJTZ2.Ca_59Chr01.tsv.sym.mat")
df7<-read.table("/home/yiliao/OLDISK/genome_assembly/hic_explorer/7_hic2cool/40k/cooltools_raw/ajust/Final/LYJP1/LJYP1.Ca_59Chr01.tsv.sym.mat")
df8<-read.table("/home/yiliao/OLDISK/genome_assembly/hic_explorer/7_hic2cool/40k/cooltools_raw/ajust/Final/LYJP2/LYJP2.Ca_59Chr01.tsv.sym.mat")

mat1<-as.matrix(df1)
mat2<-as.matrix(df2)
mat3<-as.matrix(df3)
mat4<-as.matrix(df4)
mat5<-as.matrix(df5)
mat6<-as.matrix(df6)
mat7<-as.matrix(df7)
mat8<-as.matrix(df8)

mat1<-mat1[,-8331]
mat2<-mat2[,-8331]
mat3<-mat3[,-8331]
mat4<-mat4[,-8331]
mat5<-mat5[,-8331]
mat6<-mat6[,-8331]
mat7<-mat7[,-8331]
mat8<-mat8[,-8331]

Matlst<-list(mat1,mat2,mat3,mat4,mat5,mat6,mat7,mat8)

range <-c()
for (i in 1:8330) {range[i]=((i-1)*40000)+1}
rangelst <- GRanges(seqnames = "Ca_59Chr01",ranges = IRanges(start = range,width = 40000))

df<-DataFrame(name=c("LJGR1","LJGR2","LJHL1","LJHL2","LJTZ1","LJTZ2","LJYP1","LJYP2"),batch=c("A","B","A","A","A","B","A","B"),row.names = c("GR1","GR2","HL1","HL2","TZ1","TZ2","YP1","YP2"))


Pepper <-ContactGroup(rowData = rangelst, contacts=Matlst, colData = df)
Pepper.cpm<-logCPM(Pepper)
batches<-colData(Pepper)$batch
Pepper.smooth<-boxSmoother(Pepper.cpm,5,mc.cores=1)
Pepper.bnbc<-bnbc(Pepper.smooth,batches,threshold=2e7,step=4e4,bstart=2,nbands=8329)

write.csv(contacts(Pepper.bnbc)[[1]],"matrixGR1.chr01.csv")
write.csv(contacts(Pepper.bnbc)[[2]],"matrixGR2.chr01.csv")
write.csv(contacts(Pepper.bnbc)[[3]],"matrixHL1.chr01.csv")
write.csv(contacts(Pepper.bnbc)[[4]],"matrixHL2.chr01.csv")
write.csv(contacts(Pepper.bnbc)[[5]],"matrixTZ1.chr01.csv")
write.csv(contacts(Pepper.bnbc)[[6]],"matrixTZ2.chr01.csv")
write.csv(contacts(Pepper.bnbc)[[7]],"matrixYP1.chr01.csv")
write.csv(contacts(Pepper.bnbc)[[8]],"matrixYP2.chr01.csv")

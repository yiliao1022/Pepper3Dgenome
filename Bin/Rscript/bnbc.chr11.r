library(bnbc)

setwd("/home/yiliao/OLDISK/genome_assembly/hic_explorer/13_bnbc/100k/Chr11")

df1<-read.table("/home/yiliao/OLDISK/genome_assembly/hic_explorer/13_bnbc/100k/LJGR1/LJGR1.Ca_59Chr11.tsv.sym.mat")
df2<-read.table("/home/yiliao/OLDISK/genome_assembly/hic_explorer/13_bnbc/100k/LJGR2/LJGR2.Ca_59Chr11.tsv.sym.mat")
df3<-read.table("/home/yiliao/OLDISK/genome_assembly/hic_explorer/13_bnbc/100k/LJHL1/LJHL1.Ca_59Chr11.tsv.sym.mat")
df4<-read.table("/home/yiliao/OLDISK/genome_assembly/hic_explorer/13_bnbc/100k/LJHL2/LJHL2.Ca_59Chr11.tsv.sym.mat")
df5<-read.table("/home/yiliao/OLDISK/genome_assembly/hic_explorer/13_bnbc/100k/LJTZ1/LJTZ1.Ca_59Chr11.tsv.sym.mat")
df6<-read.table("/home/yiliao/OLDISK/genome_assembly/hic_explorer/13_bnbc/100k/LJTZ2/LJTZ2.Ca_59Chr11.tsv.sym.mat")
df7<-read.table("/home/yiliao/OLDISK/genome_assembly/hic_explorer/13_bnbc/100k/LJYP1/LJYP1.Ca_59Chr11.tsv.sym.mat")
df8<-read.table("/home/yiliao/OLDISK/genome_assembly/hic_explorer/13_bnbc/100k/LJYP2/LJYP2.Ca_59Chr11.tsv.sym.mat")

mat1<-as.matrix(df1)
mat2<-as.matrix(df2)
mat3<-as.matrix(df3)
mat4<-as.matrix(df4)
mat5<-as.matrix(df5)
mat6<-as.matrix(df6)
mat7<-as.matrix(df7)
mat8<-as.matrix(df8)

mat1<-mat1[,-2761]
mat2<-mat2[,-2761]
mat3<-mat3[,-2761]
mat4<-mat4[,-2761]
mat5<-mat5[,-2761]
mat6<-mat6[,-2761]
mat7<-mat7[,-2761]
mat8<-mat8[,-2761]

Matlst<-list(mat1,mat2,mat3,mat4,mat5,mat6,mat7,mat8)

range <-c()
for (i in 1:2760) {range[i]=((i-1)*100000)+1}
rangelst <- GRanges(seqnames = "Ca_59Chr11",ranges = IRanges(start = range,width = 100000))

df<-DataFrame(name=c("LJGR1","LJGR2","LJHL1","LJHL2","LJTZ1","LJTZ2","LJYP1","LJYP2"),batch=c("A","B","A","A","A","B","A","B"),row.names = c("GR1","GR2","HL1","HL2","TZ1","TZ2","YP1","YP2"))


Pepper <-ContactGroup(rowData = rangelst, contacts=Matlst, colData = df)
Pepper.cpm<-logCPM(Pepper)
batches<-colData(Pepper)$batch
Pepper.smooth<-boxSmoother(Pepper.cpm,5,mc.cores=1)
Pepper.bnbc<-bnbc(Pepper.smooth,batches,threshold=2e7,step=1e5,bstart=2,nbands=2759)

write.csv(contacts(Pepper.bnbc)[[1]],"matrixGR1.Chr11.csv")
write.csv(contacts(Pepper.bnbc)[[2]],"matrixGR2.Chr11.csv")
write.csv(contacts(Pepper.bnbc)[[3]],"matrixHL1.Chr11.csv")
write.csv(contacts(Pepper.bnbc)[[4]],"matrixHL2.Chr11.csv")
write.csv(contacts(Pepper.bnbc)[[5]],"matrixTZ1.Chr11.csv")
write.csv(contacts(Pepper.bnbc)[[6]],"matrixTZ2.Chr11.csv")
write.csv(contacts(Pepper.bnbc)[[7]],"matrixYP1.Chr11.csv")
write.csv(contacts(Pepper.bnbc)[[8]],"matrixYP2.Chr11.csv")
 

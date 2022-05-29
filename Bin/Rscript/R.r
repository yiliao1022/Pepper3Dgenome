setwd("/home/yiliao/Documents/Pepper_2021/ABcompartments/draw/FINAL/All")

library("edgeR")

rawRNA <- read.table("RNA.40k.bed.bed.bed.bed", header=T, row.names=1, com='', quote='',
                   check.names=F, sep="\t")
rawAB <- read.table()
cpmdata <- edgeR::cpm(rawdata)
cpmdata_rd <- round(cpmdata)
write.table(cpmdata_rd,"RNA.40K.CPM.bed",sep="\t")

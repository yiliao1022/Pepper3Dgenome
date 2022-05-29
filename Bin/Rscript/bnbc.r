#BiocManager::install("bnbc", configure.args="--disable-threading")
library(bnbc)
data(cgEx)
cgEx.cpm<-logCPM(cgEx)
batches<-colData(cgEx)$Batch
cgEx.smooth<-boxSmoother(cgEx.cpm,5,mc.cores=1)
cgEx.bnbc<-bnbc(cgEx.smooth,batches,threshold=1e7,step=4e4,bstart=2,nbands=4)

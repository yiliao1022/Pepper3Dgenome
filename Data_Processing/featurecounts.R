#!/usr/bin/env Rscript
#parse parameter
library(argparser,quietly = TRUE)
#Create a parser
p<- arg_parser("run featureCounts and calculate FPKM/TPM")
p<-add_argument(p,"--bam",help="input: bam file",type="character")
p<-add_argument(p,"--gtf",help="input: gtf file",type="character")
p<-add_argument(p,"--output",help="out prefix",type="character")
p<-add_argument(p,"--threads",help="threads",type="character")
#Parse the command line arguments
argv<-parse_args(p)
library(Rsubread)
library(limma)
library(edgeR)
bamFile<- argv$bam
gtfFile<- argv$gtf
nthreads<- argv$threads
outFilePref<- argv$output
outStatsFilePath<- paste(outFilePref, '.log',sep='');
outCountsFilePath<- paste(outFilePref,'.count',sep='');
fCountsList=featureCounts(bamFile,annot.ext=gtfFile,isGTFAnnotationFile=TRUE,nthreads=nthreads,isPairedEnd=TRUE,tmpDir='/mnt/memorydisk')
dgeList=DGEList(counts=fCountsList$counts,genes=fCountsList$annotation)
fpkm=rpkm(dgeList,dgeList$genes$Length)
tpm=exp(log(fpkm)-log(sum(fpkm))+log(1e6))
write.table(fCountsList$stat,outStatsFilePath,sep="\t",col.names = FALSE,row.names = FALSE,quote=FALSE)
featureCounts=cbind(fCountsList$annotation[,1],fCountsList$counts,fpkm,tpm)
colnames(featureCounts)=c('gene_id','counts','fpkm','tpm')
write.table(featureCounts,outCountsFilePath,sep="\t",col.names = TRUE,row.names = FALSE,quote = FALSE)

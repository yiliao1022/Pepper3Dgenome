
#Step1: contact matrices construction using HiCExplorer




#Step2: converting HiCExplorer output matrices .h5 (40kb resolution) to symmetric matrices which are required for bnbc-package (R) in performing normalization and batch correction across samples
### (hicexplorer3.5.3)
#1. only keep individual chromosome each time if you just need to look for compartments or subcomparts or TADs;

for i in Ca_59Chr01 Ca_59Chr02 Ca_59Chr03 Ca_59Chr04 Ca_59Chr05 Ca_59Chr06 Ca_59Chr07 Ca_59Chr08 Ca_59Chr09 Ca_59Chr10 Ca_59Chr11 Ca_59Chr12; do hicAdjustMatrix -m ../../../../../LJYPCF_mapping/matrix/hic_matrix_40000_3.53.h5 --chromosomes $i --action keep --outFileName LYJP2.$i.h5;done
for i in Ca_59Chr01 Ca_59Chr02 Ca_59Chr03 Ca_59Chr04 Ca_59Chr05 Ca_59Chr06 Ca_59Chr07 Ca_59Chr08 Ca_59Chr09 Ca_59Chr10 Ca_59Chr11 Ca_59Chr12; do hicConvertFormat -m LYJP2.$i.h5 --inputFormat h5 --outputFormat ginteractions --outFileName LYJP2.$i;done

#2. make the matrix square
perl GetMat.pl LYJP2.Ca_59Chr01.tsv 8330 40000  ## 8330 are the bin number of chromosome 1 which obtained by dividing chromosome size with 40000;
perl GetMat.pl LYJP2.Ca_59Chr02.tsv 4463 40000
perl GetMat.pl LYJP2.Ca_59Chr03.tsv 7234 40000
perl GetMat.pl LYJP2.Ca_59Chr04.tsv 6266 40000
perl GetMat.pl LYJP2.Ca_59Chr05.tsv 6361 40000
perl GetMat.pl LYJP2.Ca_59Chr06.tsv 6278 40000
perl GetMat.pl LYJP2.Ca_59Chr07.tsv 6689 40000
perl GetMat.pl LYJP2.Ca_59Chr08.tsv 4351 40000
perl GetMat.pl LYJP2.Ca_59Chr09.tsv 7006 40000
perl GetMat.pl LYJP2.Ca_59Chr10.tsv 6212 40000
perl GetMat.pl LYJP2.Ca_59Chr11.tsv 6900 40000
perl GetMat.pl LYJP2.Ca_59Chr12.tsv 6551 40000

#3 Using bnbc for normalizing and batch correction

Rscript bnbc.chr1.r


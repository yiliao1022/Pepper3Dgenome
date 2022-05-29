####Run G.C.E (Shell)
trimmomatic PE survey_R1.fq.gz survey_R2.fq.gz  -baseout survey.clean.fq.gz ILLUMINACLIP:adapters.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36
ls *P.fq.gz > reads_list_file
kmer_number=323426998344
prefix=Pepper_survey
cov= 100
kmerfreq -k 17 -t 48 -r 10000 -p Pepper_survey reads_list_file
less ${prefix}.kmer.freq.stat | perl -ne 'next if(/^#/ || /^\s/); print; ' | awk '{print $1"\\t''$2}' > ${prefix}.kmer.freq.stat.2colum
gce -f ${prefix}.kmer.freq.stat.2colum -g ${kmer_number} -m 1 -D 8 -b 0 -H 1 -c ${cov} 1> ${prefix}.table 2 > ${prefix}.gce.result


####Run CANU (Shell)
canu -trim-assemble -p Capsicum -d Capsicum GenomeSize=3000m corMhapFilterThreshold=0.0000000002 corMhapOptions="""--threshold 0.80 --num-hashes 512 --num-min-matches 3 --ordered-sketch-size 1000 --ordered-kmer-size 17 --min-olap-length 2000 --repeat-idf-scale 50""" mhapBlockSize=500 ovlMerThreshold=500 minReadLength=30000 minOverlapLength=2000   -pacbio-corrected cns_final.fasta &>>canu.log


####Run Pilon and BUSCO (Shell)
bowtie2 --threads 40 -x contigs.bowtie2 -1 Survey_R1.fastq.gz -2 Survey_R1.fastq.gz 2> bowtie2.log |samtools view -hF 256 -  |  samtools sort -@ 10 -m 20G -o contigs.bowtie2mapping.bam
samtools index -@ 10 contigs.bowtie2mapping.bam
java -Xmx 900G -jar pilon.jar --genome contigs.fa --frags contigs.bowtie2mapping.bam --fix snps,indels --output pilon_polished.fa &> pilon.log
run_BUSCO.py -i pilon_polished.fa -l ~/database/software_database/Busco_database/embryophyta_odb9 -o pilon_polished.checkresult -m genome -c ${threads} -f
bowtie2-build Final_polish.dna.fa Final_polish.bowtie2
bowtie2 --threads 40 -x Final_polish.bowtie2 -1 Survey_R1.fastq.gz -2 Survey_R1.fastq.gz 2> bowtie2.log |samtools view -hF 256 -  |  samtools sort -@ 10 -m 20G -o QV_mapping.sorted.bam 
bcftools mpileup -Ou QV_mapping.sorted.bam -f Ca_59.dna.fa | bcftools call -mv -o QV.vcf



###Run Juicer and 3D-DNA (Shell)
python3 juicer/misc/generate_site_positions.py MboI hic polished_contigs.fa
samtools faidx polished_contigs.fa ; cut -f1,2 polished_contigs.fa.fai > polished_contigs.fa.size
bwa index $contigs
bash juicer.sh -g contig_ -d `pwd` -s MboI -z polished_contigs.fa -t 40 -y hic_MboI.txt -p polished_contigs.fa.size
mkdir ../3ddna && cd ../3ddna
bash 3d-dna/run-asm-pipeline.sh -r 0 ../polished_contigs.fa ../aligned/merged_nodups.txt
###The 3D-DNA output files, polished_contigs.final.assembly and polished_contigs.final.hic, were loaded into the Juicebox Assembly Tools (JBAT) for manually verifying and correcting the above 3D-DNA scaffolds to produce the final corrected assembly. 
bash 3d-dna/run-asm-pipeline-post-review.sh -r pepper_new_map.review.assembly ../polished_contigs.fa ../aligned/merged_nodups.txt


###Run EDTA (Shell)
EDTA.pl -specie others -threads 48 -overwrite 1 -genome ${genome_file}





###Build Kimura distance (Shell)
/PATH/TO/EDTA/EDTA.pl -specie others -threads $threads -overwrite 1 -genome Ca_59.dna.fa
grep LTR_retrotransposon$'\t' Ca_59.dna.fa.mod.EDTA.intact.gff3 | cut -f 1,4,5 | sed 's;\t;:;' | sed 's;\t;-;' > LTR.region
samtools faidx Ca_59.dna.fa -r LTR.region> LTR.fa
paste <(grep LTR_retrotransposon$'\t' Ca_59.dna.fa.mod.EDTA.intact.gff3 | cut -f 1,4,5 | sed 's;\t;:;' | sed 's;\t;-;' ) <(grep LTR_retrotransposon$'\t' Ca_59.dna.fa.mod.EDTA.intact.gff3 | cut -f 9 | cut -d\; -f 4 | cut -d= -f 2) > rename.table
IFS=$'\n';
for i in `cat rename.table`;do
raw_name=`echo $i | cut -f 1`;
type=`echo $i | cut -f 2`;
sed -i "s;$raw_name;${raw_name}#${type};" LTR.fa;
done
/PATH/TO/REPEATMASKER/RepeatMasker  -norna -nolow -pa 24 -s -xsmall -gff -lib LTR.fa Ca_59.dna.fa
/PATH/TO/REPEATMASKER/util/calcDivergenceFromAlign.pl -s ltr.div Ca_59.dna.fa.cat.gz
genome_size=3000000000
output_pdf=output.pdf
tail -n 72 ltr.div > genome.Kimura.distance
/usr/bin/env R --no-save <<EOF
library(reshape)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(tidyverse)
library(gridExtra)
sessionInfo()
#R version 3.6.1 (2019-07-05)
#Platform: x86_64-conda_cos6-linux-gnu (64-bit)
#Running under: Ubuntu 20.04.1 LTS
#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     
#other attached packages:
#  [1] gridExtra_2.3     forcats_0.5.0     stringr_1.4.0     dplyr_0.8.5       purrr_0.3.4      
#  [6] readr_1.3.1       tidyr_1.1.0       tibble_3.0.1      tidyverse_1.3.0   hrbrthemes_0.8.0 
#  [11] viridis_0.5.1     viridisLite_0.3.0 ggplot2_3.3.0     reshape_0.8.8    
KimuraDistance <- read.csv("genome.Kimura.distance",sep=" ")
pdf("$output_pdf")
#add here the genome size in bp
genomes_size=$genome_size
kd_melt = melt(KimuraDistance,id="Div")
kd_melt\$norm = kd_melt\$value/genomes_size * 100
ggplot(kd_melt, aes(fill=variable, y=norm, x=Div)) + 
  geom_bar(position="stack", stat="identity",color="black") +
  scale_fill_viridis(discrete = T) +
  theme_classic() +
  xlab("Kimura substitution level") +
  ylab("Percent of the genome") + 
  labs(fill = "") +
  coord_cartesian(xlim = c(0, 55)) +
  theme(axis.text=element_text(size=11),axis.title =element_text(size=12))
dev.off()
EOF


###Iso-seq transcriptome data process.(Shell)
for prefix in pulp placenta roots leaves bud; do
ccs ${prefix}.subreads.bam ${prefix}.ccs.bam --noPolish --minPasses 1
lima ${prefix}.ccs.bam primers.fa ${prefix}.demux.ccs.bam --isoseq --peek-guess
isoseq3 refine --require-polya ${prefix}.demux.ccs.F1_5p--R1_3p.bam primers.fa ${prefix}.flnc.bam
samtools fasta ${prefix}.flnc.bam > ${prefix}.flnc.fa
minimap2 -x splice -a Ca_59.dna.fa ${prefix}.flnc.fa | samtools sort -@ 15 -m 10G -o ${prefix}.flnc.mapping.sorted.bam
done
for prefix in pulp placenta roots leaves bud; do
python tama/tama_collapse.py -b BAM -s ../${prefix}.flnc.mapping.sorted.out.bam -f ${genomefa} -x no_cap -e common_ends -p ${prefix}_nocap_isoform > ${prefix}.nocap_tama.collapse.log
done
find ./ -mindepth 1 -iname "*_nocap_isoform.bed" > tama_nocap.list.txt
python tama/tama_merge.py -f tama_nocap.list.txt -p tama_nocap
bedToGenePred tama_nocap_isoform.bed tama_nocap_isoform.genepred
genePredToGtf 'file' tama_nocap_isoform.genepred tama_nocap_isoform.gtf
gffread -g Ca_59.dna.fa -w Ca_59.isoseq.mRNA.fa tama_nocap_isoform.gtf
transcriptomefa=Ca_59.isoseq.mRNA.fa
TransDecoder.LongOrfs -t ${transcriptomefa} -m 100 -G universal -O ${transcriptomefa}.transdecoder
hmmscan --domtblout ${transcriptomefa}.hmmscan pfam.hmm ${transcriptomefa}.transdecoder/longest.pep
blastp -query ${transcriptomefa}.transdecoder/longest.pep -db swissprot.blastdb -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 12 > ${transcriptomefa}.blastres
TransDecoder.Predict -t ${transcriptomefa}.transdecoder --retain_blastp_hits ${transcriptomefa}.blastres --retain_pfam_hits ${transcriptomefa}.hmmscan -O ${transcriptomefa}.transdecoder --single_best_only






### Mapping RNA-reads to the CA59 genome (Shell)
hisat2-build -p ${threads} ${genome_file} tmp/tmpidx
hisat2 --dta --rg-id hisat2 --rg SM:${samplename} --threads ${threads}  -x tmp/tmpidx -1 ${reads_R1} -2 ${reads_R2} | samtools view -Shb - > hisat2/${samplename}.unsort.bam
samtools sort -@ ${threads} hisat2/${samplename}.unsort.bam > hisat2/${samplename}.sorted.bam
### Reconstruction of a new transcriptome using StringTie (Shell) 
stringtie hisat2/${samplename}.sorted.bam  -p ${threads} -o stringtie/${samplename}.gtf -A stringtie/${samplename}.tab
find stringtie/ -iname  "*.gtf" > assemblies.txt
stringtie --merge -l Unigene -p ${threads} -G ${gtffile} -o stringtie_merged.gtf assemblies.txt
## Quantification of gene/transcripts expression level (Shell)
Rsciprt featurecounts.R --bam ${samplename}.sorted.bam --gtf $gtffile --threads 20 --output ${samplename}.countmatrix

####featurecounts.R was modified from https://www.jianshu.com/p/46fa4e699187 and provided as below(Rscript):
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

###The R package ggbiplot and following R script were used to build PCA plot:
###PCA result of RNA seq(CITATION ggbiplot : https://github.com/vqv/ggbiplot)
library("ggbiplot")
pca <- read.table("featureCounts_quantify.result",header=T,row.names=1)
pca <-  pca[rowSums(pca != 0) != 0, ]
fpca <- t(as.matrix(pca))
fp.pca <- prcomp(fpca, scale = TRUE)
fp.class<-factor(c("Pulp ","Pulp ","Pulp ","Flower bud ","Flower bud ","Flower bud ","Roots","Roots","Roots","Placenta ","Placenta ","Placenta ","Leaves","Leaves","Leaves")) 
PCA3 <- ggbiplot(fp.pca, obs.scale = 1, var.scale = 1,groups = fp.class, ellipse = TRUE, circle = TRUE,var.axes=F) +scale_color_discrete(name = '') + theme(legend.direction = 'horizontal', legend.position = 'top')
ggsave("pca.png",PCA3)




##########MAKER command(Shell):
mpiexec -n ${threads} maker -fix_nucleotides 1> maker.log 2> maker.err


#############maker_bopts.ctl config file in round 1 shows as following:
#-----Genome (these are always required)
genome=Ca_59.dna.fa #genome sequence (fasta file or fasta embedded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic
 
#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anything else in maker_gff: 1 = yes, 0 = no
 
#-----EST Evidence (for best results provide a file for at least one)
est=isoseq_clusted.fa #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closely related species in GFF3 format
 
#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=Ca_zunla.pep.fa  #protein sequence file in fasta format (i.e. from multiple organisms)
protein_gff=  #aligned protein homology evidence from an external GFF3 file
 
#-----Repeat Masking (leave values blank to skip repeat masking)
model_org=all #select a model organism for RepBase masking in RepeatMasker
rmlib=Ca_59.dna.fa.mod.EDTA.TElib.fa #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein=~/database/miniconda3/envs/maker/te_proteins.fasta #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff= #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmasker prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)
 
#-----Gene Prediction
snaphmm= #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species= #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no
 
#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file
 
#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of CPUs to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)
 
#-----MAKER Behavior Options
max_dna_len=1000000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=10000 #skip genome contigs below this length (under 10kb are often useless)
 
pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)
 
split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes
 
tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP=/mnt/memorydisk #specify a directory other than the system default temporary directory for temporary files

#############maker_bopts.ctl config file in round 2 and 3 shows as following:
#-----Genome (these are always required)
genome=Ca_59.dna.fa #genome sequence (fasta file or fasta embedded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic
 
#-----Re-annotation Using MAKER Derived GFF3
maker_gff=Ca_59.dna.all.gff.round1 #MAKER derived GFF3 file, which was generated from the last round and will be changed to Ca_59.dna.all.gff.round2 in round 3.
est_pass=1 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=1 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=1 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=1 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=1 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=1 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=1 #passthrough anything else in maker_gff: 1 = yes, 0 = no
 
#-----EST Evidence (for best results provide a file for at least one)
est=isoseq_clusted.fa #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closely related species in GFF3 format
 
#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=Ca_zunla.pep.fa  #protein sequence file in fasta format (i.e. from multiple organisms)
protein_gff=  #aligned protein homology evidence from an external GFF3 file
 
#-----Repeat Masking (leave values blank to skip repeat masking)
model_org=all #select a model organism for RepBase masking in RepeatMasker
rmlib=Ca_59.dna.fa.mod.EDTA.TElib.fa #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein=~/maker/te_proteins.fasta #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff= #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmasker prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)
 
#-----Gene Prediction
snaphmm= #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species= #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no
 
#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file
 
#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)
 
#-----MAKER Behavior Options
max_dna_len=1000000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=10000 #skip genome contigs below this length (under 10kb are often useless)
 
pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)
 
split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes
 
tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP=/mnt/memorydisk #specify a directory other than the system default temporary directory for temporary files

##########Run HiCexplorer (Shell):
bwa mem -t 24 -A1 -B4 -E50 -L0 $reference_fasta ${prefix}_mapping/${prefix}_R1.fastq.gz 2> ${prefix}_mapping/${prefix}_R1.log | samtools view -Shb - > ${prefix}_mapping/${prefix}_R1.bam
bwa mem -t 24 -A1 -B4 -E50 -L0 $reference_fasta ${prefix}_mapping/${prefix}_R2.fastq.gz 2> ${prefix}_mapping/${prefix}_R2.log | samtools view -Shb - > ${prefix}_mapping/${prefix}_R2.bam
mkdir mapped_files QC bamfile matrix
hicFindRestSite --fasta $reference_fasta --searchPattern GATC -o $reference_fasta.rest_site_positions.bed
hicBuildMatrix --danglingSequence GATC --samFiles mapped_files/HiC_R1.bam mapped_files/HiC_R2.bam --binSize ${binSize} --restrictionSequence GATC --threads 8 --inputBufferSize 100000 -o matrix/hic_matrix_${binSize}_3.53.h5  --restrictionCutFile ../Ca_59.dna.fa_rest_site_positions.bed --QCfolder QC/${binSize}_3.53

###########Run Juicer (Shell):
#The Hi-C fastq files should be placed in the directory called “fastq”. All of the details of the process running could be found on the Github page: https://github.com/aidenlab/juicer.
python juicer/misc/generate_site_positions.py MboI hic $genome_file
samtools faidx $genome_file ; cut -f1,2 ${genome_file}.fai > ${genome_file}.size
bwa index $genome_file
bash juicer/CPU/juicer.sh -g contig_ -d `pwd` -s MboI -z $genome_file -t ${threads} -y hic_MboI.txt -p ${genome_file}.size

#############Run Bismark (Shell):
bismark_genome_preparation --path_to_aligner miniconda3/envs/methylation/bin/ --verbose . --parallel 20
bismark --gzip --parallel 30 --genome . -1 $R1.fq.gz -2 $R2.fq.gz
deduplicate_bismark --bam $R1_bismark_bt2_pe.bam
bismark_methylation_extractor --gzip --bedGraph $R1_bismark_bt2_pe.deduplicated.bam
bismark2report

###########Run MACS2 (Shell):
bwa mem -t 24 -M -R “@RG\\tID:${sample}\\tLB:${sample}\\tSM:${sample}\\tPL:ILLUMINA” ${genome_file} ${file1} ${file2} |samtools sort -@ 20 -m 10G > /mnt/memorydisk/${sample}/${sample}.sort.bam
macs2 callpeak -t $prefix.sort.bam -c ${prefix}input.sort.bam -f BAMPE -g 3e9 -n $prefix.contain_input -q 0.05   --shift -100 --extsize 200 --nomodel -B

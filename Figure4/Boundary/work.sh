rsync -auP -e 'ssh -p 10089' yiliao@403lib.f3322.net:/home/yiliao/OLDISK/genome_assembly/hic_explorer/16_TADclusters/YP2/Final/TAD_TADboundaries/Boundaries/*.cov.bed.bed ./
for i in `ls *.bed.bed`; do paste cor.bed $i >> $i.bed;done

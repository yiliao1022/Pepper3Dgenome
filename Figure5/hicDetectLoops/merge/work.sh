cat ./LJHL/LJHL.10_15_20_25.merged.sort.csv ./LJGR/LJGR.10_15_20_25.merged.sort.csv ./LJYP/LJYP.10_15_20_25.merged.sort.csv > LJHL.LJGR.LJYP.10_15_20_25.csv
cat ./LJHL/LJHL.10_15_20_25.merged.sort.csv ./LJTZ/10_15_20_25.merge.sort.bed ./LJYP/LJYP.10_15_20_25.merged.sort.csv > LJHL.LJTZ.LJYP.10_15_20_25.csv
cat ./LJGR/LJGR.10_15_20_25.merged.sort.csv ./LJTZ/10_15_20_25.merge.sort.bed ./LJYP/LJYP.10_15_20_25.merged.sort.csv > LJGR.LJTZ.LJYP.10_15_20_25.csv
cat ./LJHL/LJHL.10_15_20_25.merged.sort.csv ./LJTZ/10_15_20_25.merge.sort.bed ./LJGR/LJGR.10_15_20_25.merged.sort.csv > LJHL.LJTZ.LJGR.10_15_20_25.csv
~/OLDISK/software/pgltools/sh/pgltools sort LJHL.LJTZ.LJGR.10_15_20_25.csv > LJHL.LJTZ.LJGR.10_15_20_25.sort.csv
~/OLDISK/software/pgltools/sh/pgltools sort LJHL.LJTZ.LJYP.10_15_20_25.csv > LJHL.LJTZ.LJYP.10_15_20_25.sort.csv
~/OLDISK/software/pgltools/sh/pgltools sort LJGR.LJTZ.LJYP.10_15_20_25.csv > LJGR.LJTZ.LJYP.10_15_20_25.sort.csv
~/OLDISK/software/pgltools/sh/pgltools intersect -a ./LJTZ/10_15_20_25.merge.sort.bed -b LJHL.LJGR.LJYP.10_15_20_25.sort.csv -wo -d 25000 | cut -f1-6 | sort -k1,1 -k2,2n | sort | uniq > LJTZ.shared.loop.csv
~/OLDISK/software/pgltools/sh/pgltools intersect -a ./LJTZ/10_15_20_25.merge.sort.bed -b LJHL.LJGR.LJYP.10_15_20_25.sort.csv -v -d 25000 | cut -f1-6 | sort -k1,1 -k2,2n | sort | uniq > LJTZ.specific.loop.csv
~/OLDISK/software/pgltools/sh/pgltools intersect -a ./LJGR/LJGR.10_15_20_25.merged.sort.csv -b LJHL.LJTZ.LJYP.10_15_20_25.sort.csv -wo -d 25000 | cut -f1-6 | sort -k1,1 -k2,2n | uniq > LJGR.shared.loop.csv
~/OLDISK/software/pgltools/sh/pgltools intersect -a ./LJGR/LJGR.10_15_20_25.merged.sort.csv -b LJHL.LJTZ.LJYP.10_15_20_25.sort.csv -v -d 25000 | cut -f1-6 | sort -k1,1 -k2,2n | uniq > LJGR.specific.loop.csv
~/OLDISK/software/pgltools/sh/pgltools intersect -a ./LJHL/LJHL.10_15_20_25.merged.sort.csv -b LJGR.LJTZ.LJYP.10_15_20_25.sort.csv -v -d 25000 | cut -f1-6 | sort -k1,1 -k2,2n | uniq > LJHL.specific.loop.csv
~/OLDISK/software/pgltools/sh/pgltools intersect -a ./LJHL/LJHL.10_15_20_25.merged.sort.csv -b LJGR.LJTZ.LJYP.10_15_20_25.sort.csv -wo -d 25000 | cut -f1-6 | sort -k1,1 -k2,2n | uniq > LJHL.shared.loop.csv
~/OLDISK/software/pgltools/sh/pgltools sort LJYP.shared.loop.csv > LJYP.shared.loop.sort.csv
~/OLDISK/software/pgltools/sh/pgltools sort LJYP.specific.loop.csv > LJYP.specific.loop.sort.csv
~/OLDISK/software/pgltools/sh/pgltools sort LJGR.shared.loop.csv > LJGR.shared.loop.sort.csv
~/OLDISK/software/pgltools/sh/pgltools sort LJGR.specific.loop.csv > LJGR.specific.loop.sort.csv
~/OLDISK/software/pgltools/sh/pgltools sort LJHL.shared.loop.csv > LJHL.shared.loop.sort.csv
~/OLDISK/software/pgltools/sh/pgltools sort LJHL.specific.loop.csv > LJHL.specific.loop.sort.csv
hicMergeLoops --inputFiles LJGR.shared.loop.sort.csv LJHL.shared.loop.sort.csv LJTZ.shared.loop.sort.csv LJYP.shared.loop.sort.csv --outFileName LJGR.LJTZ.LJYP.LJHL.merge.shared.csv --lowestResolution 25000
~/OLDISK/software/pgltools/sh/pgltools intersect -a LJGR.LJTZ.LJYP.LJHL.merge.shared.sort.csv -b ~/OLDISK/genome_assembly/hic_explorer/4_hicFindTADs/10k/LJGR1_hic_10k_corrected_domains.bed.sort.bed -d 25000 | sort | uniq | wc -l
cat LJYP.shared.loop.sort.csv | awk '{print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3}' > LJYP.shared.loop.sort.rev.csv
cat LJTZ.shared.loop.sort.csv | awk '{print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3}' > LJTZ.shared.loop.sort.rev.csv
cat LJHL.shared.loop.sort.csv | awk '{print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3}' > LJHL.shared.loop.sort.rev.csv
cat LJGR.shared.loop.sort.csv | awk '{print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3}' > LJGR.shared.loop.sort.rev.csv
cat LJGR.specific.loop.sort.csv | awk '{print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3}' > LJGR.specific.loop.sort.rev.csv
cat LJHL.specific.loop.sort.csv | awk '{print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3}' > LJHL.specific.loop.sort.rev.csv
cat LJTZ.specific.loop.sort.csv | awk '{print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3}' > LJTZ.specific.loop.sort.rev.csv
cat LJYP.specific.loop.sort.csv | awk '{print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3}' > LJYP.specific.loop.sort.rev.csv
cat LJGR.specific.loop.csv LJTZ.specific.loop.csv LJYP.specific.loop.csv LJHL.specific.loop.csv > LJGR_LJTZ_LJYP_LJHL.all.specific.csv
~/OLDISK/software/pgltools/sh/pgltools sort LJGR_LJTZ_LJYP_LJHL.all.specific.csv > LJGR_LJTZ_LJYP_LJHL.all.specific.sort.csv
cat LJYP.shared.loop.anchor.tau.bed | awk '{print "YP""\t""Shared""\t"$2"\n""YP""\t""Shared""\t"$3}' | grep -v \-  > LJYP.shared.loop.anchor.tau.bed.bed
cat LJYP.specific.loop.anchor.tau.bed | awk '{print "YP""\t""Specific""\t"$2"\n""YP""\t""Specific""\t"$3}' | grep -v \-  > LJYP.specific.loop.anchor.tau.bed.bed
cat LJHL.specific.loop.anchor.tau.bed | awk '{print "HL""\t""Specific""\t"$2"\n""HL""\t""Specific""\t"$3}' | grep -v \-  > LJHL.specific.loop.anchor.tau.bed.bed
cat LJHL.shared.loop.anchor.tau.bed | awk '{print "HL""\t""Shared""\t"$2"\n""HL""\t""Shared""\t"$3}' | grep -v \-  > LJHL.shared.loop.anchor.tau.bed.bed
cat LJGR.shared.loop.anchor.tau.bed | awk '{print "GR""\t""Shared""\t"$2"\n""GR""\t""Shared""\t"$3}' | grep -v \-  > LJGR.shared.loop.anchor.tau.bed.bed
cat LJGR.specific.loop.anchor.tau.bed | awk '{print "GR""\t""Specific""\t"$2"\n""GR""\t""Specific""\t"$3}' | grep -v \-  > LJGR.specific.loop.anchor.tau.bed.bed
cat LJTZ.specific.loop.anchor.tau.bed | awk '{print "TZ""\t""Specific""\t"$2"\n""TZ""\t""Specific""\t"$3}' | grep -v \-  > LJTZ.specific.loop.anchor.tau.bed.bed
cat LJTZ.shared.loop.anchor.tau.bed | awk '{print "TZ""\t""Shared""\t"$2"\n""TZ""\t""Shared""\t"$3}' | grep -v \-  > LJTZ.shared.loop.anchor.tau.bed.bed
cat LJGR_LJTZ_LJYP_LJHL.loop.anchor.tau.shared.bed | awk '{print "All""\t""Shared""\t"$2"\n""All""\t""Shared""\t"$3}' | grep -v \-  >  LJGR_LJTZ_LJYP_LJHL.loop.anchor.tau.shared.bed.bed
cat LJGR_LJTZ_LJYP_LJHL.loop.anchor.tau.specific.bed | awk '{print "All""\t""Specific""\t"$2"\n""All""\t""Specific""\t"$3}' | grep -v \-  >  LJGR_LJTZ_LJYP_LJHL.loop.anchor.tau.specific.bed.bed
cat LJYP.shared.loop.anchor.tau.bed.bed LJYP.specific.loop.anchor.tau.bed.bed LJHL.specific.loop.anchor.tau.bed.bed LJHL.shared.loop.anchor.tau.bed.bed LJGR.shared.loop.anchor.tau.bed.bed LJGR.specific.loop.anchor.tau.bed.bed LJTZ.specific.loop.anchor.tau.bed.bed LJTZ.shared.loop.anchor.tau.bed.bed LJGR_LJTZ_LJYP_LJHL.loop.anchor.tau.shared.bed.bed LJGR_LJTZ_LJYP_LJHL.loop.anchor.tau.specific.bed.bed > Figure6.Final.tau.loop.tissues.bed
cat LJGR.LJTZ.LJYP.LJHL.merge.shared.sort.csv LJGR_LJTZ_LJYP_LJHL.all.specific.sort.csv > Union.all.csv
~/OLDISK/software/pgltools/sh/pgltools sort Union.all.csv > Union.all.sorted.csv

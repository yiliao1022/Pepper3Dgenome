for i in CA59 HQ RH89A RH89B;do sort -k6,6 -k8,8n SL4.$i.long.paf > SL4.$i.long.sort.paf;done
for i in CA59 HQ RH89A RH89B;do ~/Pepper/PopASSYSV/utilities/k8 ~/Softwares/minimap2/misc/paftools.js view -f maf SL4.$i.long.sort.paf | sed 's/^a /a score=/g' > SL4.
$i.long.sort.maf;done
for i in CA59 HQ RH89A RH89B;do ~/Softwares/linux.x86_64/mafToAxt SL4.$i.long.sort.maf SL4 $i SL4.$i.axt;done
sed -i 's/s E/s HQ.E/g' SL4.HQ.long.sort.maf
~/Softwares/linux.x86_64/mafToAxt SL4.HQ.long.sort.maf SL4 HQ SL4.HQ.axt
~/Softwares/linux.x86_64/mafToAxt SL4.RH89A.long.sort.maf SL4 RH89 SL4.RH89A.axt
~/Softwares/linux.x86_64/mafToAxt SL4.RH89B.long.sort.maf SL4 RH89 SL4.RH89B.axt
cat RH89B.chain.filter.tnet.synnet | grep 'net\|top\|inv\|syn\|non' | awk '$3>50000' > RH89B.chain.filter.tnet.synnet.bed
cat HQ.chain.filter.tnet.synnet | grep 'net\|top\|inv\|syn\|non' | awk '$3>50000' > HQ.chain.filter.tnet.synnet.bed
perl ../../GetSynBreaks.pl  RH89B.chain.filter.tnet.synnet.bed
perl ../../GetSynBreaks.pl  HQ.chain.filter.tnet.synnet.bed

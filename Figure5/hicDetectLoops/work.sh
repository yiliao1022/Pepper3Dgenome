for i in 10k 15k 20k 25k;do cat LJGR.all.sum.$i.loop.8M.bed | awk '{print "LJGR""\t"$3-$2"\t"$5-$3}' > LJGR.all.sum.$i.loop.8M.bed.bed;done
for i in 10k 15k 20k 25k;do cat LJHL.all.sum.$i.loop.8M.bed | awk '{print "LJHL""\t"$3-$2"\t"$5-$3}' > LJHL.all.sum.$i.loop.8M.bed.bed;done
for i in 10k 15k 20k 25k;do cat LJTZ.all.sum.$i.loop.8M.bed | awk '{print "LJTZ""\t"$3-$2"\t"$5-$3}' > LJTZ.all.sum.$i.loop.8M.bed.bed;done
for i in 10k 15k 20k 25k;do cat LJYP.all.sum.$i.loop.8M.bed | awk '{print "LJYP""\t"$3-$2"\t"$5-$3}' > LJYP.all.sum.$i.loop.8M.bed.bed;done
cat *.bed.bed >> All.bed

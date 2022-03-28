hicMergeLoops --inputFiles LJHL.all.sum.10k.loop.8M.bed LJHL.all.sum.15k.loop.8M.bed LJHL.all.sum.20k.loop.8M.bed LJHL.all.sum.25k.loop.8M.bed --outFileName LJHL.10_15_20_25.merged.csv --lowestResolution 25000
~/OLDISK/software/pgltools/sh/pgltools sort LJHL.10_15_20_25.merged.csv > LJHL.10_15_20_25.merged.sort.csv

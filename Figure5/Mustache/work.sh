cat ../LJGRsum/Mustache/All.combined.loop.bed.sort.merge.d0.bed ../LJTZsum/Mustache/All.combined.loop.bed.sort.merge.d0.bed.sort.bed ../LJHLsum/Mustache/Four.combined.sort.merge.csv ../LJYPsum/Mustache/Four.combined.sort.merge.csv | grep -v \# > Four.tissues.loops.combined.csv
~/OLDISK/software/pgltools/sh/pgltools sort Four.tissues.loops.combined.csv > Four.tissues.loops.combined.sort.csv
~/OLDISK/software/pgltools/sh/pgltools merge -a Four.tissues.loops.combined.sort.csv > Four.tissues.loops.combined.sort.merge.csv

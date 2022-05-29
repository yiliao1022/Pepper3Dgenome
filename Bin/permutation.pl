#!/usr/bin/perl
use strict;
use warnings;

my $i;

for ($i=1;$i<=100;$i++) {
system "bedtools shuffle -i $ARGV[0] -g $ARGV[1] -chrom -seed $i >$ARGV[0].shuffle.$i.bed";
system "perl Peak.ins.pl $ARGV[0].shuffle.$i.bed $ARGV[2]";
system "cat $ARGV[2].peak.breatpoints.bed >> $ARGV[2].background.bed";
system "rm $ARGV[0].shuffle.$i.bed $ARGV[2].peak.breatpoints.bed";  
}

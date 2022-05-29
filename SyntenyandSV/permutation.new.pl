#!/usr/bin/perl

use strict;
use warnings;

my $i;

for ($i=1;$i<=100;$i++) {
system "bedtools shuffle -i $ARGV[0] -g $ARGV[1] -chrom -seed $i >$ARGV[0].shuffle.$i.bed";
system "perl break2bin.pl $ARGV[0].shuffle.$i.bed $ARGV[2] >> $ARGV[2].background.bed";
system "perl break2bin.pl $ARGV[0].shuffle.$i.bed $ARGV[3] >> $ARGV[3].background.bed";
system "perl break2bin.pl $ARGV[0].shuffle.$i.bed $ARGV[4] >> $ARGV[4].background.bed";
system "perl break2bin.pl $ARGV[0].shuffle.$i.bed $ARGV[5] >> $ARGV[5].background.bed";

system "perl Scale.pl $ARGV[0].shuffle.$i.bed $ARGV[6]";
system "rm *_K.bed";
system "perl Scale.pl $ARGV[0].shuffle.$i.bed $ARGV[7]";
system "rm *_K.bed";
system "perl Scale.pl $ARGV[0].shuffle.$i.bed $ARGV[8]";
system "rm *_K.bed";
system "perl Scale.pl $ARGV[0].shuffle.$i.bed $ARGV[9]";
system "rm *_K.bed";

system "rm $ARGV[0].shuffle.$i.bed";  
}

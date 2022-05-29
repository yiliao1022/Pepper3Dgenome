#!/usr/bin/perl

use strict;
use warnings;

my $i;

for ($i=0;$i<20;$i++) {
 open In, "$ARGV[0]" or die "$!";
 open Out, ">$ARGV[0].flank.$i\_K.bed" or die "$!";
  
  while (<In>) {
   chomp;
   my @temp = split(/\t/,$_);
   my $len = $temp[2]-$temp[1];
   my $beg0 = int ($temp[1]-($len/2));
   my $bin = int ($len/10);
       
   my $start = $beg0+$bin*$i;
   my $end   = $beg0+($i+1)*$bin;;
  
    if ($start<0) {
         $start =0;
         $end = $bin;
      } 
    if ($end <0) {
         $start =0;
         $end = $bin;
     }
     print Out "$temp[0]\t$start\t$end\n";
   }
 print "$i:"; 
 system "cat $ARGV[0].flank.$i\_K.bed | awk '{sum+=(\$3-\$2+1)} END {print sum}'";
 print "\n";
 system "bedtools intersect -a $ARGV[1] -b $ARGV[0].flank.$i\_K.bed | awk '{sum+=(\$3-\$2+1)} END {print sum}' >> $ARGV[1].peak.cov.bed";
 # system "echo "\n" >> $ARGV[1].peak.cov.bed";
}


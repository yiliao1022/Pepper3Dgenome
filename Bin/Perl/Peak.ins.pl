#!/usr/bin/perl

use strict;
use warnings;

my $i;
for($i=-100;$i<101;$i++) {
   open In, "$ARGV[0]" or die "$!";
   open Out, ">$ARGV[0].flank.$i\_K.bed" or die "$!";
   
 while (<In>) {
   chomp;
   my @temp = split(/\t/,$_);
   my $start = $temp[1]+20000*$i;
   my $end   = $temp[2]+20000*$i;
  
    if ($start<0) {
         $start =0;
      } 
    if ($end <0) {
         $end =40000;
     }
  print Out "$temp[0]\t$start\t$end\n";
 }

  system "bedtools intersect -a $ARGV[1] -b $ARGV[0].flank.$i\_K.bed -u | wc -l >> $ARGV[1].peak.breatpoints.bed";
  system "rm $ARGV[0].flank.$i\_K.bed";
}





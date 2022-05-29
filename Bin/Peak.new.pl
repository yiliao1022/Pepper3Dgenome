#!/usr/bin/perl

use strict;
use warnings;

my %hash;

open SIZE, "$ARGV[2]" or die "$!";
while (<SIZE>) {
chomp;
my @tmp = split(/\t/,$_);
$hash{$tmp[0]} = $tmp[1];
}


my $i;
my $m=0;
my $n=0;
for($i=-100;$i<101;$i++) {
   open In, "$ARGV[0]" or die "$!";
   open Out, ">$ARGV[0].flank.$i\_K.bed" or die "$!";
   
 while (<In>) {
   chomp;
   my @temp = split(/\t/,$_);
   my $start = $temp[1]+5000*$i;
   my $end   = $temp[2]+5000*$i;
  
  
    if ($start<0) {
         $n++;
         $start =$n*10000;
         $end = 10000+$n*10000;
    } 
    if ($end > $hash{$temp[0]}) {
         $m++;
         $start = $hash{$temp[0]} - 10000*$m;
         $end =$hash{$temp[0]} -10000*$m+10000;
     }
       
  print Out "$temp[0]\t$start\t$end\n"; 
  
}

  system "bedtools intersect -a $ARGV[1] -b $ARGV[0].flank.$i\_K.bed -u | wc -l >> $ARGV[1].peak.breakponits.bed";
  system "rm $ARGV[1].$ARGV[0].flank.$i\_K.bed";
}

system "paste cor.bed $ARGV[1].peak.breakpoints.bed > $ARGV[1].peak.breakpoints.bed.bed";

#!/usr/bin/perl
use strict;
use warnings;

open In, "$ARGV[0]" or die "$!";
open Out, ">$ARGV[0].sym.mat" or die "$!";
my @mat;
my $n = $ARGV[1];
my $bin = $ARGV[2];

my $i;
my $j;

for ($i=0;$i<$n;$i++) {
   for ($j=0;$j<$n;$j++) {
     $mat[$i][$j]=0;
     }
}

while (<In>) {
chomp;
my @temp = split(/\t/,$_);
my $m = $temp[1]/$bin;
my $k = $temp[4]/$bin;
$mat[$m][$k]=$temp[6];
$mat[$k][$m]=$temp[6];
}


my $row;
my $col;

for ($row=0;$row<$n;$row++) {
  for ($col=0;$col<$n;$col++) {
       print Out "$mat[$row][$col]\t";
   }
   print Out "\n";
}  

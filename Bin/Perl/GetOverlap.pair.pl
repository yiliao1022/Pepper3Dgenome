#!/bin/perl
use strict;
use warnings;

my @temp;
open In, "$ARGV[0]" or die "$!";
while (<In>) {
chomp;
push (@temp,$_);
}

my $i;
my $j;

for ($i=0;$i<=$#temp;$i++){
 for ($j=0;$j<=$#temp;$j++) {
  my $inum = `cat $temp[$i] | wc -l`;
  my $jnum = `cat $temp[$j] | wc -l`;
  my $tmp = `~/OLDISK/software/bedtools intersect -wa -wb -a $temp[$i] -b $temp[$j] | wc -l`;
  my $fre = $tmp/($inum);
  print "$fre\t";
 }
  print "\n";
}

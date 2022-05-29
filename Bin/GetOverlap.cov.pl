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
  my $inum = `cat $temp[$i] | awk '{sum+=(\$3-\$2+1)}END{print sum}'`;
  my $jnum = `cat $temp[$j] | awk '{sum+=(\$3-\$2+1)}END{print sum}'`;
  my $tmp = `~/OLDISK/software/bedtools intersect -wa -wb -a $temp[$i] -b $temp[$j] -f 0.8 -F 0.8 | awk '{sum+=(\$3-\$2+1)}END{print sum}'`;
  my $fre = $tmp/($inum);
  print "$fre\t" unless ($j=$#temp);
 }
  print "\n";
}

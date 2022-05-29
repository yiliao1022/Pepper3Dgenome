#!/usr/bin/perl
use strict;
use warnings;


open In, "$ARGV[0]" or die "$!";
open OUT, ">$ARGV[0].out" or die "$!";

while (<In>) {
chomp;
my @temp = split(/\t/,$_);
if (($temp[2]-$temp[1])/100000 == 1) {
print OUT "$_\n";
} else {
my $n = ($temp[2]-$temp[1])/100000;
my $i;
for ($i=0;$i<$n;$i++) {
my $beg = $temp[1] + 100000*$i;
my $end = $temp[1] + 100000*($i+1);
print OUT "$temp[0]\t$beg\t$end\t$temp[3]\n";
  }
 }
}

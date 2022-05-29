#!/usr/bin/perl 

use strict;
use warnings;


open In, "$ARGV[0]" or die "$!";
open OUT, ">$ARGV[0].out" or die "$!";

while (<In>) {
next if ($_=~/genename/);
chomp;
my @temp = split (/\t/,$_);

my $LJGR = 0.1 + ($temp[1]+$temp[2]+$temp[3])/3;
my $LJYP = 0.1 + ($temp[4]+$temp[5]+$temp[6])/3;
my $root = 0.1 + ($temp[7]+$temp[8]+$temp[9])/3;
my $LJTZ = 0.1 + ($temp[10]+$temp[11]+$temp[12])/3;
my $LJHL = 0.1 + ($temp[13]+$temp[14]+$temp[15])/3;


if ($temp[13]> 1 && $temp[14] > 1 && $temp[15] > 1 && $LJHL/$LJGR > 3 && $LJHL/$root > 3 && $LJHL/$LJTZ > 3 && $LJHL/$LJYP > 3) {

 print OUT "$_\n";
} 


}

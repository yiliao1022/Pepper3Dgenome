#!/usr/bin/perl

use strict;
use warnings;

open In, "$ARGV[0]" or die "$!";
open Out, ">$ARGV[0].bed" or die "$!";

print Out "LJGR1\tLJGR2\tLJGR3\tLJHL1\tLJHL2\tLJHL3\tLJTZ1\tLJTZ2\tLJTZ3\tLJYP1\tLJYP2\tLJYP3\tRoot1\tRoot2\tRoot3\n";
while (<In>) {
chomp;
my @temp = split(/\t/,$_);
my $n=$temp[1]/$ARGV[1];
my $ID=join("_",($temp[0],$n));
print Out "$ID\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\t$temp[8]\t$temp[9]\t$temp[10]\t$temp[11]\t$temp[12]\t$temp[13]\t$temp[14]\t$temp[15]\t$temp[16]\t$temp[17]\n";
}

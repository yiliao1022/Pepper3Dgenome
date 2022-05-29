#!/usr/bin/perl
use strict;
use warnings;

open In, "$ARGV[0]" or die "$!";
open ALL, ">$ARGV[0].all.out" or die "$!";
open MB, ">$ARGV[0].2m.out" or die "$!";
print ALL "GR1\tGR2\tHL1\tHL2\tTZ1\tTZ2\tYP1\tYP2\n";
print MB "GR1\tGR2\tHL1\tHL2\tTZ1\tTZ2\tYP1\tYP2\n";

while (<In>) {
chomp;
my @temp = split (/\t/,$_);
if ($temp[2]<2000 and $temp[3]<2000 and $temp[4]<2000 and $temp[5]<2000 and $temp[6]<2000 and $temp[7]<2000 and $temp[8]<2000 and $temp[9]<2000) {
#print ALL "$temp[0]\_$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\t$temp[8]\t$temp[9]\n";
} else {
print ALL "$temp[0]\_$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\t$temp[8]\t$temp[9]\n";
if ($temp[1]-$temp[0]>2000000 and $temp[1]-$temp[0]<20000000) {
print MB "$temp[0]\_$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\t$temp[8]\t$temp[9]\n";
}

} #elsif ($temp[1]-$temp[0]>2000000 and $temp[1]-$temp[0]<20000000) {
#print MB "$temp[0]\_$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\t$temp[8]\t$temp[9]\n";
# }
}

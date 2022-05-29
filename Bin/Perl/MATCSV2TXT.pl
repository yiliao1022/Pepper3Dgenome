#!/usr/bin/perl
use strict;
use warnings;

open In, "$ARGV[0]" or die "$!";
open Out, ">$ARGV[0].txt" or die "$!";
my $bin = $ARGV[1];

while (<In>) {
next if ($_=~/V1/);
chomp;
my @temp = split (/,/,$_);
$temp[0]=~s/\"//g;
my $i;
for ($i=$temp[0];$i<=$#temp;$i++) {
my $entry = 1000000*(10**$temp[$i]);
my $beg = ($temp[0]-1)*$bin;
my $end = ($i-1)*$bin;
print Out "$beg\t$end\t$entry\n"
 }
}

#!/usr/bin/perl
use strict;
use warnings;

open In, "$ARGV[0]" or die "$!";
open Out, ">$ARGV[0].TopDom.matrix" or die "$!";
my $bin = $ARGV[1];
my $chr = $ARGV[2];
while (<In>) {
next if ($_=~/V1/);
chomp;
my @temp = split (/,/,$_);
$temp[0]=~s/\"//g;
my $beg = ($temp[0]-1)*$bin;
my $end = ($temp[0])*$bin;
print Out "$chr\t$beg\t$end\t";
my $i;
for ($i=1;$i<$#temp;$i++) {
my $entry = 1000000*(10**$temp[$i]);
print Out "$entry\t"
}
my $end_entry = 1000000*(10**$temp[$#temp]);
print Out "$end_entry\n";
}

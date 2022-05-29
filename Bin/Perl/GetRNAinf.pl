#!/bin/perl
use strict;
use warnings;

my %hash;

open IN, "$ARGV[0]" or die "$!";
open OUT, ">$ARGV[0].RNA.out" or die "$!";

while (<IN>) {
next if ($_=~/logFC/);
my @temp = split (/\t/,$_);
my @unit = split (/_/,$temp[0]);
my $beg = $unit[1]*40000;
my $end = $beg + 40000;
my $logFC= abs($temp[1]);
print OUT "$unit[0]\t$beg\t$end\t$logFC\t$temp[5]\n"
}

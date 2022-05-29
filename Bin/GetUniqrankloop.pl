#!/bin/perl
use strict;
use warnings;

my %hash;

open IN, "$ARGV[0]" or die "$!";
open OUT, ">$ARGV[0].out" or die "$!";

while (<IN>) {
next if ($_=~/logFC/);
my @temp = split (/\t/,$_);
my @unit = split (/_/,$temp[0]);
my $anchor1 = join ("_",($unit[0],$unit[1]));
my $anchor2 = join ("_",($unit[2],$unit[3]));
if (exists $hash{$anchor1}) {
} else {
my $end1 = $unit[1]+40000;
print OUT "$unit[0]\t$unit[1]\t$end1\n"
}

if (exists $hash{$anchor2}) {
} else {
my $end2 = $unit[3]+40000;
print OUT "$unit[0]\t$unit[3]\t$end2\n"
}

}

#!/bin/perl
use strict;
use warnings;
my %hash;
open In, "$ARGV[0]" or die "$!";
while (<In>) {
chomp;
$hash{$_}=1;
}
open GENE, "$ARGV[1]" or die "$!";
while (<GENE>) {
chomp;
my @temp = split (/\t/,$_);
if (exists $hash{$temp[0]}) {
  print "$temp[16]\t$temp[17]\t$temp[18]\n";
}
}

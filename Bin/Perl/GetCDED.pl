#!/usr/bin/perl
use strict;
use warnings;

my %hash;

open In, "$ARGV[0]" or die "$!";

while (<In>) {
next if ($_=~/logFC/);
my @temp = split (/\t/,$_);
my @unit = split (/_/,$temp[0]);
my $bin1 = $unit[7]/40000;
my $bin2 = $unit[9]/40000;
my $key1 = join ("_",($unit[6],$bin1));
my $key2 = join ("_",($unit[8],$bin2));

$hash{$key1} = $temp[1];
$hash{$key2} = $temp[1];
}


open In1, "$ARGV[1]" or die "$!";

while (<In1>) {
next if ($_=~/logFC/);
my @temp1 = split (/\t/,$_);
if (exists $hash{$temp1[0]}) {
print "$temp1[0]\t$hash{$temp1[0]}\t$temp1[1]\n";
}
}


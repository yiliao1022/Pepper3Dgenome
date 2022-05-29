#!/usr/bin/perl
use strict;
use warnings;

my $short = 0;
my $medium = 0;
my $long =0;
my $inter=0;
open In, "$ARGV[0]" or die "$!";

while (<In>) {
chomp;
my @temp = split (/\t/,$_);
if ($temp[0] eq $temp[3]) {
if ($temp[4]-$temp[1]<=2000000) {
$short+=$temp[6];
} elsif ($temp[4]-$temp[1]>10000000) {
$long+=$temp[6];
} else {
$medium++;
} } else {
$inter++
}
}

print "Short: $short\nMedium: $medium\nLong: $long\nInterChr: $inter\n";

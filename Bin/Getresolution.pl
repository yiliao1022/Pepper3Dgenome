#!/bin/perl 

use strict;
use warnings;

my %hash;

open In, "$ARGV[0]" or die "$!";

while (<In>) {
my @temp=split (/\t/,$_);
my $bin1 = join ("_",($temp[0],$temp[1]));
my $bin2 = join ("_",($temp[3],$temp[4]));
if(exists $hash{$bin1} and ($bin1!~$bin2)) {
$hash{$bin1}+=$temp[6];
} else {
if ($bin1!~$bin2) {
$hash{$bin1}=$temp[6];
}
}
}
my $n;
my $m;
foreach my $key (keys %hash) {
$n++;
if ($hash{$key}>=1000) {
$m++
}
}

print "$m/$n\n";

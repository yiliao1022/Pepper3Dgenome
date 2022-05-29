#!/usr/bin/perl
use strict;
use warnings;

my %hash;

open In, "$ARGV[0]" or die "$!";

while (<In>) {
chomp;
my @temp = split (/\t|\s/,$_);
my $beg1 = $temp[1]-20000;
my $beg2 = $temp[1]+20000;
my $end1 = $temp[2]-20000;
my $end2 = $temp[2]+20000;

my $beg = join ("_",($temp[0],$beg1,$beg2));
my $end = join ("_",($temp[0],$end1,$end2));

if (exists $hash{$beg}) {
print "$beg\t$hash{$beg}\_$temp[3]\n";
} 

$hash{$end} = $temp[3];

}

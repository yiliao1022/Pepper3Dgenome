#!/bin/perl

use strict;
use warnings;

my %hash1;
my %hash2;

open In1, "$ARGV[0]" or die "$!";

while (<In1>) {
chomp;
my @temp = split (/\t/,$_);
my $key = join ("_",($temp[0],$temp[1],$temp[2],$temp[3],$temp[4],$temp[5]));
if (exists $hash1{$key}) {
  if ($hash1{$key}=-0.1) {
$hash1{$key} = $temp[9];
print "1: haha!\n";
  } 
} else {
$hash1{$key} = $temp[9];
 
}
}
open In2, "$ARGV[1]" or die "$!";

while (<In2>) {
chomp;
my @temp1 = split (/\t/,$_);
my $key1 = join ("_",($temp1[3],$temp1[4],$temp1[5],$temp1[0],$temp1[1],$temp1[2]));
if (exists $hash2{$key1}) {
  if ($hash2{$key1}=-0.1) {
$hash2{$key1} = $temp1[9];
print "haha!\n";
} 
} else {
$hash2{$key1} = $temp1[9];
}
}

open OUT, ">$ARGV[2]" or die "$!";
foreach my $loop (keys %hash1) {
if (exists $hash2{$loop}) {
print OUT "$loop\t$hash1{$loop}\t$hash2{$loop}\n";
} else {
print OUT "$loop\t$hash1{$loop}\t-0.1\n";

}

}

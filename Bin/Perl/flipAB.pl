#!/usr/bin/perl

use strict;
use warnings;

open In, "$ARGV[0]" or die "$!";
my %hash = (
             Ca_59Chr01 => 1, Ca_59Chr02 => 2, Ca_59Chr05 =>3, Ca_59Chr08=>4, Ca_59Chr09=>5, Ca_59Chr10=>6
            );

while (<In>) {
chomp;
my @temp = split (/\t/,$_);
if (exists $hash{$temp[0]}) {
my $value;
   if ($temp[3]<=0) {
      $value = abs ($temp[3]);
   } else {
       $value = $temp[3]-2*$temp[3];
   }
print "$temp[0]\t$temp[1]\t$temp[2]\t$value\n";
} else {
print "$_\n";;
}
}

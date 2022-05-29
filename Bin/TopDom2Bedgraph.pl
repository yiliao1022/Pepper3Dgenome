#!/usr/bin/perl
use strict;
use warnings;

open In, "$ARGV[0]" or die "$!";
open Out, ">$ARGV[0].bedgraph" or die "$!";
my $n=0;
while (<In>) {
chomp;
$n++;
my @temp = split (/\t/,$_);
next if ($_=~/from/);
my $m = int ($n/2);

if ($m==$n/2) {
print Out "$temp[0]\t$temp[2]\t$temp[4]\tdomain\t$temp[1]\t\.\t$temp[2]\t$temp[4]\t\#FF4848\n";
} else {
print Out "$temp[0]\t$temp[2]\t$temp[4]\tdomain\t$temp[1]\t\.\t$temp[2]\t$temp[4]\t\#4848FF\n";

}


}

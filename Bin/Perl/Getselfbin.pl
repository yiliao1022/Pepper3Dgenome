#!/bin/perl 
use strict;
use warnings;

open In, "$ARGV[0]" or die "$!";
while (<In>) {
my @temp = split (/\t/,$_);
my @unit = split (/_/,$temp[0]);
if ( ($unit[0] eq $unit[2]) and ($unit[1]==$unit[3])) {
print "$_";
}
}

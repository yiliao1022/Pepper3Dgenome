#!/bin/perl 

use strict;
use warnings;

my %hash;
open In, "$ARGV[0]" or die "$!";

while (<In>) {
next if ($_=~/log/);
my @temp = split (/\t/,$_);
#my @unit = split (/_/,$temp[0]);
my $key = join("_",($temp[0],$temp[1],$temp[2],$temp[3],$temp[4],$temp[5]));
if (!exists $hash{$key}) {
print "$_";
$hash{$key}=0;
}else {
}
}

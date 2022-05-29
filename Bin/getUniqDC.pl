#!/bin/perl 

use strict;
use warnings;

my %hash;
open In, "$ARGV[0]" or die "$!";

while (<In>) {
next if ($_=~/log/);
my @temp = split (/\t/,$_);
my $key1 = join("_",($temp[0],$temp[1]));
my $key2 = join("_",($temp[2],$temp[3]));
if (!exists $hash{$key1} or !exists $hash{$key2}) {
print "$_";
$hash{$key1}=0;
$hash{$key2}=0;
}else {
}
}

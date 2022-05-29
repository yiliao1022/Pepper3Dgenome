#!/usr/bin/perl
use strict;
use warnings;

open In, "$ARGV[0]" or die "$!";
open Out, ">$ARGV[0].bed" or die "$!";
while (<In>) {
chomp;
my @temp = split (/\t/,$_);
my $value = &log10($temp[3]+1);
print Out "$temp[0]\t$temp[1]\t$temp[2]\t$value\n";
}

sub log10 
{
    my $n = shift;     
    # using pre-defined log function
    return log($n) / log(10);
}

#!/usr/bin/perl
use strict;
use warnings;

open In, "$ARGV[0]" or die "$!";
open Out, ">$ARGV[0].bed" or die "$!";
my $n=1;

while (<In>) {
    chmod;
    my @temp = split(/\s/,$_);
	my $is = $n % 201;
	print "$is\n";    
    if ($is != 0) {
        print Out "$temp[0]\t";
    } else {
        print Out "$temp[0]\n";
    }
$n++;
}

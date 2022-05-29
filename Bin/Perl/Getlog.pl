#!/usr/bin/perl

use strict;
use warnings;

open In, "$ARGV[0]" or die "$!";
open Out, ">$ARGV[0].bed" or die "$!";

while (<In>) { 
  my $log;
  my @temp = split (/\t/,$_);
  if ($temp[3]<=10) {
  $log = 0;
  } else {
  $log = &log10($temp[3]);
  }
print Out "$temp[0]\t$temp[1]\t$temp[2]\t$log\n";
}
sub log10 {
    my $n = shift;
    return log($n)/log(10);
}

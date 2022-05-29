#!/usr/bin/perl
use strict;
use warnings;


my $fa = &readSeq($ARGV[0]);


open In, "$ARGV[1]" or die "$!";
open Out, ">$ARGV[1].out" or die "$!";

while (<In>) {
chomp;
my @temp = split (/\t/,$_);
my $beg = $temp[2]-10;
my $end = $temp[3];
my $five = substr $$fa{$temp[0]}, $beg, 10;
my $three = substr $$fa{$temp[0]}, $end, 10;


print Out "$temp[0]\t$temp[1]\t$five\t$three\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\t$temp[8]\t$temp[9]\t$temp[10]\t$temp[11]\n";
}







sub readSeq {
my $seq = shift;
my %hash;
local $/="\n>";
open FA, "$seq" or die "$!";
while (my $sequence = <FA>) {
         chomp $sequence;
         my ($id) = $sequence =~ /^>*(\S+)/;  # parse ID as first word in FASTA header
         $sequence =~ s/^>*.+\n//;  # remove FASTA header
         $sequence =~ s/\n//g;  # remove endlines
         $hash{$id}=$sequence;
    }
my $ref=\%hash;
return $ref;
}



#!/usr/bin/perl
use strict;
use warnings;

my $fa = &readSeq($ARGV[0]);

$/="###\n";
open In, "$ARGV[1]" or die "$!";
open Out, ">$ARGV[1].ltr.fa" or die "$!";

while (<In>) {
next if ($_=~/^\s*$/);
my @temp = split (/\n/,$_);
my $seq1;
my $seq2;
my $seq3;
my $chr;
my $coor;
my $orient;
if ($temp[0]=~/tsd\=NA/) {
my @unit1 = split (/\t/,$temp[1]);
$chr=$unit1[0];
$coor=$unit1[3];
$orient=$unit1[6];
my $len1 = $unit1[4] - $unit1[3] + 1;
$seq1 = substr $$fa{$unit1[0]}, $unit1[3]-1, $len1;

my $beg = $unit1[4]+1;
my @unit2 = split (/\t/,$temp[3]);
my $len2 = $unit2[4] - $unit2[3] + 1;
$seq2 = substr $$fa{$unit1[0]}, $unit2[3]-1, $len2;

my $end = $unit2[3]-1;
my $len3 = $end - $beg + 1;
$seq3 = substr $$fa{$unit1[0]}, $beg-1, $len3;

} else {
my @unit1 = split (/\t/,$temp[2]);
$chr=$unit1[0];
$coor=$unit1[3];
$orient=$unit1[6];
my $len1 = $unit1[4] - $unit1[3] + 1;
$seq1 = substr $$fa{$unit1[0]}, $unit1[3]-1, $len1;

my $beg = $unit1[4]+1;
my @unit2 = split (/\t/,$temp[4]);
my $len2 = $unit2[4] - $unit2[3] + 1;
$seq2 = substr $$fa{$unit1[0]}, $unit2[3]-1, $len2;

my $end = $unit2[3]-1;
my $len3 = $end - $beg + 1;
$seq3 = substr $$fa{$unit1[0]}, $beg-1, $len3;
}





if ($orient=~/\-/) {
print "reverse\n";
my $revcomp1 = reverse $seq1;
my $revcomp2 = reverse $seq2;
my $revcomp3 = reverse $seq3;
$revcomp1 =~ tr/ATGCatgc/TACGtacg/;
$revcomp2 =~ tr/ATGCatgc/TACGtacg/;
$revcomp3 =~ tr/ATGCatgc/TACGtacg/;
print Out ">$chr\_$coor\_5LTR\n$revcomp1\n";
print Out ">$chr\_$coor\_3LTR\n$revcomp2\n";
print Out ">$chr\_$coor\_int\n$revcomp3\n";
} else {
print Out ">$chr\_$coor\_5LTR\n$seq1\n";
print Out ">$chr\_$coor\_3LTR\n$seq2\n";
print Out ">$chr\_$coor\_int\n$seq3\n";
 }
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




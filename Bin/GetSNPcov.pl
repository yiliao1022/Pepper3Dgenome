#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);

my ($ref,$query,$reffolder,$singfolder,$bed,$maf2fasta,$help);

# establish parameters
GetOptions( 'ref=s'=>\$ref,
            'query=s'=>\$query,
            'reffolder=s'=>\$reffolder,
            'sing=s'=>\$singfolder,
            'bed=s'=>\$bed,
            'maf2fasta=s'=>\$maf2fasta, 
            'help'=>\$help
           );

open In, "$bed" or die "$!";

while (<In>) {
chomp;
my @temp = split (/\t/,$_);
my @unit = split (/\./,$temp[0]);
print "$unit[1]\n";
system "${maf2fasta}maf2fasta $reffolder/$unit[1]/$ref $singfolder/$unit[1]/$ref.$query.sing.maf $temp[1] $temp[2] > $unit[1].$temp[1].$temp[2].fa";
system "perl 9_MAF2VCF_SNPnew.pl -input $unit[1].$temp[1].$temp[2].fa -fasize CA59.sizes -spec pepper -glst glst.bed >> $bed.bed ";
system "rm $unit[1].$temp[1].$temp[2].fa";
}

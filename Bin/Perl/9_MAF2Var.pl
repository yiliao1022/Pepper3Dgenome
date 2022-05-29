#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use List::MoreUtils qw(uniq);

####################################################################################
# MAF2VCF_SNP.pl
#
# Authors: Yi Liao (07/11/2020)
#
# Copyright (C) Not for commercial use
#
# Report SNP from .maf file to .vcf file
#
# Prerequisite : TBA tools; MUSCLE
#
# Usage: perl $0 -sv SV.txt -svtype INS -muscle /path/to/muscle -tba /path/to/tba
# -sv         raw merged SV file                        [REQUIRED]
# -svtype     SV types [INS|DEL|CNV|INV]                [REQUIRED]
# -gseq       path to all genome sequences              [REQUIRED]
# -muscle     path for muscle program
# -tba        path for tba program
# -help       print this information
#
####################################################################################

my ($input,$fasize,$chr,$ref,$query,$glst,$output,$bed,$Help);

GetOptions (
"input:s" =>\$input,
"fasize:s" =>\$fasize,
"ref:s" => \$ref,
"query:s" => \$query,
"glst:s"=>\$glst,
"output:s" =>\$output,
"bed:s"=>\$bed,
"chr:s"=>\$chr,
"help" =>\$Help
);

if ($Help) {
print << "END.";
Useage: perl $0 -input maf2fasta.fa -spec rice/human/fruitfly/maize -ref IRGSP1 -fasize chr.sizes -glst glst.txt -head

-input      : maf2fasta files
-fasize     : reference sizes by chromosome
####-spec       : species name, e.g. rice, human,maize
-ref        : reference name, e.g. hg38, ISO1, IRGSP1
-glst       : individuals in the population, the first one should be the reference
-head       : if defined then print the vcf header eitherwise not 

END.
exit;
}


open GLST, "$glst" or die "$!";
chomp (my @genome=<GLST>);
close GLST;

my ($sequence,$coordinate,$index)=&Parse_MAF2FastaOutput($input);

my $ref1 = shift @genome;
my $a;
my $size = @$coordinate;
my %conserved;
my %missing;
my %SNP;
my %DEL;


foreach my $ing (@genome) {
$conserved{$ing}=0;
$missing{$ing}=0;
$SNP{$ing}=0;
$DEL{$ing}=0;
}


open BED,"$bed" or die "$!";
while (<BED>) {
my @unit = split(/\t/,$_);
next if ($unit[0]!~/$chr/);
for ($a=$unit[1];$a<=$unit[2];$a++) {
   foreach my $ind (@genome) {
     $$sequence{$ref1}->[$$coordinate[$a]] = uc ($$sequence{$ref1}->[$$coordinate[$a]]);
     $$sequence{$ind}->[$$coordinate[$a]] = uc ($$sequence{$ind}->[$$coordinate[$a]]);

     if (($$sequence{$ref1}->[$$coordinate[$a]] eq $$sequence{$ind}->[$$coordinate[$a]]) and $$sequence{$ref1}->[$$coordinate[$a]] !~/N\|-/) {
	     $conserved{$ind}++;
     } elsif ($$sequence{$ref1}->[$$coordinate[$a]]!~/N\|-/ and $$sequence{$ind}->[$$coordinate[$a]] eq "N") {
             $missing{$ind}++;
     } elsif ($$sequence{$ref1}->[$$coordinate[$a]]!~/N\|-/ and $$sequence{$ind}->[$$coordinate[$a]] eq "-") {
             $DEL{$ind}++;
     } elsif ($$sequence{$ref1}->[$$coordinate[$a]]!~/N\|-/ and $$sequence{$ind}->[$$coordinate[$a]] !~/N\|-/ and ($$sequence{$ref1}->[$$coordinate[$a]] ne $$sequence{$ind}->[$$coordinate[$a]])) {
	     $SNP{$ind}++;
     }
   }
}


foreach my $k1 (keys %conserved) {
print "$k1\t$conserved{$k1}\t";
}

foreach my $k2 (keys %missing) {
print "$k2\t$missing{$k2}\t";
}

foreach my $k3 (keys %SNP) {
print "$k3\t$SNP{$k3}\t";
}

foreach my $k4 (keys %DEL) {
print "$k4\t$DEL{$k4}\n";
}

}
##################################
############ Subroutine ##########
##################################

sub Parse_MAF2FastaOutput {
   my ($file) = @_;
   my @data;
   my %hash;
   $/="\n";
   open(DATA, "$file") or die "Unable to open file maf2fasta file, $!";
   while (<DATA>) {
         chomp;
       push (@data, $_);
   }
   close DATA;
## load reference sequence into hash
   my @seq_num = split (/\s/,$data[0]);
   my @ref_name = split (/:/,$data[1]);
   my $ref = $ref_name[0];
   my @ref_seq = split(//,$data[$seq_num[0]+1]);
   $hash{$ref}=\@ref_seq;
 ##
## load query sequence into hash
   my $i;
   for ($i=2;$i<$seq_num[0]+1;$i++) {
   my @query_seq = split(//,$data[$seq_num[0]+$i]);
      $hash{$data[$i]} = \@query_seq;
   }
 ##
## relate the sequence coordinates to Maf2fasta reference coordinates
   my $len_n=1;
   my %index_proxy;
   my $t;
   my @coordinate_proxy;
   for ($t=0;$t<$seq_num[1];$t++) {
       if ($hash{$ref}->[$t] =~/[ATGCatgc]/) {
       $index_proxy{$t}=$len_n;
       $len_n++;
       push (@coordinate_proxy,$t);
       }
   }
##
my $seq=\%hash;
my $cor_proxy=\@coordinate_proxy;
my $index_convert=\%index_proxy;
   return ($seq,$cor_proxy,$index_convert);
}

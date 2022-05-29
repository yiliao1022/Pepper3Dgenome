#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use List::MoreUtils qw(uniq);

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
spec       : species name, e.g. rice, human,maize
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


my %hash;
open SIZE, "$fasize" or die "$!";
while (<SIZE>) {
chomp;
my @tmp = split(/\t/,$_);
$hash{$tmp[0]} = $tmp[1];
}


my $i;
my $m=0;
my $n=0;

for($i=-100;$i<101;$i++) {
   open In, "$bed" or die "$!";
   open Out, ">$bed.flank.$i\_K.bed" or die "$!";
   open Out1, ">$bed.flank.$i\_K.bed.bed" or die "$!";
   
 while (<In>) {
   chomp;
   my @temp = split(/\t/,$_);
   my $start = $temp[1]+5000*$i;
   my $end   = $temp[2]+5000*$i;
  
  
    if ($start<0) {
         $n++;
         $start =$n*40000;
         $end = 40000+$n*40000;
    } 
    if ($end > $hash{$temp[0]}) {
         $m++;
         $start = $hash{$temp[0]} - 40000*$m;
         $end =$hash{$temp[0]} -40000*$m+40000;
     }
       
  print Out "$temp[0]\t$start\t$end\n"; 
  
 }


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

open BED,"$bed.flank.$i\_K.bed" or die "$!";
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
}

foreach my $k1 (keys %conserved) {
print Out1 "$k1\t$conserved{$k1}\t";
}

foreach my $k2 (keys %missing) {
print Out1 "$k2\t$missing{$k2}\t";
}

foreach my $k3 (keys %SNP) {
print Out1 "$k3\t$SNP{$k3}\t";
}

foreach my $k4 (keys %DEL) {
print Out1 "$k4\t$DEL{$k4}\n";
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

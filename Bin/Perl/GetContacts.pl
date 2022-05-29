#!/usr/bin/perl
use strict;
use warnings;

my %interactions;
my $n=0;
open IN, "$ARGV[0]" or die "$!";
open OUT, ">$ARGV[0].out" or die "$!";
while (<IN>) {
print "Haha $n\t$_\n";
chomp;
  open TSV, "$_" or die "$!";
  while (<TSV>) {
   chomp;
   
   my @temp = split (/\t/,$_);
   my $cor1=$temp[1]/100000;
   my $cor2=$temp[4]/100000;
   $temp[0]=~s/Ca_59//;
   $temp[3]=~s/Ca_59//;
   my $contact = join ("_",($temp[0],$cor1,$temp[3],$cor2));
   if(exists $interactions{$contact}) {
   $interactions{$contact}[$n]=$temp[6];
   } else {
   my @array = (0,0,0,0,0,0,0,0);
   $interactions{$contact}=\@array;
   $interactions{$contact}[$n]=$temp[6];
   }
  }
$n++;
}

foreach my $key (keys %interactions) {
  print OUT "$key\t";
  print OUT join("\t",@{$interactions{$key}}),"\n";
}

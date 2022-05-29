#!/usr/bin/perl
use strict;
use warnings;

my %size;
open In, "$ARGV[0]" or die "$!";
while (<In>) {
chomp;
my @unit = split (/\t/,$_);
$size{$unit[0]} = $unit[1];
}

my $chr = $ARGV[2];
my $window = $ARGV[3];
my %hash;
open TSV, "$ARGV[1]" or die "$!";
while (<TSV>) {
chomp;
my @temp = split (/\t/,$_);
my $key = join ("_",($temp[1],$temp[2],$temp[4],$temp[5]));
$hash{$key} =$temp[6];
}

my $i;
my $j;
my $bins = int($size{$chr}/$window);
#print "bins: $bins\n";
open Out, ">$ARGV[1].TopDpm.matrix" or die "$!";

for ($i=0;$i<$bins;$i++) {
   my $beg = $i*$window;
   my $end = ($i+1)*$window;
   print Out "$chr\t$beg\t$end\t";
   for ($j=0;$j<$bins;$j++) {
   my $beg1 = $j*$window;
   my $end1 = ($j+1)*$window;
   my $ele = join ("_",($beg,$end,$beg1,$end1));
   if ($j!=$bins-1) {
      if (exists $hash{$ele}) {
      print Out "$hash{$ele}\t";
      } else {
      print Out "0\t";
      }
   } else {
      if (exists $hash{$ele}) {
      print Out "$hash{$ele}\n";
      } else {
      print Out "0\n";
      } 
   }   
  }
}

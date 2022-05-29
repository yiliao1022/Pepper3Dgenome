#!/usr/bin/perl
use strict;
use warnings;

open In, "$ARGV[0]" or die "$!";
my @cov;
while (<In>) {
chomp;
push (@cov,$_);
}

open In1, "$ARGV[1]" or die "$!";
my $j=0;
while (<In1>) {
chomp;
my @tmp;
my $k;
for ($k=$j*20;$k<$j*20+20;$k++) {
push (@tmp, $cov[$k]);
}

my $refavg = &Average(@tmp);

my @unit = split (/\t/,$_);
my $m;
for ($m=0;$m<19;$m++) {
my $rate = int ($unit[$m]/$$refavg[$m]);
print "$rate\t";
}
my $ratio = int ($unit[19]/$$refavg[19]);
print "$ratio\n";

$j++;
}


sub Average {
   # get total number of arguments passed.
   my @list = @_;
   my $n = scalar(@list);
   my $sum = 0;

   foreach my $item (@_) {
      $sum += $item;
   }
  my $average = $sum / $n;

   my @avglist;
   my $i;
   for ($i=0;$i<$n;$i++) {
   my $new = $list[$i]/$average;
   push (@avglist,$new);
   }
   my $listref = \@avglist;
   return $listref;
}

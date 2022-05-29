#!/usr/bin/perl 
use strict;
use warnings;

open In, "$ARGV[0]" or die "$!";

my %TwoD;
my $bin = 10000;

while (<In>) {
chomp;
my @unit = split (/\t/,$_);
#$sizes{$unit[0]} = $unit[1];
my $num = int ($unit[1]/$bin);
my $i;
for ($i=0;$i<=$num;$i++) {
my $key = join ("\.",($unit[0],$i));
$TwoD{$key}=0;
 }
}



open Contact, "$ARGV[1]" or die "$!";
open OUT, ">$ARGV[1].fre.out" or die "$!";

while (<Contact>) {
chomp;
my @temp = split (/\s/,$_);
my $n1 = int ($temp[2]/$bin);
my $n2 = int ($temp[6]/$bin);
my $key1 = join ("\.",($temp[1],$n1));
my $key2 = join ("\.",($temp[5],$n2));
$TwoD{$key1}++;
$TwoD{$key2}++;
}

foreach my $window (sort keys %TwoD) {
my @tmp = split (/\./,$window);
my $beg = $tmp[1]*$bin;
my $end = ($tmp[1]+1)*$bin;
print OUT "$tmp[0]\t$beg\t$end\t$TwoD{$window}\n";
}

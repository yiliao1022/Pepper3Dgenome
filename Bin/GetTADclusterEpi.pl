#!/usr/bin/perl
use strict;
use warnings;

my %hash;

open In, "$ARGV[0]" or die "$!";
while (<In>) {
chomp;
next if ($_=~/AChr/);
my @temp = split (/\t/,$_);
$temp[0]=~s/Ca_59//;
my $start = $temp[1];
my $end = $temp[2];
my $key = join ("_",($temp[0],$start,$end));
my @value = ($temp[13],$temp[14],$temp[15],$temp[16],$temp[17],$temp[18],$temp[19],$temp[20]);
$hash{$key}=\@value;
}

open TAD, "$ARGV[1]" or die "$!";
while (<TAD>) {
chomp;
my @rank=(0,0,0,0,0,0,0,0);
my @unit = split (/\t/,$_);
$unit[0]=~s/Ca_59//;
my $n = ($unit[2]-$unit[1])/10000;
my $i;

for ($i=0;$i<$n;$i++) {
my $beg = $unit[1]+$i*10000;
my $end1 = $unit[1]+($i+1)*10000;
my $key1 = join ("_",($unit[0],$beg,$end1));


$rank[0] += ${$hash{$key1}}[0];
$rank[1] += ${$hash{$key1}}[1];
$rank[2] += ${$hash{$key1}}[2];
$rank[3] += ${$hash{$key1}}[3];
$rank[4] += ${$hash{$key1}}[4];
$rank[5] += ${$hash{$key1}}[5];
$rank[6] += ${$hash{$key1}}[6];
$rank[7] += ${$hash{$key1}}[7];

=pod  
   if ($hash{$key1}==0.125) {
       $rank[0]++;
   } elsif ($hash{$key1}==0.25) {
       $rank[1]++;
   } elsif ($hash{$key1}==0.375) {
       $rank[2]++;
   } elsif ($hash{$key1}==0.5) {
       $rank[3]++;
   } elsif ($hash{$key1}==0.625) {
       $rank[4]++;
   } elsif ($hash{$key1}==0.75) {
       $rank[5]++;
   } elsif ($hash{$key1}==0.875) {
       $rank[6]++;
   } elsif ($hash{$key1}==1) {
       $rank[7]++;
   }
=cut 

}

my $tad = "$unit[0]\_$unit[1]\_$unit[2]";
my $m1 = $rank[0]/$n;
my $m2 = $rank[1]/$n;
my $m3 = $rank[2]/$n;
my $m4 = $rank[3]/$n;
my $m5 = $rank[4]/$n;
my $m6 = $rank[5]/$n;
my $m7 = $rank[6]/$n;
my $m8 = $rank[7]/$n;
print "$tad\t$m1\t$m2\t$m3\t$m4\t$m5\t$m6\t$m7\t$m8\n"
}

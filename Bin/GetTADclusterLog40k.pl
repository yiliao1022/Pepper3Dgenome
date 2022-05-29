#!/usr/bin/perl
use strict;
use warnings;

my %hash;

open In, "$ARGV[0]" or die "$!";
while (<In>) {
chomp;
next if ($_=~/AChr/);
my @temp = split (/\s/,$_);
$temp[12]=~s/Ca_59//;
my $start = $temp[13];
my $end = $temp[14];
my $key = join ("_",($temp[12],$start,$end));
my @value = ($temp[1],$temp[2],$temp[3],$temp[4],$temp[5],$temp[6],$temp[7],$temp[8],$temp[9],$temp[10],$temp[11]);
$hash{$key}=\@value;
}

open TAD, "$ARGV[1]" or die "$!";
while (<TAD>) {
chomp;
my @rank=(0,0,0,0,0,0,0,0,0,0,0);
my @unit = split (/\t/,$_);
$unit[0]=~s/Ca_59//;
my $n = ($unit[2]-$unit[1])/40000;
my $i;
my @NA = (0,0,0,0,0,0,0,0,0,0,0);
for ($i=0;$i<$n;$i++) {
my $beg = $unit[1]+$i*40000;
my $end1 = $unit[1]+($i+1)*40000;
my $key1 = join ("_",($unit[0],$beg,$end1));
#my @NA = (0,0,0,0,0,0,0,0,0,0,0);

if (${$hash{$key1}}[0]=~/NA/) {
$NA[0]++;
} else{
$rank[0] += ${$hash{$key1}}[0];
}

if (${$hash{$key1}}[1]=~/NA/) {
$NA[1]++;
} else{
$rank[1] += ${$hash{$key1}}[1];
}

if (${$hash{$key1}}[2]=~/NA/) {
$NA[2]++;
} else{
$rank[2] += ${$hash{$key1}}[2];
}

if (${$hash{$key1}}[3]=~/NA/) {
$NA[3]++;
} else{
$rank[3] += ${$hash{$key1}}[3];
}

if (${$hash{$key1}}[4]=~/NA/) {
$NA[4]++;
} else{
$rank[4] += ${$hash{$key1}}[4];
}

if (${$hash{$key1}}[5]=~/NA/) {
$NA[5]++;
} else{
$rank[5] += ${$hash{$key1}}[5];
}

if (${$hash{$key1}}[6]=~/NA/) {
$NA[6]++;
} else{
$rank[6] += ${$hash{$key1}}[6];
}


if (${$hash{$key1}}[7]=~/NA/) {
$NA[7]++;
} else{
$rank[7] += ${$hash{$key1}}[7];
}

if (${$hash{$key1}}[8]=~/NA/) {
$NA[8]++;
} else{
$rank[8] += ${$hash{$key1}}[8];
}

if (${$hash{$key1}}[9]=~/NA/) {
$NA[9]++;
} else{
$rank[9] += ${$hash{$key1}}[9];
}

if (${$hash{$key1}}[10]=~/NA/) {
$NA[10]++;
} else{
$rank[10] += ${$hash{$key1}}[10];
}

#$rank[1] += ${$hash{$key1}}[1];
#$rank[2] += ${$hash{$key1}}[2];
#$rank[3] += ${$hash{$key1}}[3];
#$rank[4] += ${$hash{$key1}}[4];
#$rank[5] += ${$hash{$key1}}[5];
#$rank[6] += ${$hash{$key1}}[6];
#$rank[7] += ${$hash{$key1}}[7];
#$rank[8] += ${$hash{$key1}}[8];
#$rank[9] += ${$hash{$key1}}[9];
#$rank[10] += ${$hash{$key1}}[10];




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
my $m1 = $rank[0]/($n-$NA[0]);
my $m2 = $rank[1]/($n-$NA[1]);
my $m3 = $rank[2]/($n-$NA[2]);
my $m4 = $rank[3]/($n-$NA[3]);
my $m5 = $rank[4]/($n-$NA[4]);
my $m6 = $rank[5]/($n-$NA[5]);
my $m7 = $rank[6]/($n-$NA[6]);
my $m8 = $rank[7]/($n-$NA[7]);
my $m9 = $rank[8]/($n-$NA[8]);
my $m10 = $rank[9]/($n-$NA[9]);
my $m11 = $rank[10]/($n-$NA[10]);
print "$tad\t$m1\t$m2\t$m3\t$m4\t$m5\t$m6\t$m7\t$m8\t$m9\t$m10\t$m11\n"
}

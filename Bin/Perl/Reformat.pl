#!/bin/perl
use strict;
use warnings;
open In,"$ARGV[0]" or die "$!";
open Out,">$ARGV[0].out" or die "$!";
while (<In>) {
chomp;
my @temp=split(/\t/,$_);
my @unit=split(/_/,$temp[0]);
my $chr=join("_",($unit[0],$unit[1]));
my $beg = $unit[2]*40000;
my $end = $unit[2]*40000+40000;
my $GR=0.1+($temp[1]+$temp[2]+$temp[3])/3;
my $HL=0.1+($temp[4]+$temp[5]+$temp[6])/3;
my $TZ=0.1+($temp[7]+$temp[8]+$temp[9])/3;
my $YP=0.1+($temp[10]+$temp[11]+$temp[12])/3;
my $RO=0.1+($temp[13]+$temp[14]+$temp[15])/3;
my $GRvsHL;
my $GRvsTZ;
my $GRvsYP;
my $HLvsTZ;
my $HLvsYP;
my $TZvsYP;

if ($GR<1 && $HL<1 && $TZ<1 && $YP<1 && $RO<1) {
$GRvsHL=0.0;
$GRvsTZ=0.0;
$GRvsYP=0.0;
$HLvsTZ=0.0;
$HLvsYP=0.0;
$TZvsYP=0.0;
} else {
$GRvsHL=&log2($GR/$HL);
$GRvsTZ=&log2($GR/$TZ);
$GRvsYP=&log2($GR/$YP);
$HLvsTZ=&log2($HL/$TZ);
$HLvsYP=&log2($HL/$YP);
$TZvsYP=&log2($TZ/$HL);
}
print Out "$chr\t$beg\t$end\t$GRvsHL\t$GRvsTZ\t$GRvsYP\t$HLvsTZ\t$HLvsYP\t$TZvsYP\n";
}

sub log2 
{
    my $n = shift;     
    # using pre-defined log function
    return log($n) / log(2);
}

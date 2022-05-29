#!/bin/perl
use strict;
use warnings;
open In,"$ARGV[0]" or die "$!";
open Out,">$ARGV[0].quandif.out" or die "$!";
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
$GRvsHL=-0.1;
$GRvsTZ=-0.1;
$GRvsYP=-0.1;
$HLvsTZ=-0.1;
$HLvsYP=-0.1;
$TZvsYP=-0.1;
} else {
my $GRvsHLmax=&max($GR,$HL);
my $GRvsTZmax=&max($GR,$TZ);
my $GRvsYPmax=&max($GR,$YP);
my $HLvsTZmax=&max($HL,$TZ);
my $HLvsYPmax=&max($HL,$YP);
my $TZvsYPmax=&max($TZ,$HL);
my $GRvsHLabs=abs($GR-$HL);
my $GRvsTZabs=abs($GR-$TZ);
my $GRvsYPabs=abs($GR-$YP);
my $HLvsTZabs=abs($HL-$TZ);
my $HLvsYPabs=abs($HL-$YP);
my $TZvsYPabs=abs($TZ-$HL);


$GRvsHL=$GRvsHLabs/$GRvsHLmax;
$GRvsTZ=$GRvsTZabs/$GRvsTZmax;
$GRvsYP=$GRvsYPabs/$GRvsYPmax;
$HLvsTZ=$HLvsTZabs/$HLvsTZmax;
$HLvsYP=$HLvsYPabs/$HLvsYPmax;
$TZvsYP=$TZvsYPabs/$TZvsYPmax;
}
print Out "$chr\t$beg\t$end\t$GRvsHL\t$GRvsTZ\t$GRvsYP\t$HLvsTZ\t$HLvsYP\t$TZvsYP\n";
}

sub log2 
{
    my $n = shift;     
    # using pre-defined log function
    return log($n) / log(2);
}

sub max {
my ($a1,$a2) = @_;
my $max;
  if ($a1>$a2) {
      $max=$a1;
  } else {
      $max=$a2;
  }
}

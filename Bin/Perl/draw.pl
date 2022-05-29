#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use SVG;

my ($size,$TEdensity,$svg2xxx,$Help);

GetOptions(
  "size:s"=>\$size,
  "density:s"=>\$TEdensity,
  "svg2xxx:s"=>\$svg2xxx,
  "help"     =>\$Help
);

##################################################################################################
if ($Help){
print <<"END.";
  Usage: perl $0 -size genome.size -density LTR.density -svg2xxx /path/to/svg2xxx
  -size           chromosome sizes                                                   [REQUIRED]
  -density        TE density                                                         [REQUIRED]
  -svg2xxx        Path to svg2xxx tools                                              [REQUIRED]
  -help           help
END.
exit;
}
###################################################################################################

my $svg= SVG -> new (width=>800, height =>400);
my $rate = 300/300000000;

open SIZE, "$size" or die "$!";


while (<SIZE>) {
chomp;
my @temp = split (/\t/,$_);
$hash{$temp[0]}=$temp[1];
}

####Chromosome
$svg->rectangle (x=>50,y=>50,width=>$hash{"Chr01"}*$rate,height=>4,style => { 'fill'=> "yellow", 'stroke'=> "yellow"});
$svg->rectangle (x=>50,y=>100,width=>$hash{"Chr02"}*$rate,height=>4,style => { 'fill'=> "yellow", 'stroke'=> "yellow"});
$svg->rectangle (x=>50,y=>150,width=>$hash{"Chr03"}*$rate,height=>4,style => { 'fill'=> "yellow", 'stroke'=> "yellow"});
$svg->rectangle (x=>50,y=>200,width=>$hash{"Chr04"}*$rate,height=>4,style => { 'fill'=> "yellow", 'stroke'=> "yellow"});
$svg->rectangle (x=>50,y=>250,width=>$hash{"Chr05"}*$rate,height=>4,style => { 'fill'=> "yellow", 'stroke'=> "yellow"});
$svg->rectangle (x=>50,y=>300,width=>$hash{"Chr06"}*$rate,height=>4,style => { 'fill'=> "yellow", 'stroke'=> "yellow"});
$svg->rectangle (x=>450,y=>50,width=>$hash{"Chr07"}*$rate,height=>4,style => { 'fill'=> "yellow", 'stroke'=> "yellow"});
$svg->rectangle (x=>450,y=>100,width=>$hash{"Chr08"}*$rate,height=>4,style => { 'fill'=> "yellow", 'stroke'=> "yellow"});
$svg->rectangle (x=>450,y=>150,width=>$hash{"Chr09"}*$rate,height=>4,style => { 'fill'=> "yellow", 'stroke'=> "yellow"});
$svg->rectangle (x=>450,y=>200,width=>$hash{"Chr10"}*$rate,height=>4,style => { 'fill'=> "yellow", 'stroke'=> "yellow"});
$svg->rectangle (x=>450,y=>250,width=>$hash{"Chr11"}*$rate,height=>4,style => { 'fill'=> "yellow", 'stroke'=> "yellow"});
$svg->rectangle (x=>450,y=>300,width=>$hash{"Chr12"}*$rate,height=>4,style => { 'fill'=> "yellow", 'stroke'=> "yellow"});

#########
open DENSITY, "$TEdensity" or die "$!";
while (<DENSITY>) {
my @unit=split(/\t/,$_);

if ($unit[0]=~/Chr01/) {
my $x = 50+$unit[1]*$rate,y=>50-50*$unit[3]/250,width=>5000000*$rate,$height=50*$unit[3]/250;style => { 'fill'=> "red", 'stroke'=> "red"});
my $x = 50+$unit[1]*$rate,y=>50-50*$unit[3]/250-50*$unit[4]/250,width=>5000000*$rate,$height=50*$unit[4]/250;style => { 'fill'=> "orange", 'stroke'=> "orange"});
my $x = 50+$unit[1]*$rate,y=>50-50*$unit[3]/250-50*$unit[4]/250-50*$unit[5]/250,width=>5000000*$rate,$height=50*$unit[5]/250;style => { 'fill'=> "green", 'stroke'=> "green"});
my $x = 50+$unit[1]*$rate,y=>50-50*$unit[3]/250-50*$unit[4]/250-50*$unit[5]/250-50*$unit[6]/250,width=>5000000*$rate,$height=50*$unit[6]/250;style => { 'fill'=> "blue", 'stroke'=> "blue"});
}

}

open OUT, ">temp.svg" or die "Can not open my file";
print OUT $svg->xmlify;
close OUT;

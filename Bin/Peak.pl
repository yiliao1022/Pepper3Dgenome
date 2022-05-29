#!/usr/bin/perl
use strict;
use warnings;

my $i;
my $j;

for ($i=1;$i<=2;$i++) {
system "perl 10_MAF2Var.pl -fasize CA59.sizes -input ./sing/Chr0$i/CA59.CM334.fa -glst glst.bed -bed All.5.40k.boundaries.bed.merge.bed.bed -chr Chr0$i";
  for ($j=-100;$j<=100;$j++) {
  system "cat All.5.40k.boundaries.bed.merge.bed.bed.flank.$j\_K.bed.bed >> Chr0$i.bed";
  }
system "rm *.flank.*";
}

=pod
for ($i=10;$i<=12;$i++) {
system "perl 10_MAF2Var.pl -fasize CA59.sizes -input ./sing/Chr$i/CA59.CM334.fa -glst glst.bed -bed All.5.40k.boundaries.bed.merge.bed.bed -chr Chr$i";
  for ($j=-100;$j<=100;$j++) {
  system "cat All.5.40k.boundaries.bed.merge.bed.bed.flank.$j\_K.bed.bed >> Chr$i.bed";
  }
system "rm *.flank.*";
}
=cut

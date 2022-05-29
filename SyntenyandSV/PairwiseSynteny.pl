use strict;
use warnings;

my $chr;

open In, "$ARGV[0]" or die "$!";
open Out, ">$ARGV[0].bed" or die "$!";

while (<In>) {
  if ($_=~/net/) {
  my @tmp = split (/\s/,$_);
  $chr = $tmp[1];
  }

  if ($_=~/fill/) {
   my @tmp = split (/\s+/,$_);
   my $start = $tmp[2];
   my $end = $tmp[2] + $tmp[3];
   my $s1 = $start - 20000;
   my $e1 = $start + 20000;
   my $s2 = $end - 20000;
   my $e2 = $end + 20000;  
   print Out "$chr\t$s1\t$e1\n";
   print Out "$chr\t$s2\t$e2\n";
  }
 
# print Out "$chr\t$start\t$end\n"; 
}

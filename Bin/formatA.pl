use strict;
use warnings;
open In, "$ARGV[0]" or die "$!";
open Out, ">$ARGV[0].bed" or die "$!";
while (<In>){
chomp;
my @temp = split (/\t/,$_);
if ($_=~/LJGR/) {
  print Out "$_\n";
} elsif ($temp[1]<10 && $temp[2]<10 && $temp[3]<10 && $temp[4]<10 && $temp[5]<10 && $temp[6]<10 && $temp[7]<10 && $temp[8]<10 && $temp[9]<10 && $temp[10]<10 && $temp[11]<10 && $temp[12]<10 )  {

} else {
print Out "$_\n";
}


}

use strict;
use warnings;

open In, "$ARGV[0]" or die "$!";
my %num;
my %counts;
while (<In>) {
my @temp = split (/\t/,$_);
my $key = $temp[4]-$temp[1];

if (exists $num{$key}){
$num{$key}++;
} else {
$num{$key}=1;
}

if (exists $num{$key}){
$counts{$key}+=$temp[6];
} else {
$counts{$key}=$temp[6];
}
}

foreach my $n (keys %num) {
my $ave = int ($counts{$n}/$num{$n});
print "$n\t$ave\n"

}

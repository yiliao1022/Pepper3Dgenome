use strict;
use warnings;

my %hash;
my $chr;

open In, "$ARGV[0]" or die "$!";
while (<In>) {
chomp;
my @temp = split (/\t/,$_);
$chr = $temp[0];
my $key1 = ($temp[1]+$temp[2])/2;
my $key2 = ($temp[4]+$temp[5])/2;

if (exists $hash{$key1}) {
$hash{$key1}+=$temp[6];
} else {
$hash{$key1} = $temp[6];
}

if (exists $hash{$key2}) {
$hash{$key2}+=$temp[6];
} else {
$hash{$key2} = $temp[6];
}


}

foreach my $mid (sort {$a <=> $b} keys %hash) {
print "$chr\t0\t$mid\t$hash{$mid}\t1\n";

}

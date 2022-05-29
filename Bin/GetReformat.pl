use strict;
use warnings;

my %hash;
open GEXP, "$ARGV[0]" or die "NO gene expression matric!";

while (<GEXP>) {
chomp;
my @temp = split (/\t/,$_);
my $bin = join ("-",($temp[0],$temp[1],$temp[2]));
shift @temp;
shift @temp;
shift @temp;
my $ref = \@temp;
$hash{$bin} = $ref;
}



open In, "$ARGV[1]" or die "No group file!";
my $n = 1;
while (<In>) {
chomp;
my @unit = split (/\t/,$_);
my %hash1;
foreach my $bin (@unit)  {
my $start = $bin - 20000;
my $end = $bin + 20000;
my $value = join("-",("Ca_59Chr08",$start,$end));
$hash1{$bin} = $value;
}

print "### Group $n:\n";
foreach my $mid (sort {$a <=> $b} keys %hash1) {
print "$hash1{$mid}\t@{$hash{$hash1{$mid}}}\n";
}
$n++;

}

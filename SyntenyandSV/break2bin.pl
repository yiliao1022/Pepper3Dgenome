#!/usr/bin/perl
use strict;
use warnings;

open In, "$ARGV[0]" or die "$!";
my @cov;
while (<In>) {
chomp;
push (@cov,$_);
}

open In1, "$ARGV[1]" or die "$!";
my $j=0;
while (<In1>) {
chomp;
my @tmp;
my $k;
for ($k=$j*20;$k<$j*20+20;$k++) {
push (@tmp, $cov[$k]);
}

my $refavg = &Average(@tmp);

my @unit = split (/\t/,$_);
my $m;
for ($m=0;$m<19;$m++) {
my $rate = int ($unit[$m]/$$refavg[$m]);
print "$rate\t";
}
my $ratio = int ($unit[19]/$$refavg[19]);
print "$ratio\n";

$j++;
}


sub Average {
   # get total number of arguments passed.
   my @list = @_;
   my $n = scalar(@list);
   my $sum = 0;

   foreach my $item (@_) {
      $sum += $item;
   }
  my $average = $sum / $n;

   my @avglist;
   my $i;
   for ($i=0;$i<$n;$i++) {
   my $new = $list[$i]/$average;
   push (@avglist,$new);
   }
   my $listref = \@avglist;
   return $listref;
}

(base) liaoy12@thoth:~/Pepper/Peper_2021/1_Synteny/align/NEW/SL4/NEW/scale$ more break2bin.pl 
#!/usr/bin/perl
use strict;
use warnings;

my %hash;
open In, "$ARGV[0]" or die "$!";
while (<In>) {
chomp;
my @temp = split(/\t/,$_);
my $coordinate = join ("_",($temp[1],$temp[2]));
push (@{$hash{$temp[0]}},$coordinate);
}

open SV, "$ARGV[1]" or die "$!";
my @breaks;
while (<SV>) {
chomp;
my @unit = split (/\t/,$_);
my $ele = join ("_",($unit[0],$unit[1],$unit[2]));
push (@breaks,$ele);
}

my @counts=(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
foreach my $key (%hash) {
   foreach my $coordinate (@{$hash{$key}}) {
       my @bins;
       my @tmp = split (/_/,$coordinate);
       my $len = $tmp[1]-$tmp[0];
       my $start = int ($tmp[0]-($len/2));
       my $end = int ($tmp[1]+($len/2));
       my $bin = int ($len/10);
       my $i;

       for ($i=0;$i<20;$i++) {
            $bins[$i]=0;
            my $s = $start+$i*$bin;
            my $e = $start+($i+1)*$bin;
             foreach my $ind (@breaks) {
                my @peak_ind = split (/_/,$ind);
                my $center = int (($peak_ind[1]+$peak_ind[2])/2);
                if ($peak_ind[0] eq $key) {
                if ($s<=$center and $center<=$e) {          
                            $bins[$i]++;
                         }       
                    } 
               }
          }

            #print "$key\t$coordinate\t";
            #print @bins;    
                
         my $j;
         for ($j=0;$j<20;$j++) {
             $counts[$j]=$counts[$j]+$bins[$j];
         }
    }            
  }

foreach my $num (@counts) {
  print "$num\t";
}
print "\n";

##!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);

####################################################################################
# ToAxt.pl
#
# Authors: Yi Liao (05/28/2020) Shujun Ou (05/31/2020)
#
# Copyright (C) Not for commercial use
#
# Prerequisite : Kent's utilities; minimap2; paftools.js, please make sure all these utilities are on your $PATH before run it.
#
# Convert genome alignment file to axt format (LAST, LASTZ, MUMmer,minimap2 et al...) and call SVs using paftools.js
#
# Usage:  perl ToAxt.pl -aligner minimap2 -wkfolder . -tname A -qname B
#
# Options: -aligner  [last|lastz|MUMmer|minimap2] aligner for whole genome pairwise alignment      [REQUIRED]
#          -wkfolder [FILE] where the alignment files put                                          [REQUIRED]
#          -tname    [FILE] the name of target genome assembly                                     [REQUIRED]
#          -qname    [FILE] the name of query genome assembly                                      [REQUIRED]
#          -help     print this information
####################################################################################

my ($aligner,$wk,$tname,$qname,$MUM320,$kentUtil,$TBA,$paftools,$k8,$Help);

# get input parameters
GetOptions( 'tname=s'=>\$tname,
	    'qname=s'=>\$qname,
            'aligner=s'=>\$aligner,
            'wkfolder=s'=>\$wk,
	    'mum320=s'=>\$MUM320,
            'kentUtil=s'=>\$kentUtil,
            'tba=s'=>\$TBA,
	    'k8=s'=>\$k8,
            'paftools=s'=>\$paftools,
            'help'=>\$Help
           );

# program path. Will use the ones in $PATH if unspecified.
$k8 ||='';
$paftools ||='';
$TBA ||='';
$kentUtil ||='';
$MUM320 ||='';
$wk ||='.';

if ($Help){
print <<"END.";
  Usage:  perl ToAxt.pl -aligner minimap2 -wkfolder /path/to/align -kentUtil /path/to/kentUtil -tba /path/to/TBA -k8 /path/to/k8 -paftools /path/to/paftools -tname A -qname B

  Options:  -aligner  [last|lastz|MUMmer|minimap2] aligner for whole genome pairwise alignment      [REQUIRED]
            -wkfolder [path] where the alignment files put                                          [Optional]
	    -kentUtil [path] where the kent utilities programs put			            [Optional]
            -tba      [path] where the TBA program put                                              [Optional]
	    -mum320   [path] where the MUMmer3.20 put                                               [Optional]
            -paftools [path] where the paftools program put                                         [Optional]
            -tname    [FILE] the name of target genome assembly                                     [REQUIRED]
            -qname    [FILE] the name of query genome assembly                                      [REQUIRED]
            -help     print this information                                                        [Optional]
END.
exit;
}

# make output directories
`mkdir $wk/Target` unless -e "$wk/Target" && -d "$wk/Target";
`mkdir $wk/Query` unless -e "$wk/Query" && -d "$wk/Query";

if ($aligner eq "minimap2") {

opendir (DIR,$wk) || die "Error in opening dir $wk\n";
while ((my $filename = readdir (DIR))) {
if ($filename =~/(.*)vs(.*).long.paf/) {
`sort -k6,6 -k8,8n $filename > $1.$2.sorted.paf;
${k8}k8 ${paftools}paftools.js call $1.$2.sorted.paf > $1.$2.sorted.var.txt;
${k8}k8 ${paftools}paftools.js view -f maf $1.$2.sorted.paf | sed 's/^a /a score=/g; s/_/./g' > $1.$2.sorted.maf;
${TBA}maf_order $1.$2.sorted.maf $2 $1 > $2.$1.sorted.maf;
${kentUtil}mafToAxt $1.$2.sorted.maf $1 $2 $1.$2.axt;
${kentUtil}mafToAxt $2.$1.sorted.maf $2 $1 $2.$1.axt;
mv $1.$2.axt $wk/Target;
mv $2.$1.axt $wk/Query;
rm $1.$2.sorted.maf $2.$1.sorted.maf;
rm $1.$2.sorted.paf;
`
}

}

=pod
`for i in \`ls $wk/*.paf|grep -v sorted\`; do
i=\$(echo \$i|perl -nle 's/.long.paf//; print \$_');
t=\$(echo \$i|perl -nle 's/vs(.*)//; print \$_');
q=\$(echo \$i|perl -nle 's/(.*)vs//; print \$_');
sort -k6,6 -k8,8n \$i.long.paf > \$i.sorted.paf;
${k8}k8 ${paftools}paftools.js call \$i.sorted.paf > \$i.sorted.var.txt;
${k8}k8 ${paftools}paftools.js view -f maf \$i.sorted.paf | sed 's/^a /a score=/g; s/_/./g' > \$i.sorted.maf;
${TBA}maf_order \$i.sorted.maf \$q \$t > \$i.\$q.\$t.maf;
${kentUtil}mafToAxt \$i.sorted.maf \$t \$q  \$i.\$t.\$q.axt;
${kentUtil}mafToAxt \$i.\$q.\$t.maf \$q \$t \$i.\$q.\$t.axt;
mv \$i.\$t.\$q.axt $wk/Target;
mv \$i.\$q.\$t.axt $wk/Query;
rm \$i.\$q.\$t.maf \$i.sorted.maf;
rm \$i.sorted.paf;
done`;
=cut


} elsif ($aligner eq "mummer") {
`for i in \`ls $wk/*.delta\`; do
 i=\$(echo \$i|perl -nle 's/.delta//; print \$_');
 ${MUM320}delta2maf -r \$i.delta > \$i.sorted.maf;
 ${TBA}maf_order \$i.sorted.maf $qname $tname > \$i.$qname.$tname.maf;
 ${kentUtil}mafToAxt \$i.sorted.maf $tname $qname \$i.$tname.$qname.axt;
 ${kentUtil}mafToAxt \$i.$qname.$tname.maf $qname $tname \$i.$qname.$tname.axt;
 mv \$i.$tname.$qname.axt $wk/Target_$tname;
 mv \$i.$qname.$tname.axt $wk/Target_$qname;
 rm \$i.$qname.$tname.maf \$i.sorted.maf;
 done`;
} elsif ($aligner eq "last") {
`for i in \`ls $wk/*.maf|grep -v sorted\`; do
i=\$(echo \$i|perl -nle 's/.maf//; print \$_');
sed -i '1i \#\#maf version=1' \$i.maf; 
${TBA}maf_sort \$i.maf $tname > \$i.sorted.maf;
${TBA}maf_order \$i.sorted.maf $qname $tname > \$i.$qname.$tname.maf;
${kentUtil}mafToAxt \$i.maf $tname $qname \$i.$tname.$qname.axt;
${kentUtil}mafToAxt \$i.$qname.$tname.maf $qname $tname \$i.$qname.$tname.axt;
mv \$i.$tname.$qname.axt $wk/Target_$tname;
mv \$i.$qname.$tname.axt $wk/Target_$qname;
rm \$i.$qname.$tname.maf \$i.sorted.maf;
done`;
} elsif ($aligner eq "lastz") {
`for i in \`ls $wk/*.maf|grep -v sorted\`; do
i=\$(echo \$i|perl -nle 's/.maf//; print \$_');
${TBA}maf_sort \$i.maf $tname > \$i.sorted.maf;
${TBA}maf_order \$i.sorted.maf $qname $tname > \$i.$qname.$tname.maf;
${kentUtil}mafToAxt \$i.maf $tname $qname \$i.$tname.$qname.axt;
${kentUtil}mafToAxt \$i.$qname.$tname.maf $qname $tname \$i.$qname.$tname.axt;
mv \$i.$tname.$qname.axt $wk/Target_$tname;
mv \$i.$qname.$tname.axt $wk/Target_$qname;
rm \$i.$qname.$tname.maf \$i.sorted.maf;
done`;
} else {
print "Please define the aligner you used! Or You may use an unsurpported or unknown aligner!"
}

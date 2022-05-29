#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use lib "$FindBin::Bin";
use lib "$FindBin::Bin/./pub";
use pub::Tree::DAG_Node;
use lib::SV qw(SynNet ExtractSeq ParseFill);

####################################################################################
# 5_PairwiseASSYSV.pl
#
# Authors: Yi Liao (06/01/2020) <yiliao1022@gmail.com>
#
# Copyright (C) Not for commercial use
#
# Preliminary SVs calling based on ChainNetSynnet output file
#
# Prerequisite : Kent's utilities; Perl module Tree::DAG_Node
#
#
####################################################################################

#### Global Variablelist
my ($inSynnet,$outdir,$refseq,$Queseq,$help);
my ($invlen,$cnvlen,$dellen,$trllen,$toplen,$secondlen,$ali_rate);
#my ($OutIndel,$OutInversion,$OutComplex,$OutCnv,$OutTranslocation);

GetOptions(
  "inSynnet:s"=>\$inSynnet,
	"outdir:s"=>\$outdir,
	"refseq:s"=>\$refseq,
  "qseq:s"=>\$Queseq,
	"invlen:i"=>\$invlen,
	"cnvlen:i"=>\$cnvlen,
	"dellen:i"=>\$dellen,
  "trllen:i"=>\$trllen,
  "rate:f"=>\$ali_rate,
	"filltoplen:i"=>\$toplen,
	"fill2len:i"=>\$secondlen,
	"help"=>\$help
);

#### Printing help information
if (defined ($help) || !defined($inSynnet) || !defined($refseq) || !defined($Queseq)) {
print <<"END.";
Usage: perl $0 --inSynNet <SynNet file> --refseq <Reference genome sequence> --qseq <Query genome sequence> --outdir rawSV

Options
--inSynNet      <s> : Chain/Net/Synnet output file
--outdir        <s> : dir for output [default: current folder]
--refseq        <s> : reference/target genome sequence
--qseq          <s> : query genome sequence

--invlen        <n> : Minimum size of inversions to report [10,000]
--cnvlen        <n> : Minimum size of CNVs (tandom duplications) to report [500]
--dellen        <n> : Maximum size of deletions to report [20,000]
--trllen        <n> : Minumum size of translocation to report [5,000]
--filltoplen    <n> : Minimum size of the top fills to be considered [20,000]
--fill2len      <n> : Minimum size of the second level fills to be considered [10,000]
--rate          <p> : Minimum proportion of bases in query region that overlaping with the target gap that are need to for local duplications [0.25]
--help              : show this help information
END.
exit;
}

#### set default parameters
$outdir ||= ".";
$invlen ||= 10000; # the minimum inversion size to report
$cnvlen ||= 500; # the minimum CNV size to report
$dellen ||= 20000; # the maximum deletion size to report
$trllen ||= 5000;
$ali_rate ||= 0.25;
$toplen ||= 20000; # the minimum top fill size to report
$secondlen ||= 10000; # the minimum seconde fill size to report
my $PAV=0; # Count the sequence content of Presence and absence variations

#### creating output files
$outdir =~ s/\/$//;
mkdir($outdir) unless(-d $outdir);
my $basename = basename ($inSynnet);
my $IndelFile="$outdir/$basename.indel";
my $InvFile="$outdir/$basename.inv";
my $CxSVFile="$outdir/$basename.cxsv";
my $CNVFile="$outdir/$basename.cnv";
my $TRLFile="$outdir/$basename.trl";
#####
open (DEL,">$IndelFile") || die "Cannot open $IndelFile\n";
open (INV,">$InvFile") || die "Cannot open $InvFile\n";
open (CXS,">$CxSVFile") || die "Cannot open $CxSVFile\n";
open (CNV,">$CNVFile") || die "Cannot open $CNVFile\n";
open (TRL,">$TRLFile") || die "Cannot open $TRLFile\n";

my $HEAD1 = "#Reference\tRef_beg\tRef_end\tSVtypes\tSample\tSample_beg\tSample_end\tSVlength\tStrand\t5_flanking\t3_flanking\tRef_5_flanking_20bp\tRef_3_flanking_20bp\tSample_5_flanking_20bp\tSample_3_flanking_20bp\tRef_SV_sequence\tSample_SV_sequence\n";

#print INV "Reference\tRef_beg\tRef_end\tSVtypes\tSample\tSample_beg\tSample_end\tSVlength\tStrand\t5'-flanking\t3'-flanking\tRef_5'-flanking_20bp\tRef_3'-flanking_20bp\tSample_5'-flanking_20bp\tSample_3'-flanking_20bp\tRef_SV_sequence\tSample_SV_sequence\n";
#open (CXS,">$CxSVFile") || die "Cannot open $CxSVFile\n";
#print CXS "Reference\tRef_beg\tRef_end\tSVtypes\tSample\tSample_beg\tSample_end\tSVlength\tStrand\t5'-flanking\t3'-flanking\tRef_5'-flanking_20bp\tRef_3'-flanking_20bp\tSample_5'-flanking_20bp\tSample_3'-flanking_20bp\tRef_SV_sequence\tSample_SV_sequence\n";
#open (CNV,">$CNVFile") || die "Cannot open $CNVFile\n";
#print CNV "Reference\tRef_beg\tRef_end\tSVtypes\tSample\tSample_beg\tSample_end\tSVlength\tStrand\t5'-flanking\t3'-flanking\tRef_5'-flanking_20bp\tRef_3'-flanking_20bp\tSample_5'-flanking_20bp\tSample_3'-flanking_20bp\tRef_SV_sequence\tSample_SV_sequence\n";
#open (TRL,">$TRLFile") || die "Cannot open $TRLFile\n";

my $HEAD2 = "#Reference\tRef_beg\tRef_end\tSample\tSample_beg\tSample_end\tStrand\tRef_len\tSample_len\tTRLtypes\tLevel\n";

###############################################################
###############################################################
my $root = SynNet($inSynnet);
my $Tseq = ExtractSeq($refseq);
my $Qseq = ExtractSeq($Queseq);

my @daughters = $root->daughters; # Load chromsome name into @Chromosome
foreach my $daughter (@daughters) { # For each chromosome or contig
        my $name1=$daughter->name;
        my @dayghters_playerA=$daughter->daughters; # Load top fill into @dayghters_playerA
        foreach my $playerA_dau (@dayghters_playerA) { # For
                my $mama=$playerA_dau->mother;
                my $mamaname=$mama->name; # Name of chromosomes/contigs
	              my $topfill = $playerA_dau->name; # The top fill line
	              my @unit = split(/_/,$topfill);
                if  ($unit[12]>=$toplen) {
		                 &Parse_Fill_branch($playerA_dau,$Tseq,$Qseq,$mamaname,$invlen,$cnvlen,$dellen,$ali_rate);
                     my @gap1=$playerA_dau->daughters;
		                 foreach my $gap (@gap1) {
                             my  @downfill=$gap->daughters;
	                           foreach my $playerB_dau (@downfill) {
				                             my $secondfill = $playerB_dau->name;
			                               my @secondunit = split (/_/,$secondfill);
			                               if ($secondunit[12]>=$secondlen) {
			                                  &Parse_Fill_branch($playerB_dau,$Tseq,$Qseq,$mamaname,$invlen,$cnvlen,$dellen,$ali_rate);
                                        my @gap2 = $playerB_dau->daughters;
                                        foreach my $gap0 (@gap2) {
                                                my @downfill2=$gap0->daughters;
                                                foreach my $playerC_dau (@downfill2) {
                                                        my $thirdfill = $playerC_dau->name;
                                                        my @thirdunit = split (/_/,$thirdfill);
                                                        if ($thirdunit[12]>=$secondlen) {
                                                           &Parse_Fill_branch($playerC_dau,$Tseq,$Qseq,$mamaname,$invlen,$cnvlen,$dellen,$ali_rate);
                                                        }
                                                }                                                                                          
                                           }
		                                 }
			                       }
	                   }
               } 
        }
}



`cat $IndelFile | sort -k1,1 -k2,2n > $IndelFile.sort`;
`cat $InvFile | sort -k1,1 -k2,2n > $InvFile.sort`;
`cat $CxSVFile | sort -k1,1 -k2,2n > $CxSVFile.sort`;
`cat $CNVFile | sort -k1,1 -k2,2n > $CNVFile.sort`;
`cat $TRLFile | sort -k1,1 -k2,2n > $TRLFile.sort`;

`sed -i '1i $HEAD1' $IndelFile.sort`;
`sed -i '1i $HEAD1' $InvFile.sort`;
`sed -i '1i $HEAD1' $CxSVFile.sort`;
`sed -i '1i $HEAD1' $CNVFile.sort`;
`sed -i '1i $HEAD2' $TRLFile.sort`;

`rm $IndelFile $InvFile $CxSVFile $CNVFile $TRLFile`;

close DEL;
close INV;
close CXS;
close CNV;
close TRL;
#print "The total presence/absence variations account for $PAV bp\n";

###########################
######  Subroutines  ######
###########################

sub Parse_Fill_branch {
 my ($fill_branch,$tseq,$qseq,$mamaname,$INVlen,$CNVlen,$DELlen,$ali_rate) = @_;
 my $fillname = $fill_branch->name; # The top fill line
 my $hash_ref = ParseFill($fillname);
 my @unit = split(/_/,$fillname);
 my $fill_beg = $unit[1];
 my $fill_end = $unit[1]+$unit[2];
 my @gaps=$fill_branch->daughters; # Load the gaps into an array

###################################### Specific for potentil translocations in the top chain
my $q_beg = $unit[5];
my $q_end = $unit[5] + $unit[6];

if ($unit[-1]=~/top/) {
  if ($unit[2] > 10*$trllen) {
    print TRL "$mamaname\t$fill_beg\t$fill_end\t$unit[3]\t$q_beg\t$q_end\t$unit[4]\t$unit[2]\t$unit[6]\tChrEND\t$unit[-1]\n";
  }  
} elsif ($unit[-1]=~/syn|inv/) {
  if ($unit[2] > $trllen && $$hash_ref{"qFar"} > 2*$unit[2]) {
    print TRL "$mamaname\t$fill_beg\t$fill_end\t$unit[3]\t$q_beg\t$q_end\t$unit[4]\t$unit[2]\t$unit[6]\tIntraChr\t$unit[-1]\n";
  }
} elsif ($unit[-1]=~/nonSyn/) {
  if ($unit[2] > $trllen) {
    print TRL "$mamaname\t$fill_beg\t$fill_end\t$unit[3]\t$q_beg\t$q_end\t$unit[4]\t$unit[2]\t$unit[6]\tInterChr\t$unit[-1]\n";
  }  
}

######################################
# TAG1
 if (@gaps) { # If exist gaps for the fill
  ## TAG2
  foreach my $i (0..$#gaps) {
          my $name0;
          my $name1;
          my $name=$gaps[$i]->name; # the gap line
   ### TAG3
     if ($#gaps>0) {
         if ($i==0) {
             $name0 ="gap_$fill_beg\_0";
             $name1 =$gaps[1]->name;
         } elsif ($i==$#gaps) {
             $name0 =$gaps[$#gaps-1]->name;
             $name1 ="gap_$fill_end\_0";
         } else {
             $name0 =$gaps[$i-1]->name;
             $name1 =$gaps[$i+1]->name;
         }
     } else {
            $name0 ="gap_$fill_beg\_0";
            $name1 ="gap_$fill_end\_0";
     }
   ### TAG3
 my @tmp0 = split (/_/,$name0);
 my @tmp = split (/_/,$name);
 my @tmp1 = split(/_/,$name1);
 my $Tend = $tmp[1]+$tmp[2];
 my $Qend = $tmp[5]+$tmp[6];
 my $forward = $tmp[1] - $tmp0[1] - $tmp0[2];
 my $backward = $tmp1[1] - $tmp[1] - $tmp[2];

 my ($seqT,$seqQ,$seqT5,$seqQ5,$seqT3,$seqQ3);

  if ($tmp[2]<1000) {
      $seqT = substr(${$tseq}{$mamaname},$tmp[1],$tmp[2]+1);
  } else {
      my $seqT1 = substr(${$tseq}{$mamaname},$tmp[1],500);
      my $seqT2 = substr(${$tseq}{$mamaname},$Tend-500,500);
      $seqT = join ("NNNNNNNNNN",($seqT1,$seqT2));
  }

  $seqT5 = substr(${$tseq}{$mamaname},$tmp[1]-20,20);
  $seqT3 = substr(${$tseq}{$mamaname},$Tend,20);
            
  if ($tmp[6]<1000) {
      $seqQ = substr(${$qseq}{$tmp[3]},$tmp[5],$tmp[6]+1);
  } else {
      my $seqQ1 = substr(${$qseq}{$tmp[3]},$tmp[5],500);
      my $seqQ2 = substr(${$qseq}{$tmp[3]},$Qend-500,500);
      $seqQ = join ("NNNNNNNNNN",($seqQ1,$seqQ2));
  }
 
  $seqQ5 = substr(${$qseq}{$tmp[3]},$tmp[5]-20,20);
  $seqQ3 = substr(${$qseq}{$tmp[3]},$Qend,20);

  if ($tmp[4]=~/\-/) {
     $seqQ = &RevCom($seqQ);
     $seqQ5 = &RevCom($seqQ5);
     $seqQ3 = &RevCom($seqQ3);
  }
                                     
#### TAG4
my @downfills=$gaps[$i]->daughters;
            
  if (@downfills){ # if gaps were filled
                    
    if ($tmp[2]/($tmp[6]+0.1)>10 and $tmp[6]<13) {
        my $ali=0;
        my @coordinates=();
        my $inv_len=0;
        my $total_cov=0;
                              
			  foreach my $fill (@downfills) {	# For each filled segments
                my $name=$fill->name;  # Filled segments name
                my $hash = ParseFill ($name);
                   $total_cov = $total_cov + $$hash{"target_len"};
            if ($$hash{"type"}=~/syn|inv/ and $$hash{"qFar"} < 2*$tmp[2]) { # Search potential tandom duplications (copy number variations).
                   $inv_len = $inv_len + $$hash{"ali"} if ($$hash{"type"} =~/inv/);
                   $ali= $ali + $$hash{"ali"};
                my $end = $$hash{"query_beg"} + $$hash{"query_len"};
                          push (@coordinates,$$hash{"query_beg"});
                          push (@coordinates,$end);
            }
        }

        my @array = sort {$a <=> $b} @coordinates;  
                          
		    if ($inv_len > ($tmp[6]+$tmp[2])/4 and $inv_len> $CNVlen) {
                print CNV "$mamaname\t$tmp[1]\t$Tend\tInvCNV\t$tmp[3]\t$array[0]\t$array[-1]\t$tmp[2]\t$tmp[4]\t$forward\t$backward\t$seqT5\t$seqT3\t$seqQ5\t$seqQ3\t$seqT\t$seqQ\n";
        } elsif ($ali/$tmp[2] > $ali_rate and $ali > $CNVlen) {
                print CNV "$mamaname\t$tmp[1]\t$Tend\tCNV\t$tmp[3]\t$array[0]\t$array[-1]\t$tmp[2]\t$tmp[4]\t$forward\t$backward\t$seqT5\t$seqT3\t$seqQ5\t$seqQ3\t$seqT\t$seqQ\n";
        } elsif ($tmp[2] <= $DELlen) {
             if ($total_cov/$tmp[2] < 0.1 and $tmp[2]> 10000) {
			          print DEL "$mamaname\t$tmp[1]\t$Tend\tDEL_PAV\t$tmp[3]\t$tmp[5]\t$Qend\t$tmp[2]\t$tmp[4]\t$forward\t$backward\t$seqT5\t$seqT3\t$seqQ5\t$seqQ3\t$seqT\t$seqQ\n";
             } else {
                print DEL "$mamaname\t$tmp[1]\t$Tend\tDEL\t$tmp[3]\t$tmp[5]\t$Qend\t$tmp[2]\t$tmp[4]\t$forward\t$backward\t$seqT5\t$seqT3\t$seqQ5\t$seqQ3\t$seqT\t$seqQ\n";
             }
                         
       } else {
             if ($total_cov/$tmp[2] < 0.1 and $tmp[2]> 10000) {
                print CXS "$mamaname\t$tmp[1]\t$Tend\tCMPLX_PAV\t$tmp[3]\t$tmp[5]\t$Qend\t$tmp[2]\t$tmp[4]\t$forward\t$backward\t$seqT5\t$seqT3\t$seqQ5\t$seqQ3\t$seqT\t$seqQ\n";
             } else {
                print CXS "$mamaname\t$tmp[1]\t$Tend\tCMPLX\t$tmp[3]\t$tmp[5]\t$Qend\t$tmp[2]\t$tmp[4]\t$forward\t$backward\t$seqT5\t$seqT3\t$seqQ5\t$seqQ3\t$seqT\t$seqQ\n";
             }    
       }

   } else {
         my $ali=0;
         my @coordinates=();
         my $inv_len=0;
         my $total_cov=0;
                          
        foreach my $fill (@downfills) { # For each filled segments
                my $name=$fill->name;  # Filled segments name
                my $hash = ParseFill ($name);
                  $total_cov = $total_cov + $$hash{"target_len"};
                  
                  if ($$hash{"type"} =~/inv|syn/ and $$hash{"qFar"} < 2*$tmp[2] ){ # Define inversions
                      $inv_len = $inv_len + $$hash{"ali"} if ($$hash{"type"} =~/inv/);
                      $ali= $ali + $$hash{"ali"};
                   my $end = $$hash{"query_beg"} + $$hash{"query_len"};
                      push (@coordinates,$$hash{"query_beg"});
                      push (@coordinates,$end);
                  }

  ######## Specially for inversions
        if ($$hash{"type"} =~/inv/ && $$hash{"qFar"} < 2*$tmp[2] && ($$hash{"qOver"} >= $INVlen or ($$hash{"target_len"} > 2*$INVlen and $$hash{"query_len"} > 2*$INVlen)) && $$hash{"ali"} >= $INVlen){  # Define inversions
                my $sub_Tbeg = $$hash{"target_beg"};
                my $sub_Tend = $$hash{"target_beg"} + $$hash{"target_len"};
                my $sub_Qbed = $$hash{"query_beg"};
                my $sub_Qend = $$hash{"query_beg"} + $$hash{"query_len"};
                my $len = $$hash{"target_len"};                                
           print INV "$mamaname\t$sub_Tbeg\t$sub_Tend\tINV\t$tmp[3]\t$sub_Qbed\t$sub_Qend\t$len\t$tmp[4]\t$forward\t$backward\t$seqT5\t$seqT3\t$seqQ5\t$seqQ3\t$seqT\t$seqQ\n";
        }
  ########
  
      }
   my @array = sort {$a <=> $b} @coordinates;

	 if ($inv_len > ($tmp[6]+$tmp[2])/4 and $tmp[2]/($tmp[6]+0.01)>10 and $inv_len > $CNVlen) {
	       print CNV "$mamaname\t$tmp[1]\t$Tend\tinvCNV\t$tmp[3]\t$array[0]\t$array[-1]\t$tmp[2]\t$tmp[4]\t$forward\t$backward\t$seqT5\t$seqT3\t$seqQ5\t$seqQ3\t$seqT\t$seqQ\n";
	 } elsif ($ali/$tmp[2] > $ali_rate and $tmp[2]/($tmp[6]+0.01)>10 and $ali >$CNVlen) {
         print CNV "$mamaname\t$tmp[1]\t$Tend\tCNV\t$tmp[3]\t$array[0]\t$array[-1]\t$tmp[2]\t$tmp[4]\t$forward\t$backward\t$seqT5\t$seqT3\t$seqQ5\t$seqQ3\t$seqT\t$seqQ\n";
	 } else {
      if ($total_cov/$tmp[2] < 0.1 and $tmp[2]> 10000) {
	       print CXS "$mamaname\t$tmp[1]\t$Tend\tCMPLX_PAV\t$tmp[3]\t$tmp[5]\t$Qend\t$tmp[2]\t$tmp[4]\t$forward\t$backward\t$seqT5\t$seqT3\t$seqQ5\t$seqQ3\t$seqT\t$seqQ\n";
      } else {
         print CXS "$mamaname\t$tmp[1]\t$Tend\tCMPLX\t$tmp[3]\t$tmp[5]\t$Qend\t$tmp[2]\t$tmp[4]\t$forward\t$backward\t$seqT5\t$seqT3\t$seqQ5\t$seqQ3\t$seqT\t$seqQ\n";
      }
   }     
 }

} else { # if gaps were not filled
  if ($tmp[2]/($tmp[6]+0.01)>10 and $tmp[6]<13) {  # Insertion
        print DEL "$mamaname\t$tmp[1]\t$Tend\tDEL_PAV\t$tmp[3]\t$tmp[5]\t$Qend\t$tmp[2]\t$tmp[4]\t$forward\t$backward\t$seqT5\t$seqT3\t$seqQ5\t$seqQ3\t$seqT\t$seqQ\n";
  } else {
        print CXS "$mamaname\t$tmp[1]\t$Tend\tCMPLX_PAV\t$tmp[3]\t$tmp[5]\t$Qend\t$tmp[2]\t$tmp[4]\t$forward\t$backward\t$seqT5\t$seqT3\t$seqQ5\t$seqQ3\t$seqT\t$seqQ\n";
  }
   } #### TAG4
  } ## TAG2
 } # TAG1
}

##################
sub RevCom {
my $seq = shift;
my $revcomp = reverse $seq;
$revcomp =~ tr/ATGCNatgcn/TACGNtacgn/;
return $revcomp;
}

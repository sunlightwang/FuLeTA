#!/usr/bin/perl -w
#
########################################
#
# File Name:
#   IsoPairwiseCmp.pl
# 
# Description:
#   
# 
# Usage:
#   
# 
# Author:
#   Xi Wang, Xi.Wang@mdc-berlin.de
# 
# Date:
#   Wed Dec 30 15:54:55 CET 2015
#
########################################

use strict;
use List::MoreUtils qw(uniq);
use List::Util qw(min max sum);
use Blocks;

my $min_5UTR_diff_len = 100;
my $min_3UTR_diff_len = 100;

my $usage = "$0 <isoform.refflat> <isoform_pairwise_comparison_results>\n";
my $infile = shift || die $usage;
my $outfile = shift || die $usage;
open(IN, $infile) || die "Can't open $infile for reading!\n";
open(OUT, ">$outfile") || die "Can't open $outfile for writing!\n";

## save new and known isoform into %isofs
my %isofs; 
my $num = 0; 
while(<IN>) {
  chomp;
  my @a = split;
  my $isof_key = "$a[0]";
  push @{$isofs{$isof_key}}, join "\t", @a[0..10];
}

foreach my $gl (sort keys %isofs) {
  my @isoforms = @{$isofs{$gl}};
  print OUT join "\t", ($gl, "Only_One_Isoform"), "\n" if(@isoforms == 1);
  for(my $i=0; $i<@isoforms; $i++) { # for each new isoform in the gene locus
    my @isof1 = split /\t/, $isoforms[$i];
    for(my $j=($i+1); $j<@isoforms; $j++) { 
      my @isof2 = split /\t/, $isoforms[$j];
      if( Blocks::refflatSize(Blocks::intersectRefflat($isoforms[$i], $isoforms[$j])) == 0 ) {
        print OUT join "\t", ($gl, "Nonoverlap_Isoforms", $isof1[1].":".$isof2[1]), "\n" ;
        next;
      }
      my $tmp = Blocks::cmpIsoform_refflat($isoforms[$i], $isoforms[$j], $min_5UTR_diff_len, $min_3UTR_diff_len); 
      my $cmp_rst = Blocks::isoformCmp_spliceType(Blocks::format_isoformCmpRst($tmp));
      print OUT join "\t", ($gl, $cmp_rst, $isof1[1].":".$isof2[1]), "\n" ;
    }
  }
}


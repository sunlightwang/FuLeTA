#!/usr/bin/perl -w
#
########################################
#
# File Name:
#   cmpIso2known_v2.pl
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

my $usage = "$0 <new_isoform.refflat> <known_isoform.refflat> <new_isoform_comments.refflat>\n";
my $infile1 = shift || die $usage;
my $infile2 = shift || die $usage;
my $outfile = shift || die $usage;
open(NEW, $infile1) || die "Can't open $infile1 for reading!\n";
open(KNOWN, $infile2) || die "Can't open $infile2 for reading!\n";
open(OUT, ">$outfile") || die "Can't open $outfile for writing!\n";

## save new and known isoform into %isofs
my %isofs; 
my $num = 0; 
while(<NEW>) {
  chomp;
  my @a = split;
  $a[2] = "chr".$a[2] unless($a[2] =~ /^chr/);
  next unless ($a[2] =~ /^chr[\dMXY]+$/);
  $num ++; 
  my $isof_key = "$a[2]_$a[3]_$a[4]_$a[5]_${num}_new";
  $isofs{$isof_key} = join "\t", @a[0..10];
}

while(<KNOWN>) {
  chomp;
  my @a = split;
  $a[2] = "chr".$a[2] unless($a[2] =~ /^chr/);
  next unless ($a[2] =~ /^chr[\dMXY]+$/);
  $num ++; 
  my $isof_key = "$a[2]_$a[3]_$a[4]_$a[5]_${num}_known";
  $isofs{$isof_key} = join "\t", @a[0..10];
}
close NEW;
close KNOWN;

## clustering gene loci
my %genes;
my ($chr, $strand, $s, $e) = ("chr0", "+", 0, 0);
my $gene_num = 0;
foreach my $isof_key (sort { my @aa=split "_",$a; my @bb=split "_",$b; $aa[0] cmp $bb[0] || $aa[1] cmp $bb[1] || $aa[2] <=> $bb[2] } keys %isofs) {
  my @a = split /_/, $isof_key; 
  if($chr ne $a[0] || $strand ne $a[1] || $a[2] > $e ) { # next gene loci
    $gene_num ++;
    $chr = $a[0];
    $strand = $a[1];
    $s = $a[2];
    $e = $a[3];
  } else {
    $e = max($e, $a[3]);
  }
  push @{$genes{$gene_num}}, $isof_key;
}

foreach $gene_num (sort {$a <=> $b} keys %genes) {
  my @isof_keys = @{$genes{$gene_num}};
  my @new_isofs = ();
  my @known_isofs = (); 
  my %known_gene_names;
  foreach my $k (@isof_keys) { 
    push @new_isofs, $isofs{$k} if($k =~ /new$/);
    if($k =~ /known$/) {
      push @known_isofs, $isofs{$k} if($k =~ /known$/);
      my @a = split /\t/, $isofs{$k}; 
      $known_gene_names{$a[0]} ++; 
    }
  }
  my $cluster_gene_names = join ",", keys %known_gene_names; 
  next if(@new_isofs == 0); ### only known

  #### new gene loci
  if(@known_isofs == 0) { 
    for(my $i=0; $i<@new_isofs; $i++) { 
      print OUT $new_isofs[$i]."\tNA:NA:-1:Novel_Gene_Locus\n"; 
    }
    next; 
  }

  ### intersect with known gene isoforms
  for(my $i=0; $i<@new_isofs; $i++) { # for each new isoform in the gene locus
    my %intersectBlockNum;
    for(my $j=0; $j<@known_isofs; $j++) { 
      next if( Blocks::refflatSize(Blocks::intersectRefflat($new_isofs[$i], $known_isofs[$j])) == 0 ); ## non-overlap with current known isoform
      my $tmp = Blocks::cmpIsoform_refflat($new_isofs[$i], $known_isofs[$j], $min_5UTR_diff_len, $min_3UTR_diff_len); 
      next if($tmp eq "NA");
      my $tmp_n = Blocks::isoformCmpRst_diffblockNum($tmp);
      my @tmp_a = split /\t/, $known_isofs[$j];
      #my $isof_cmp = join ":", ($cluster_gene_names, $tmp_a[1], Blocks::format_isoformCmpRst($tmp));
      my $isof_cmp = join ":", ($cluster_gene_names, $tmp_a[1], Blocks::isoformCmp_spliceType(Blocks::format_isoformCmpRst($tmp)));
      $intersectBlockNum{$tmp_n} = $isof_cmp;
    }
    ### new isoform with a known gene locus
    if(scalar(keys %intersectBlockNum)==0) {
      print OUT $new_isofs[$i]."\t".$cluster_gene_names.":NA:-1:Nonoverlap_Isoform\n";
      next;
    }
    ### the nearset know isoform
    my $min_s = min(keys %intersectBlockNum); 
    my $isof_cmp = $intersectBlockNum{$min_s};
    print OUT $new_isofs[$i]."\t".$isof_cmp."\n";
  }
}


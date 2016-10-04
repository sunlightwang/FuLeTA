#!/usr/bin/perl -w
#
########################################
#
# File Name:
#   compare.annot.pl
# 
# Description:
#   add comparison comments to the last column
# 
# Usage:
#   
# 
# Author:
#   Xi Wang, Xi.Wang@mdc-berlin.de
# 
# Date:
#   Mon Dec  2 10:20:39 CET 2013
#
########################################

use strict;
use List::Util qw(max min);
use List::MoreUtils qw(uniq);

my $usage = "$0 <FLT.refflat> <ref.refflat> <outfile>\n";
my $infile = shift || die $usage;
my $reffile = shift || die $usage;
my $outfile = shift || die $usage;
open(IN, $infile) || die "Can't open $infile for reading!\n";
open(REF, $reffile) || die "Can't open $reffile for reading!\n";
open(OUT, ">$outfile") || die "Can't open $outfile for writing!\n";

my %isofs;
my %juncs;
my %pasites; 
my %caps;
my $pasite_fuzzy = 50;
my $cap_fuzzy = 50;
my $sj_fuzzy = 10;
my $min_3UTR_diff_len = 100;
my $min_5UTR_diff_len = 100;

my $num = 0;
while(<REF>) {
  chomp;
  my @a = split;
  $a[2] = "chr".$a[2] unless($a[2] =~ /^chr/);
  next unless ($a[2] =~ /^chr[\dMXY]+$/);
  $num ++; 
  my $isof_key = "$a[2]_$a[3]_$a[4]_$a[5]_${num}_ref";
  #print $isof_key,"\n";
  $isofs{$isof_key} = join "\t", @a;
  # junctions
  my @s = split ',', $a[9];
  my @e = split ',', $a[10];
  for(my $i=1; $i<$a[8]; $i++) {
    my $junc_key = join "_", ($a[2], $a[3], $e[$i-1], $s[$i]);
    $juncs{$junc_key} ++;
  }
  # PA sites
  if($a[3] eq '+') {
    for(my $i=($a[5]-$pasite_fuzzy); $i<=($a[5]+$pasite_fuzzy); $i++) {
      my $pakey = "$a[2]_$a[3]_$i";
      $pasites{$pakey} ++ ;
    }
  } else { 
    for(my $i=($a[4]-$pasite_fuzzy); $i<=($a[4]+$pasite_fuzzy); $i++) {
      my $pakey = "$a[2]_$a[3]_$i";
      $pasites{$pakey} ++ ;
    }
  }
  # cap 
  if($a[3] eq '-') {
    for(my $i=($a[5]-$cap_fuzzy); $i<=($a[5]+$cap_fuzzy); $i++) {
      my $pakey = "$a[2]_$a[3]_$i";
      $caps{$pakey} ++ ;
    }
  } else { 
    for(my $i=($a[4]-$cap_fuzzy); $i<=($a[4]+$cap_fuzzy); $i++) {
      my $pakey = "$a[2]_$a[3]_$i";
      $caps{$pakey} ++ ;
    }
  }
}

while(<IN>) {
  chomp;
  my @a = split;
  $a[2] = "chr".$a[2] unless($a[2] =~ /^chr/);
  next unless ($a[2] =~ /^chr[\dMXY]+$/);
  $num ++; 
  my $isof_key = "$a[2]_$a[3]_$a[4]_$a[5]_${num}_FLT";
  #print $isof_key,"\n";
  $isofs{$isof_key} = join "\t", @a;
}
close IN;
close REF;

## clustering gene loci
my %genes;
my ($chr, $strand, $s, $e);
my $first = 1;
my $geneN = 1;
foreach my $isof_key (sort { my @aa=split "_",$a; my @bb=split "_",$b; $aa[0] cmp $bb[0] || $aa[1] cmp $bb[1] || $aa[2] <=> $bb[2] } keys %isofs) {
  #print $isof_key."\n";
  my @a = split "\t", $isofs{$isof_key};
  if($first) {
    $first = 0;
    $chr = $a[2];
    $strand = $a[3];
    $s = $a[4];
    $e = $a[5];
    push @{$genes{$geneN}}, $isof_key;
  } else {
    if($chr ne $a[2] || $strand ne $a[3] || $a[4] > $e ) { # new gene loci
      $geneN ++;
      $chr = $a[2];
      $strand = $a[3];
      $s = $a[4];
      $e = $a[5];
    } else {
      $e = max($e, $a[5]);
    }
    push @{$genes{$geneN}}, $isof_key;
  }
}

## comparing FLT and ref
foreach my $g_n (sort {$a <=> $b} keys %genes) {
  my @isof_keys = @{$genes{$g_n}};
  my @ref_isofs = ();
  my @FLT_isofs = (); 
  foreach my $k (@isof_keys) { 
    push @ref_isofs, $isofs{$k} if($k =~ /ref$/);
    push @FLT_isofs, $isofs{$k} if($k =~ /FLT$/);
  }
  next if($#FLT_isofs == -1); ### only ref
  if($#ref_isofs == -1) { ### new gene loci
    for(my $i=0; $i<@FLT_isofs; $i++) { 
      print OUT $FLT_isofs[$i]."\tNovel_Gene_Loci\n"; 
    }
    next; 
  }

  for(my $i=0; $i<@FLT_isofs; $i++) {
    my $tag = "";  
    # first check 5'cap, 3'pas, splicing
    # if novel, add to $tag
    # if no novel, then check if annotated
    # if annotated, add to $tag
    # if length $tag == 0, then novel combination. 

    my @a = split "\t",$FLT_isofs[$i];
    # for caps
    if($a[3] eq '-') {
      my $capkey = "$a[2]_$a[3]_$a[5]";
      if(! exists $caps{$capkey}) { 
        $tag .= "Novel_Cap|";
      }
    } else {
      my $capkey = "$a[2]_$a[3]_$a[4]";
      if(! exists $caps{$capkey}) {
        $tag .= "Novel_Cap|";
      }
    }
    # for PA sites
    if($a[3] eq '+') {
      my $pakey = "$a[2]_$a[3]_$a[5]";
      if(! exists $pasites{$pakey}) { 
        $tag .= "Novel_PAS|";
      }
    } else {
      my $pakey = "$a[2]_$a[3]_$a[4]";
      if(! exists $pasites{$pakey}) {
        $tag .= "Novel_PAS|";
      }
    }
    # for junctions 
    my @s = split ',', $a[9];
    my @e = split ',', $a[10];
    for(my $ii=1; $ii<$a[8]; $ii++) {
      my $junc_key; 
      my $existing = 0; 
      for(my $j=0; $j<=$sj_fuzzy;$j++) { 
        $junc_key = join "_", ($a[2], $a[3], $e[$ii-1]+$j, $s[$ii]+$j);
        if(exists $juncs{$junc_key}) { 
          $existing = 1; 
          last;
        }
        $junc_key = join "_", ($a[2], $a[3], $e[$ii-1]-$j, $s[$ii]-$j);
        if(exists $juncs{$junc_key}) { 
          $existing = 1; 
          last;
        }
      }
      if(! $existing) {
        $tag .= "Novel_SJC|";
        last;
      }
    }
    unless(length $tag == 0) {
      print OUT $FLT_isofs[$i]."\t".$tag."\n"; 
      next;
    }

    ### for annotated isoforms
    my $j; 
    for($j=0; $j<@ref_isofs; $j++) { 
      my @iso1 = split "\t", $FLT_isofs[$i];
      my @iso2 = split "\t", $ref_isofs[$j];
      next unless($iso1[8] == $iso2[8]); 
      if($iso1[3] eq '+') {
        next unless(abs($iso1[4] - $iso2[4]) < $min_5UTR_diff_len);
        next unless(abs($iso1[5] - $iso2[5]) < $min_3UTR_diff_len);
      } else {
        next unless(abs($iso1[4] - $iso2[4]) < $min_3UTR_diff_len);
        next unless(abs($iso1[5] - $iso2[5]) < $min_5UTR_diff_len);
      }
      my @es1 = split ",", $iso1[9];
      my @ee1 = split ",", $iso1[10];
      my @es2 = split ",", $iso2[9];
      my @ee2 = split ",", $iso2[10];
      my $ii;
      for($ii=1; $ii<$iso1[8]; $ii++) {
        last if(abs($ee1[$ii-1] - $ee2[$ii-1]) > $sj_fuzzy || abs($es1[$ii] - $es2[$ii]) > $sj_fuzzy || ($ee1[$ii-1] - $ee2[$ii-1] - $es1[$ii] + $es2[$ii])!=0); 
      }
      next unless($ii == $iso1[8]);  ## unmatched
      last; 
    }
    if($j<@ref_isofs) {
      print OUT $FLT_isofs[$i]."\tAnnotated\n"; 
      next;
    }
    ## remaining: comb of known RNA processing events
    print OUT $FLT_isofs[$i]."\tNovel_Comb\n"; 
  }
}

close OUT;

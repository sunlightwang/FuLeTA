#!/usr/bin/perl -w
#
########################################
#
# File Name:
#   ORFpred2Bed12.pl
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
#   Thu Oct  1 13:47:12 CEST 2015
#
########################################

use strict;
my $usage = "$0 <isoform.bed12> <predicted_ORF.start_end_on_isoform> <isoform_wORF.bed12> <ORF.bed12>\n";
my $infile1 = shift || die $usage;
my $infile2 = shift || die $usage;
my $outfile1 = shift || die $usage;
my $outfile2 = shift || die $usage;

my %RNA; 
open(IN, $infile1) || die "Can't open $infile1 for reading!\n";
while(<IN>) { 
  chomp; 
  my @a = split; 
  my @b = split /\|/, $a[3];
  my $k = join "|", @b[1..1];
  $RNA{$k} = $_; 
} 
close IN; 

open(IN, $infile2) || die "Can't open $infile2 for reading!\n";
open(OUT1, ">$outfile1") || die "Can't open $outfile1 for writing!\n";
open(OUT2, ">$outfile2") || die "Can't open $outfile2 for writing!\n";

while(<IN>) { 
  chomp;
  my @x = split; 
  my @y = split /\|/, $x[0];
  if(! exists $RNA{$y[0]}) { print STDERR "WARNING: Isoform ID $y[0] cannot be found in $infile1\n"; next} 
  my @a = split /\t/, $RNA{$y[0]}; 
  my @b = split /\t/, $RNA{$y[0]}; 
  my @sizes = split ",", $a[10];
  my @starts = split ",", $a[11];
  my ($cur_s_bl, $cur_e_bl, $th_s, $th_e);
  if($a[5] eq '+') {  # strand +
    $cur_s_bl = 0;
    $cur_e_bl = 0; 
    $th_s = $a[1] + $x[3]; 
    $th_e = $a[1] + $x[4];
    # for start
    while(($th_s - $a[1]) >= ($starts[$cur_s_bl] + $sizes[$cur_s_bl])) {
      $cur_s_bl ++;
      last if($cur_s_bl == $a[9]);
      $th_s = $th_s + $starts[$cur_s_bl] - $starts[$cur_s_bl-1] - $sizes[$cur_s_bl-1];
    }
    # for end
    while(($th_e - $a[1]) > ($starts[$cur_e_bl] + $sizes[$cur_e_bl])) {
      $cur_e_bl ++;
      last if($cur_e_bl == $a[9]);
      $th_e = $th_e + $starts[$cur_e_bl] - $starts[$cur_e_bl-1] - $sizes[$cur_e_bl-1];
    }
  } else { # strand -
    $cur_s_bl = $a[9] - 1;
    $cur_e_bl = $a[9] - 1;
    $th_e = $a[2] - $x[3]; 
    $th_s = $a[2] - $x[4];
    # for start
    while(($th_s - $a[1]) < $starts[$cur_s_bl]) {
      $cur_s_bl --;
      last if($cur_s_bl == -1);
      $th_s = $th_s - $starts[$cur_s_bl+1] + $starts[$cur_s_bl] + $sizes[$cur_s_bl];
    }
    # for end
    while(($th_e - $a[1]) <= $starts[$cur_e_bl]) {
      $cur_e_bl --;
      last if($cur_e_bl == -1);
      $th_e = $th_e - $starts[$cur_e_bl+1] + $starts[$cur_e_bl] + $sizes[$cur_e_bl];
    }
  }
  if($cur_e_bl < $cur_s_bl || $cur_e_bl == $a[9] || $cur_s_bl == $a[9] || $cur_e_bl == -1 || $cur_s_bl == -1) {
    print STDERR "ORF start/end outside the gene region: $_\n"; 
    next;
  }

  ## OUT1
  $a[6] = $th_s;
  $a[7] = $th_e; 
  print OUT1 join "\t", @a;
  print OUT1 "\n";

  ## OUT2
  $b[1] = $th_s;
  $b[2] = $th_e;
  $b[6] = $th_s;
  $b[7] = $th_e;
  $b[9] = $cur_e_bl - $cur_s_bl + 1; 
  my (@abs_starts, @abs_ends, @b_start, @b_size); 
  for(my $i=$cur_s_bl; $i<=$cur_e_bl; $i++) { 
    $abs_starts[$i-$cur_s_bl] = $a[1] + $starts[$i];
    $abs_ends[$i-$cur_s_bl] = $a[1] + $starts[$i] + $sizes[$i];
  } 
  $abs_starts[$cur_s_bl-$cur_s_bl] = $th_s;
  $abs_ends[$cur_e_bl-$cur_s_bl] = $th_e;
  for(my $i=0; $i<$b[9]; $i++) { 
    $b_start[$i] = $abs_starts[$i] - $b[1];
    $b_size[$i] = $abs_ends[$i] - $abs_starts[$i];
  } 
  $b[10] = join ",", @b_size; 
  $b[11] = join ",", @b_start; 
  print OUT2 join "\t", @b;
  print OUT2 "\n";
}

close IN;
close OUT1;
close OUT2;

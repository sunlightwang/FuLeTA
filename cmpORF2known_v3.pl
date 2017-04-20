#!/usr/bin/perl -w
#
########################################
#
# File Name:
#   cmpIso2known.pl
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

my $NtermTruncBaseAllow = 30;  #1. output isoform with thick; maybe also output ncRNA

my $usage = "$0 <FLT_wORF.bed> <ORF_annotated.bed> <out.bed>\n";
my $infile1 = shift || die $usage;
my $infile2 = shift || die $usage;
my $outfile = shift || die $usage;

my (%GL_ORF, %GL_FLT);
open(IN, $infile1) || die "Can't open $infile1 for reading!\n";
while(<IN>) { 
  chomp; 
  my @a = split; 
  next if($a[6] == $a[7]);
  my @b = split /\|/, $a[3]; 
  my $k = $b[0]."_".$a[0]."_".$a[5]; 
  my $ORF = Blocks::thickBed2CDSbed($_);
  push @{$GL_FLT{$k}}, $_;
  push @{$GL_ORF{$k}}, $ORF;
}
close IN;

my %GL_start_end; 
foreach my $g (keys %GL_FLT) { 
  my ($gl, $chr, $strand) = split /_/, $g; 
  my @FLTs = @{$GL_FLT{$g}}; 
  my ($min_s, $max_e, @s, @e); 
  for(my $i=0; $i<@FLTs; $i++) { 
    my @a = split /\t/, $FLTs[$i];
    push @s, $a[1]; 
    push @e, $a[2]; 
  }
  $min_s = min @s;
  $max_e = max @e; 
  $GL_start_end{$g} = join "_", ($chr, $min_s, $max_e); 
}

my %aORF; 
open(IN, $infile2) || die "Can't open $infile2 for reading!\n";
while(<IN>) {
  chomp; 
  my @a = split;
  my $k = $a[0]."_".$a[5]; 
  push @{$aORF{$k}}, $_;
} 
close IN;

# new (no aORF), up/down/intron (no intersection with aORF), fusion, find the nearest and cmp to it
# fusion: should check two intersect isoform cannot be overlapped with any aORF 
# priority: new, fusion, same, up-stream/down-stream/intron, overlap: N, C, in-frame/diff-frame isoform
open(OUT, ">$outfile") || die "Can't open $outfile for writing!\n";
foreach my $g (keys %GL_ORF) { 
  my ($gl, $chr, $strand) = split /_/, $g; 
  my @ORFs = @{$GL_ORF{$g}}; 
  # get aORFs overlaping with ORFs from a GL
  my $k = $chr."_".$strand;
  my @aORFs = @{$aORF{$k}}; 
  my @aORFs_in = Blocks::bedLinesInRegion(\@aORFs, $GL_start_end{$g});

  #### novel: no aORFs_in
  if(scalar(@aORFs_in) == 0) {
    for(my $i=0; $i<@ORFs; $i++) { 
      print OUT $ORFs[$i], "\tNA:NA:Novel-ORF\n"; 
    }
    next; 
  }
  
  my $aORF_range = Blocks::bedLinesRange(\@aORFs_in); 

  my %aORF_gene_names;
  foreach my $k (@aORFs_in) { 
    my @a = split /\t/, $k; 
    $aORF_gene_names{$a[3]} ++; 
  }
  my $cluster_gene_names = join ",", keys %aORF_gene_names; 

  ### intersect with known gene isoforms
  for(my $i=0; $i<@ORFs; $i++) { # for each new ORF in the gene locus
    my $ORF_size_FL = Blocks::bed12Size($ORFs[$i]);
    my %intersectBlockNum;
    my %intersectBed; 
    for(my $j=0; $j<@aORFs_in; $j++) { 
      my $tmp_bed = Blocks::intersectBed($ORFs[$i], $aORFs_in[$j]);
      my $tmp_size = Blocks::bed12Size($tmp_bed);
      next if( $tmp_size == 0 ); ## non-overlap with current known isoform
      my $tmp = Blocks::cmpIsoform_refflat(Blocks::bedToRefflat($ORFs[$i]), Blocks::bedToRefflat($aORFs_in[$j])); 
      next if($tmp eq "NA");
      #my $tmp_n = Blocks::isoformCmpRst_blockNum($tmp);
      my $tmp_n = Blocks::isoformCmpRst_diffblockNum($tmp);
      my $tmp_k = $tmp_size / ($tmp_n + 1);
      $intersectBlockNum{$tmp_k} = $aORFs_in[$j];
      $intersectBed{$tmp_bed} = $aORFs_in[$j]; 
    }
    ### new isoform with a known gene locus
    if(scalar(keys %intersectBlockNum)==0) {
      my ($tmp_chr, $tmp_s, $tmp_e) = split /_/, $aORF_range;
      my @tmp_a = split /\t/, $ORFs[$i];
      if(($tmp_a[5] eq "+" && $tmp_a[2] < $tmp_s) || ($tmp_a[5] eq "-" && $tmp_a[1] > $tmp_e)) { 
        print OUT $ORFs[$i]."\t".$cluster_gene_names.":NA:ORF-5UTR_".$ORF_size_FL."\n"; #5UTR
      } elsif (($tmp_a[5] eq "-" && $tmp_a[2] < $tmp_s) || ($tmp_a[5] eq "+" && $tmp_a[1] > $tmp_e)) { 
        print OUT $ORFs[$i]."\t".$cluster_gene_names.":NA:ORF-3UTR_".$ORF_size_FL."\n"; #3UTR
      } else {  # lies in intron
        print OUT $ORFs[$i]."\t".$cluster_gene_names.":NA:ORF-Intron_".$ORF_size_FL."\n"; #intron
      }
      next;
    }
    ### check if the detected ORF is the fusion of two ORFs
    my @intersectBeds = keys %intersectBed; 
    my $fusion_flag = 0;
    for(my $q=0; $q<@intersectBeds; $q++) { 
      for(my $p=($q+1); $p<@intersectBeds; $p++) { 
        if(Blocks::bed12Size(Blocks::intersectBed($intersectBeds[$q], $intersectBeds[$p])) == 0) { 
          my $type1 = &cmpORF($ORFs[$i], $intersectBed{$intersectBeds[$q]}); 
          my $type2 = &cmpORF($ORFs[$i], $intersectBed{$intersectBeds[$p]}); 
          #print "$type1\t$type2\n";
          next if($type1 eq "diff-frame-isoform" || $type2 eq "diff-frame-isoform");
          my $oo = 0;
          for(; $oo<@intersectBeds; $oo++) { 
            last if($oo != $q && $oo != $p && 
              Blocks::bed12Size(Blocks::intersectBed($intersectBeds[$oo], $intersectBeds[$p])) * Blocks::bed12Size(Blocks::intersectBed($intersectBeds[$oo], $intersectBeds[$q])) > 0);
          }
          $fusion_flag = 1 if($oo == @intersectBeds);
        }
        last if($fusion_flag);
      }
      last if($fusion_flag);
    }
    if($fusion_flag) { 
      print OUT $ORFs[$i]."\t".$cluster_gene_names.":NA:ORF-Fusion_".$ORF_size_FL."\n"; #fusion
      next;
    }
    ### the nearset know isoform
    my $max_k = max(keys %intersectBlockNum); 
    my $nearest_ORF = $intersectBlockNum{$max_k};
    my @nearest_ORF_a = split /\t/, $nearest_ORF;
    # same ORF
    if(Blocks::bed12Same($nearest_ORF, $ORFs[$i])) { 
      print OUT $ORFs[$i]."\t".$cluster_gene_names.":".$nearest_ORF_a[3].":ORF-Same_".$ORF_size_FL."\n"; # same
      next;
    }
    # others
    my $type = &cmpORF($ORFs[$i], $nearest_ORF); 
    my $ORF_size_annot = Blocks::bed12Size($nearest_ORF);
    $type = "ORF-Same" if($type eq "N-term-truncation" && ($ORF_size_annot - $ORF_size_FL) <= $NtermTruncBaseAllow);
    print OUT $ORFs[$i]."\t".$cluster_gene_names.":".$nearest_ORF_a[3].":".$type."_".$ORF_size_FL."_".$ORF_size_annot."\n"; #other types
  }
} 
close OUT;

sub cmpORF { ## cmp two overlap ORFs || overlap: N, C, in-frame/diff-frame/partial-inframe
  my $a1 = $_[0]; # ORF in study
  my $a2 = $_[1]; # annotated
  #print STDERR "\n$a1\n$a2\n";
  my @aa1 = split /\t/, $a1;
  my @aa2 = split /\t/, $a2;
  my $int_bed = Blocks::intersectBed($a1, $a2);
  #print STDERR "$int_bed\n";
  # N & C-term the same
  if($aa1[1] == $aa2[1] && $aa1[2] == $aa2[2]) { 
    return "internal-deletion" if( Blocks::bed12Same($int_bed, $a1) ) ;
    return "internal-insertion" if( Blocks::bed12Same($int_bed, $a2) ) ;
    return "internal-divergence"; 
  }
  my @iii = split /\t/, $int_bed; 
  # only C-term the same
  if(($aa1[5] eq "+" && $aa1[2] == $aa2[2]) || ($aa1[5] eq "-" && $aa1[1] == $aa2[1])) { 
    return "N-term-extension" if( ( ($aa1[5] eq "+" && $aa1[1] < $aa2[1]) || ($aa1[5] eq "-" && $aa1[2] > $aa2[2] ) ) && Blocks::bed12Same($int_bed, $a2) ); 
    return "N-term-truncation" if( ( ($aa1[5] eq "+" && $aa1[1] > $aa2[1]) || ($aa1[5] eq "-" && $aa1[2] < $aa2[2] ) ) && Blocks::bed12Same($int_bed, $a1) ); 
    my $int_start = ($aa1[5] eq "+")?$iii[1]:($iii[2]-1);
    return "N-term-divergence" if(&pos2bed12end($int_start, $a1) == &pos2bed12end($int_start, $a2));
  }
  # only N-term the same
  if(($aa1[5] eq "+" && $aa1[1] == $aa2[1]) || ($aa1[5] eq "-" && $aa1[2] == $aa2[2])) {
    my $int_end = ($aa1[5] eq "+")?($iii[2]-1):$iii[1];
    #print STDERR join "\n", ($a1,$a2,$int_end), "\n";
    return "C-term-divergence" if(&pos2bed12start($int_end, $a1) == &pos2bed12start($int_end, $a2));
  }
  # isoforms
  my @iii_s = split /,/, $iii[11];
  my $n_diff_frame = 0; 
  for(my $i=0; $i<$iii[9]; $i++) { 
    my $pos = $iii[1] + $iii_s[$i]; 
    my $pos_dist1 =  &pos2bed12start($pos, $a1);
    my $pos_dist2 =  &pos2bed12start($pos, $a2);
    #print STDERR "$pos_dist1 \t $pos_dist2\n";
    $n_diff_frame ++ unless(($pos_dist1 % 3) == ($pos_dist2 % 3));
  }
  if($n_diff_frame == 0) { 
    return "inframe-isoform"; 
  } elsif ($n_diff_frame == $iii[9])  { 
    return "diff-frame-isoform";
  } else { 
    return "partial-inframe-isoform";
  }
}

sub ORF_align_cigar { ## cmp two overlap ORFs || overlap: N, C, in-frame/diff-frame/partial-inframe
  my $a1 = $_[0]; # ORF in study
  my $a2 = $_[1]; # annotated
  #print STDERR "\n$a1\n$a2\n";
  my @aa1 = split /\t/, $a1;
  my @aa2 = split /\t/, $a2;
  my $int_bed = Blocks::intersectBed($a1, $a2);
  #print STDERR "$int_bed\n";
  # isoforms
  my @iii_s = split /,/, $iii[11];
  my $n_diff_frame = 0; 
  for(my $i=0; $i<$iii[9]; $i++) { 
    my $pos = $iii[1] + $iii_s[$i]; 
    my $pos_dist1 =  &pos2bed12start($pos, $a1);
    my $pos_dist2 =  &pos2bed12start($pos, $a2);
    #print STDERR "$pos_dist1 \t $pos_dist2\n";
    $n_diff_frame ++ unless(($pos_dist1 % 3) == ($pos_dist2 % 3));
  }
  if($n_diff_frame == 0) { 
    return "inframe-isoform"; 
  } elsif ($n_diff_frame == $iii[9])  { 
    return "diff-frame-isoform";
  } else { 
    return "partial-inframe-isoform";
  }
}

sub pos2bed12start { 
  my $pos = $_[0]; 
  my @a = split /\t/, $_[1]; 
  my @s = split /,/, $a[11]; 
  my @l = split /,/, $a[10]; 
  my $len = sum(@l); 
  my $i; 
  for($i=0; $i<$a[9]; $i++) { 
    last if($pos>=($a[1]+$s[$i]) && $pos<($a[1]+$s[$i]+$l[$i])); 
  }
  return "NA" if($i == $a[9]); 
  $l[$i] = $pos - ($a[1]+$s[$i]);
  my $dist = sum(@l[0..$i]);
  $dist = $len - $dist if($a[5] eq "-"); 
  return $dist; 
}

sub pos2bed12end { 
  my $pos = $_[0]; 
  my @a = split /\t/, $_[1]; 
  my @s = split /,/, $a[11]; 
  my @l = split /,/, $a[10]; 
  my $len = sum(@l); 
  my $i; 
  for($i=0; $i<$a[9]; $i++) { 
    last if($pos>=($a[1]+$s[$i]) && $pos<($a[1]+$s[$i]+$l[$i])); 
  }
  return "NA" if($i == $a[9]); 
  $l[$i] = $pos - ($a[1]+$s[$i]);
  my $dist = sum(@l[0..$i]);
  $dist = $len - $dist if($a[5] eq "+"); 
  return $dist; 
}


#!/usr/bin/perl -w

package Blocks;

use strict;
use List::MoreUtils qw(uniq);
use List::Util qw(min max sum);

##############
# bedToRefflat
# input: bedline, option: separator in the name field (default: '|')
# output: refflat_line
sub bedToRefflat { 
  my @a = split /\t/, $_[0];
  my ($chr, $s, $e, $genename, $transname, $strand, $ss, $se, $num, $ends, $starts);
  my $sep =  '\|'; 
  $sep = $_[1] if(@_ > 1); 
  $chr = $a[0];
  $s = $a[1];
  $e = $a[2];
  my @b = split $sep, $a[3];
  if(@b > 1) {
    $genename = shift @b;
    $a[3] =~ s/$genename//;
    $a[3] =~ s/$sep//;
  } else { 
    $genename = $a[3];
  }
  $transname = $a[3];
  $strand = $a[5];
  $ss = $a[6];
  $se = $a[7];
  $num = $a[9];
  $starts = "";
  $ends = "";
  my @l_s = split ',',$a[10];
  my @s_s = split ',',$a[11];
  for (my $i=0;$i<$num;$i++) {
    my $cur_s = $s_s[$i] + $s;
    my $cur_e = $cur_s + $l_s[$i];
    $starts = $starts.$cur_s.',';
    $ends = $ends.$cur_e.',';
  }
  return(join "\t", ($genename, $transname, $chr, $strand, $s, $e, $ss, $se, $num, $starts, $ends));
}

##############
# refflatToBed
# input: refflat_line, option: separator in the name field (default: '|')
# output: bed_line
sub refflatToBed { 
  my @a = split /\t/, $_[0];
  my $sep =  "|"; 
  $sep = $_[1] if(@_ > 1); 
  my ($chr, $s, $e, $name, $score, $strand, $ss, $se, $col, $num, $lengths, $starts);
  $score = 0;
  $col = 0;
  $chr = $a[2];
  $s = $a[4];
  $e = $a[5];
  $name = $a[1];
  $name = $a[0].$sep.$a[1] unless($a[1] eq $a[0]);
  $strand = $a[3];
  $ss = $a[6];
  $se = $a[7];
  $num = $a[8];
  $lengths = "";
  $starts = "";
  my @s_s = split /,/, $a[9];
  my @e_s = split /,/, $a[10];
  for (my $i=0;$i<$num;$i++) {
    my $cur_s = $s_s[$i] - $s;
    my $cur_l = $e_s[$i] - $s_s[$i];
    $lengths = $lengths.$cur_l.',';
    $starts = $starts.$cur_s.',';
  }
  return join "\t", ($chr,$s,$e,$name,$score,$strand,$ss,$se,$col,$num,$lengths,$starts);
}

###############
# thickRefflat2CDSrefflat
# input: refflat with thick part representing CDS
# output: refflat of the CDS
sub thickRefflat2CDSrefflat { 
  my @a = split /\t/, $_[0]; 
  my @s = split /,/, $a[9];
  my @e = split /,/, $a[10];
  my $n = $a[8]; 
  my ($i, $j);
  for($i=0; $i<$n; $i++) { 
    last if($a[6] >= $s[$i] && $a[6] < $e[$i]);
  }
  for($j=$n-1; $j>=0; $j--) { 
    last if($a[7] > $s[$j] && $a[7] <= $e[$j]);
  }
  $s[$i] = $a[6];
  $e[$j] = $a[7]; 
  $a[4] = $a[6]; 
  $a[5] = $a[7]; 
  $a[8] = $j - $i + 1;
  $a[9] = join ",", @s[$i..$j];
  $a[10] = join ",", @e[$i..$j];
  return join "\t", @a;
}

###############
sub thickBed2CDSbed { 
  return refflatToBed(thickRefflat2CDSrefflat(bedToRefflat($_[0])));
}

###############
# intersectBed
# input: bed_1, bed_2
# output: intersect bed
sub intersectBed { 
  my $a1 = $_[0];
  my $a2 = $_[1]; 
  my @aa1 = split /\t/, $a1;
  my @aa2 = split /\t/, $a2;
  
  my @o = ('-') x 12; 
  my $ret = join "\t", @o; 
  return $ret unless($aa1[0] eq $aa2[0] && $aa1[5] eq $aa2[5] && $aa1[1] < $aa2[2] && $aa1[2] > $aa2[1]); 

  $o[0] = $aa1[0]; 
  $o[5] = $aa1[5]; 
  my @s1 = split ",",$aa1[11];
  my @s2 = split ",",$aa2[11];
  my @l1 = split ",",$aa1[10];
  my @l2 = split ",",$aa2[10];
  my (@es1, @ee1, @es2, @ee2, @es, @ee); 
  for(my $i=0; $i<$aa1[9]; $i++) { 
    $es1[$i] = $aa1[1] + $s1[$i]; 
    $ee1[$i] = $es1[$i] + $l1[$i]; 
  } 
  for(my $i=0; $i<$aa2[9]; $i++) { 
    $es2[$i] = $aa2[1] + $s2[$i]; 
    $ee2[$i] = $es2[$i] + $l2[$i]; 
  } 
  my @uniq_bnd = sort {$a <=> $b} uniq (@es1, @ee1, @es2, @ee2);
  my (@overlap1, @overlap2, @overlapI);
  my $p1 = 0;
  my $p2 = 0;
  for(my $k=0; $k<$#uniq_bnd; $k++) {
    #p1
    if($p1 >= @es1 || $es1[$p1]>=$uniq_bnd[$k+1]) {
      $overlap1[$k] = 0;
    } else {
      $overlap1[$k] = 1;
      if($ee1[$p1] == $uniq_bnd[$k+1]) {
        $p1 ++;
      }
    }
    #p2
    if($p2 >= @es2 || $es2[$p2]>=$uniq_bnd[$k+1]) {
      $overlap2[$k] = 0;
    } else {
      $overlap2[$k] = 1;
      if($ee2[$p2] == $uniq_bnd[$k+1]) {
        $p2 ++;
      }
    }
    $overlapI[$k] = $overlap1[$k] & $overlap2[$k];
  }
  push @es, $uniq_bnd[0] if($overlapI[0] == 1);
  for(my $k=1; $k<$#uniq_bnd; $k++) {
    push @ee, $uniq_bnd[$k] if($overlapI[$k-1] == 1 && $overlapI[$k] == 0);
    push @es, $uniq_bnd[$k] if($overlapI[$k-1] == 0 && $overlapI[$k] == 1);
  }
  push @ee, $uniq_bnd[$#uniq_bnd] if($overlapI[$#overlapI] == 1);
  $o[9] = scalar @es; 
  return $ret unless($o[9] > 0); 
  $o[1] = $es[0]; 
  $o[2] = $ee[$#ee]; 
  my (@s, @l);
  for(my $i=0; $i<$o[9]; $i++) {
    $s[$i] = $es[$i] - $o[1];
    $l[$i] = $ee[$i] - $es[$i];
  }
  $o[11] = join ",", @s;
  $o[10] = join ",", @l;
  return join "\t", @o; 
}

###############
# intersectRefflat
# input: refflat_1, refflat_2
# output: intersect refflat
sub intersectRefflat { 
  my $a1 = $_[0];
  my $a2 = $_[1]; 
  my @aa1 = split /\t/, $a1;
  my @aa2 = split /\t/, $a2;
  
  my @o = ('-') x 11; 
  my $ret = join "\t", @o; 
  return $ret unless($aa1[2] eq $aa2[2] && $aa1[3] eq $aa2[3] && $aa1[4] < $aa2[5] && $aa1[5] > $aa2[4]); 

  $o[2] = $aa1[2]; 
  $o[3] = $aa1[3]; 
  my @es1 = split ",",$aa1[9];
  my @ee1 = split ",",$aa1[10];
  my @es2 = split ",",$aa2[9];
  my @ee2 = split ",",$aa2[10];
  my (@es, @ee); 
  my @uniq_bnd = sort {$a <=> $b} uniq (@es1, @ee1, @es2, @ee2);
  my (@overlap1, @overlap2, @overlapI);
  my $p1 = 0;
  my $p2 = 0;
  for(my $k=0; $k<$#uniq_bnd; $k++) {
    #p1
    if($p1 >= @es1 || $es1[$p1]>=$uniq_bnd[$k+1]) {
      $overlap1[$k] = 0;
    } else {
      $overlap1[$k] = 1;
      if($ee1[$p1] == $uniq_bnd[$k+1]) {
        $p1 ++;
      }
    }
    #p2
    if($p2 >= @es2 || $es2[$p2]>=$uniq_bnd[$k+1]) {
      $overlap2[$k] = 0;
    } else {
      $overlap2[$k] = 1;
      if($ee2[$p2] == $uniq_bnd[$k+1]) {
        $p2 ++;
      }
    }
    $overlapI[$k] = $overlap1[$k] & $overlap2[$k];
  }
  push @es, $uniq_bnd[0] if($overlapI[0] == 1);
  for(my $k=1; $k<$#uniq_bnd; $k++) {
    push @ee, $uniq_bnd[$k] if($overlapI[$k-1] == 1 && $overlapI[$k] == 0);
    push @es, $uniq_bnd[$k] if($overlapI[$k-1] == 0 && $overlapI[$k] == 1);
  }
  push @ee, $uniq_bnd[$#uniq_bnd] if($overlapI[$#overlapI] == 1);
  $o[8] = scalar @es; 
  return $ret unless($o[8] > 0); 
  $o[4] = $es[0]; 
  $o[5] = $ee[$#ee]; 
  $o[9] = join ",", @es;
  $o[10] = join ",", @ee;
  return join "\t", @o; 
}

###############
# bedSizeDiff
# input: bed_1, bed_2
# output: length absolute value
sub bedSizeDiff { 
  my $a1 = $_[0];
  my $a2 = $_[1]; 
  my @aa1 = split /\t/, $a1;
  my @aa2 = split /\t/, $a2;
  my @l1 = split ",",$aa1[10];
  my @l2 = split ",",$aa2[10];
  return(abs(sum(@l1) - sum(@l2)));
}   

###############
# bedLinesInRegion
# input: bedlines, bed_region (chrN_start_end)
# output: bedlines in the region 
sub bedLinesInRegion { 
  my @ret; 
  my @all = @{$_[0]}; 
  my ($chr, $s, $e) = split /_/, $_[1]; 
  for(my $i=0; $i<@all; $i++) { 
    my @a = split /\t/, $all[$i]; 
    unless($chr ne $a[0] || $a[2] < $s || $a[1] > $e) { 
      push @ret, $all[$i]; 
    }
  }
  return @ret; 
}

###############
# bedLinesInStrandedRegion
# input: bedlines, bed_region (chrN_start_end_strand)
# output: bedlines in the region 
sub bedLinesInStrandedRegion { 
  my @ret; 
  my @all = @{$_[0]}; 
  my ($chr, $s, $e, $str) = split /_/, $_[1]; 
  for(my $i=0; $i<@all; $i++) { 
    my @a = split /\t/, $all[$i]; 
    unless($chr ne $a[0] || $str ne $a[5] || $a[2] < $s || $a[1] > $e) { 
      push @ret, $all[$i]; 
    }
  }
  return @ret; 
}

################
# bedLinesRange
# input: bedlines
# output: chrN_start_end
sub bedLinesRange { 
  my @all = @{$_[0]};
  my (@s, @e); 
  my $chr; 
  for(my $i=0; $i<@all; $i++) {
    my @a = split /\t/, $all[$i];
    $chr = $a[0] if($i==0);
    return "NA" if($i > 0 && $chr ne $a[0]);
    push @s, $a[1];
    push @e, $a[2]; 
  } 
  return join "_", ($chr, min(@s), max(@e));
}


################
# bed12Same
# input bed_1, bed_2
# output 1 the same; 0 different
sub bed12Same {
  my $a1 = $_[0];
  my $a2 = $_[1];
  my @aa1 = split /\t/, $a1;
  my @aa2 = split /\t/, $a2;
  return 0 unless ($aa2[9] =~ /[\d]+/ && $aa1[9] =~ /[\d]+/);
  foreach my $i (0,1,2,5,9,10,11) {
    if($i == 10 || $i == 11) {
      $aa1[$i] =~ s/,$//;
      $aa2[$i] =~ s/,$//;
    }
    return 0 if($aa1[$i] ne $aa2[$i]);
  }
  return 1;
}

#################
# bed12Size
# input bed
# output the length of bed
sub bed12Size { 
  my @a = split /\t/, $_[0];
  return 0 if($a[9] eq "-" || $a[9] == 0);
  my @b = split /,/, $a[10]; 
  return sum(@b); 
}

#################
# refflatSize
# input bed
# output the length of bed
sub refflatSize { 
  my @a = split /\t/, $_[0];
  return 0 if($a[8] eq "-" || $a[8] == 0);
  my @s = split /,/, $a[9]; 
  my @e = split /,/, $a[10]; 
  my $l = 0; 
  for(my $i=0; $i<@s; $i++) { 
    $l += $e[$i] - $s[$i]; 
  }
  return $l; 
}

####################
# cmpIsoform_refflat
# input: isoform_new, isoform_known
# output: [overlap1]:[overlap2]:[len1_len2_len3_len4_len5]
## output: [string: 1o0o1i0o1i]:[len1_len2_len3_len4_len5]
sub cmpIsoform_refflat {
  my $i1 = $_[0];
  my $i2 = $_[1];
  my $min_5UTR_diff = 0; 
  $min_5UTR_diff = $_[2] if(@_ > 2);
  my $min_3UTR_diff = 0;
  $min_3UTR_diff = $_[3] if(@_ > 3);
  my @iso1 = split "\t", $i1;
  my @iso2 = split "\t", $i2; 
  return "NA" unless ($iso1[2] eq $iso2[2] && $iso1[3] eq $iso2[3]); 
  return "NA" if($iso1[4] > $iso2[5] || $iso1[5] < $iso2[4]);
  my @es1 = split ",", $iso1[9];
  my @ee1 = split ",", $iso1[10];
  my @es2 = split ",", $iso2[9];
  my @ee2 = split ",", $iso2[10];
 
  ## block boundaries 
  my @uniq_bnd = sort {$a <=> $b} uniq (@es1, @ee1, @es2, @ee2);
  my (@overlap1, @overlap2, @block_size); 
  my $p1 = 0;
  my $p2 = 0;
  for(my $k=0; $k<$#uniq_bnd; $k++) { 
    push @block_size, $uniq_bnd[$k+1] - $uniq_bnd[$k]; 
    # p1 for isoform 1
    if($p1 >= @es1 || $es1[$p1]>=$uniq_bnd[$k+1]) {
      $overlap1[$k] = 0; 
    } else {
      $overlap1[$k] = 1;
      if($ee1[$p1] == $uniq_bnd[$k+1]) {
        $p1 ++;
      }
    }
    # p2 for isoform 2 
    if($p2 >= @es2 || $es2[$p2]>=$uniq_bnd[$k+1]) {
      $overlap2[$k] = 0; 
    } else {
      $overlap2[$k] = 1;
      if($ee2[$p2] == $uniq_bnd[$k+1]) {
        $p2 ++;
      }
    }
  }
  if($iso1[3] eq '-') { 
    @overlap1 = reverse @overlap1;
    @overlap2 = reverse @overlap2;
    @block_size = reverse @block_size;
  }
  ## for 5' and 3' ends: allow a few bases differnece
  if($block_size[0] < $min_5UTR_diff && $overlap1[0] == 0 && $overlap2[0] == 1 && $overlap1[1] == 1 && $overlap2[1] == 1) { 
    shift @overlap1;
    shift @overlap2; 
    shift @block_size;
  }
  if($block_size[0] < $min_5UTR_diff && $overlap1[0] == 1 && $overlap2[0] == 0 && $overlap1[1] == 1 && $overlap2[1] == 1) { 
    shift @overlap1;
    shift @overlap2; 
    $block_size[1] += $block_size[0];
    shift @block_size;
  }
  if($block_size[$#block_size] < $min_3UTR_diff && $overlap1[$#block_size] == 0 && $overlap2[$#block_size] == 1 && $overlap1[$#block_size-1] == 1 && $overlap2[$#block_size-1] == 1) { 
    pop @overlap1;
    pop @overlap2; 
    pop @block_size;
  }
  if($block_size[$#block_size] < $min_3UTR_diff && $overlap1[$#block_size] == 1 && $overlap2[$#block_size] == 0 && $overlap1[$#block_size-1] == 1 && $overlap2[$#block_size-1] == 1) { 
    pop @overlap1;
    pop @overlap2; 
    $block_size[$#block_size-1] += $block_size[$#block_size];
    pop @block_size;
  }
  ## @overlap1, @overlap2 have their values: 1 for exon, 0 for intron 
  return join ":", (join("", @overlap1), join("", @overlap2), join("_",@block_size));
}

#################
# cmpBlocks
# input: bed12_new, bed12_known
# output: [overlap1]:[overlap2]:[len1_len2_len3_len4_len5]
sub cmpBlocks {
  my $i1 = $_[0];
  my $i2 = $_[1];
  my @iso1 = split "\t", $i1;
  my @iso2 = split "\t", $i2; 
  return "NA" unless ($iso1[0] eq $iso2[0] && $iso1[5] eq $iso2[5]); 
  return "NA" if($iso1[1] > $iso2[2] || $iso1[2] < $iso2[1]);
  my @s1 = split ",", $iso1[11];
  my @l1 = split ",", $iso1[10];
  my @s2 = split ",", $iso2[11];
  my @l2 = split ",", $iso2[10];
  my (@es1, @ee1, @es2, @ee2); 
  for(my $i=0; $i<$iso1[9]; $i++) { 
    $es1[$i] = $s1[$i] + $iso1[1];
    $ee1[$i] = $s1[$i] + $l1[$i] + $iso1[1];
  }
  for(my $i=0; $i<$iso2[9]; $i++) { 
    $es2[$i] = $s2[$i] + $iso2[1];
    $ee2[$i] = $s2[$i] + $l2[$i] + $iso2[1];
  }
  ## block boundaries 
  my @uniq_bnd = sort {$a <=> $b} uniq (@es1, @ee1, @es2, @ee2);
  my (@overlap1, @overlap2, @block_size); 
  my $p1 = 0;
  my $p2 = 0;
  for(my $k=0; $k<$#uniq_bnd; $k++) { 
    push @block_size, $uniq_bnd[$k+1] - $uniq_bnd[$k]; 
    # p1 for isoform 1
    if($p1 >= @es1 || $es1[$p1]>=$uniq_bnd[$k+1]) {
      $overlap1[$k] = 0; 
    } else {
      $overlap1[$k] = 1;
      if($ee1[$p1] == $uniq_bnd[$k+1]) {
        $p1 ++;
      }
    }
    # p2 for isoform 2 
    if($p2 >= @es2 || $es2[$p2]>=$uniq_bnd[$k+1]) {
      $overlap2[$k] = 0; 
    } else {
      $overlap2[$k] = 1;
      if($ee2[$p2] == $uniq_bnd[$k+1]) {
        $p2 ++;
      }
    }
  }
  ## @overlap1, @overlap2 have their values: 1 for exon, 0 for intron 
  if($iso1[5] eq '-') { 
    @overlap1 = reverse @overlap1;
    @overlap2 = reverse @overlap2;
    @block_size = reverse @block_size;
  }
  return join ":", (join("", @overlap1), join("", @overlap2), join("_",@block_size));
}

#################
# cmpBlocks_cigar
# input: output from cmpBlocks: [overlap1]:[overlap2]:[len1_len2_len3_len4_len5]
# output: mM_iI_dD (small letters: integer)
sub cmpBlocks_cigar {
  my $in = $_[0];
  my @e = split /:/, $in; 
  my @overlap1 = split //, $e[0];
  my @overlap2 = split //, $e[1];
  my @size = split /_/, $e[2];
  my (@cigar_size, @cigar_ch);
  for(my $i=0; $i<@overlap1; $i++) { 
    next unless($overlap1[$i] || $overlap2[$i]); 
    push @cigar_size, "$size[$i]";
    push @cigar_ch, "M" if($overlap1[$i] && $overlap2[$i]);
    push @cigar_ch, "I" if($overlap1[$i] && !$overlap2[$i]);
    push @cigar_ch, "D" if(!$overlap1[$i] && $overlap2[$i]);
  } 
  ##tidy cigar
  my (@cigar_tidy_size, @cigar_tidy_ch, @cigar); 
  my $idx = 0;
  $cigar_tidy_size[$idx] = $cigar_size[0]; 
  $cigar_tidy_ch[$idx] = $cigar_ch[0]; 
  for(my $i=1; $i<@cigar_size; $i++) { 
    if($cigar_tidy_ch[$idx] eq $cigar_ch[$i]) { 
      $cigar_tidy_size[$idx] += $cigar_size[$i]; 
    } else {
      $idx ++; 
      $cigar_tidy_ch[$idx] = $cigar_ch[$i];
      $cigar_tidy_size[$idx] = $cigar_size[$i];
    }
  } 
  for(my $i=0; $i<@cigar_tidy_size; $i++) {
    $cigar[$i] = "$cigar_tidy_size[$i]$cigar_tidy_ch[$i]";
  }
  return(join "_", @cigar); 
}
#################
# isoformCmpRst_blockNum
# input: [overlap1]:[overlap2]:[len1_len2_len3_len4_len5]
# output: number of len blocks
sub isoformCmpRst_blockNum {
  return "NA" if($_[0] eq "NA");
  my @a = split /:/, $_[0];
  return length $a[0];
}

#################
# isoformCmpRst_diffblockNum
# input: [overlap1]:[overlap2]:[len1_len2_len3_len4_len5]
# output: number of diff blocks
sub isoformCmpRst_diffblockNum {
  return "NA" if($_[0] eq "NA");
  my @a = split /:/, $_[0];
  my @a0 = split //, $a[0]; 
  my @a1 = split //, $a[1]; 
  my @tmp; 
  for my $i (0..$#a0) { 
    push @tmp, abs($a1[$i] - $a0[$i]);
  }
  return sum(@tmp);
}

#################
# format_isoformCmpRst
# input: [overlap1]:[overlap2]:[len1_len2_len3_len4_len5]
# output: [string: 1o0o1i0o1i]:[len1_len2_len3_len4_len5]
sub format_isoformCmpRst { 
  my $inline = $_[0]; 
  my @s = split /:/, $inline; 
  my @s1_org = split //, $s[1]; 
  $s[1] =~ tr/01/oi/;
  my @s0 = split //, $s[0]; 
  my @s1 = split //, $s[1]; 
  my @tmp; 
  my @out; 
  for my $i (0..$#s0) { 
    push @tmp, abs($s1_org[$i] - $s0[$i]);
    push @out, $s0[$i];
    push @out, $s1[$i];
  }
  return join ":", (0, "Annotated", $s[2]) if(sum(@tmp) == 0);
  return join ":", (sum(@tmp), join("", @out), $s[2]);
}

###############
# splice code
###############
#Exon skipping:
#1i0o0i0o1i
#1i0o1o0o1i (annot exon skipping)
#
#Alt 3'SS:
#1i0o0i1i
#1i0o1o1i
#
#Alt 5'SS:
#1i0i0o1i
#1i1o0o1i
#
#Intron retention:
#1i1o1i
#1i0i1i (annot intron retention)
#
#Mut excl exons:
#1i0o1o0o0i0o1i
#1i0o0i0o1o0o1i
#
#Alt first exons:
#^1o0o
#^0i0o
#
#Alt last exons:
#0o1o$
#0o0i$
#
#Tandem 5'UTR
#^1o1i
#^0i1i
#
#Tandem 3'UTR
#1i1o$
#1i0i$
#
####
#Last 2nd exon skipping
#1i0o1o0o1i$
#1i0o1o0o1i[10][io]$
#
#Last intron retention 
#1o1o1i$
#1o1o1i[10][io]$
#
#################
# isoformCmp_spliceType
# input: [string: 1o0o1i0o1i]:[len1_len2_len3_len4_len5]
# output: string of alt splicing type, separated by :
sub isoformCmp_spliceType { 
  my $inline = $_[0]; 
  my @s = split /:/, $inline; 
  return join ":", ($s[0], $s[1]) if($s[1] =~ /Annotated/);
  my @type; 
  my @n;
  ### Tandem UTR can only be matched once
  push @type, "Tandem5UTR" if($s[1] =~ /^1o1i/ || $s[1] =~ /^0i1i/); 
  push @type, "Tandem3UTR" if($s[1] =~ /1i1o$/ || $s[1] =~ /1i0i$/); 
  #
  ### Alt first/last exons
  if($s[1] =~ /^1o0o/ || $s[1] =~ /^0i0o/) { 
    push @type, "AltFirstExons"; 
    # adjust $s[0]
    if($s[1] =~ /^1o0o/) { 
      while($s[1] =~ /^1o0o([10io]+)$/) {$s[0] --; $s[1] = $1}
    } else { 
      while($s[1] =~ /^0i0o([10io]+)$/) {$s[0] --; $s[1] = $1}
    }
    $s[0] ++ if($s[1] =~ /^1i/); # adjust back if start postion at splice site
  }
  if($s[1] =~ /0o1o$/ || $s[1] =~ /0o0i$/) {
    push @type, "AltLastExons" ; 
    # adjust $s[0]
    if($s[1] =~ /0o1o$/) {
      while($s[1] =~ /^([10io]+)0o1o$/) {$s[0] --; $s[1] = $1}
    } else {
      while($s[1] =~ /^([10io]+)0o0i$/) {$s[0] --; $s[1] = $1}
    }
    $s[0] ++ if($s[1] =~ /1i$/); # adjust back if end postion at splice site
  }
  ### ATTENTION: $s[1] may have new value below 
  #
  ### Exon-skipping, intron-retention, alt-splice-site can be matched multi times
  # Exon-skipping
  #push @type, ("ExonSkipping") x scalar(@n) if(@n = $s[1] =~ /1i[10]?[io]?0o0i0o[10]?[io]?1i/g); #  next line: (?=...) is a lookahead assertion
  push @type, ("ExonSkipping") x scalar(@n) if(@n = $s[1] =~ /(?=1i[10]?[io]?0o0i0o[10]?[io]?1i)/g);
  if(@n = $s[1] =~ /(?=1i[10]?[io]?0o0i0o0i0o[10]?[io]?1i)/g) { push @type, ("BiExonSkipping") x scalar(@n); $s[0] = $s[0] - scalar(@n); }
  if(@n = $s[1] =~ /(?=1i[10]?[io]?0o0i0o0i0o0i0o[10]?[io]?1i)/g) { push @type, ("TriExonSkipping") x scalar(@n); $s[0] = $s[0] - 2 * scalar(@n); }
  push @type, ("AnnotExonSkipping") x scalar(@n) if(@n = $s[1] =~ /(?=1i[10]?[io]?0o1o0o[10]?[io]?1i)/g);
  if(@n = $s[1] =~ /(?=1i[10]?[io]?0o1o0o1o0o[10]?[io]?1i)/g) { push @type, ("BiAnnotExonSkipping") x scalar(@n); $s[0] = $s[0] - scalar(@n); } 
  if(@n = $s[1] =~ /(?=1i[10]?[io]?0o1o0o1o0o1o0o[10]?[io]?1i)/g) { push @type, ("TriAnnotExonSkipping") x scalar(@n); $s[0] = $s[0] - 2 * scalar(@n); }
  # Alt splice site
  push @type, ("Alt3SS") x scalar(@n) if(@n = $s[1] =~ /0o[10][io]1i/g); 
  push @type, ("Alt5SS") x scalar(@n) if(@n = $s[1] =~ /1i[10][io]0o/g);
  # Intron retention
  push @type, ("IntronRetention") x scalar(@n) if(@n = $s[1] =~ /(?=1i1o1i)/g); 
  push @type, ("AnnotIntronRetention") x scalar(@n) if(@n = $s[1] =~ /(?=1i0i1i)/g); 
  # Mutually exclusive exons
  if(@n = $s[1] =~ /1i[10]?[io]?0o1o0o0i0o[10]?[io]?1i/g) { push @type, ("MutExclExons") x scalar(@n); $s[0] = $s[0] - scalar(@n); } 
  if(@n = $s[1] =~ /1i[10]?[io]?0o0i0o1o0o[10]?[io]?1i/g) { push @type, ("MutExclExons") x scalar(@n); $s[0] = $s[0] - scalar(@n); } 
  return join ":", ($s[0], "ComplexSplicing", @type) if(@type < $s[0]);
  #print STDERR join ":", ($s[0], @type), "\n" if(@type > $s[0]);
  return join ":", ($s[0], @type);
}

1;

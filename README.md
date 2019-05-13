# FuLeTA

Tookit for full length transcriptome annotation analysis, including isoform diversity, alternative splicing among transcripts, comparison to known transcriptome, splicing coupling, and their consequences on ORF diveristy. 

* This is part of the rat hippocampus transcriptome annotation project, where hybrid sequencing technology (PacBio + Illumina) was applied to get full-length transcriptome landscape for rat hippocampus. 

## Content

  1. Transcript isoform diversity analysis

  2. Alternative splicing events 

  3. Comparison to known transcriptome annotation 

  4. Co-occurrence of alternative RNA processing events analysis

  5. ORF diversity analysis

## Usage

* __IsoDiv.pl__  
  _Description_: comparing each isoform to its nearest isoform in the annotation collection  
  _Usage_: IsoDiv.pl <isoform.refflat> <isoform.refflat w/ results in an additional column>

* IsoPairwiseCmp.pl  -- pairwise comparison of isoforms in the same gene locus 
  Usage: IsoPairwiseCmp.pl <isoform.refflat> <comparing results>


## Contact:
    xi (dot) wang (at) dkfz (dot) de

# FuLeTA

Tookit for full length transcriptome annotation analysis, including isoform diversity, alternative splicing among transcripts, comparison to known transcriptome, splicing coupling, and their consequences on ORF diveristy. 

* This is part of the rat hippocampus transcriptome annotation project, where hybrid sequencing technology (PacBio + Illumina) was applied to get full-length transcriptome landscape for rat hippocampus.  
* All the executables of different functions are based on the exon block analysis, which is very useful in comparison of isoform structure annotation. 

## Content

  1. Transcript isoform diversity analysis

  2. Alternative splicing events 

  3. Comparison to known transcriptome annotation 

  4. Co-occurrence of alternative RNA processing events analysis

  5. ORF diversity analysis

## Usage

* __IsoDiv.pl__   
  _Description_: Comparing each isoform to its nearest isoform in the annotation collection    
  _Usage_: IsoDiv.pl <isoform.refflat> <isoform.refflat w/ results in an additional column>  

* __IsoPairwiseCmp.pl__   
  _Description_: Pairwise comparison for isoform differences in the same gene locus   
  _Usage_: IsoPairwiseCmp.pl <isoform.refflat> <comparing results>  

* __cmpIso2known.pl__    
  _Description_: Comparing one collection of transcript isoforms to an existing (e.g. RefSeq) collection   
  _Usage_: cmpIso2known.pl <new_collection_isoform.refflat> <known_isoform.refflat> <new_collection_isoform.refflat w/ results in an additional column>   
  _Note_: The comparison results would be useful for investigating co-occurrence of alternative RNA processing events  

* __cmpORF2known.pl__    
  _Description_: Comparing one collection of open reading frames (ORFeome) to an existing (e.g. RefSeq) collection   
  _Usage_: cmpORF2known.pl <ORF.bed> <known_ORF.bed> <ORF.bed w/ results in an additional column>     


## Contact
    xi (dot) wang (at) dkfz (dot) de   
    for bug reporting or requiring additional functionality


## Citation
1.  Xi Wang, Xintian You, Jingyi Hou, Julian D. Langer, Fiona Rupprecht, Irena Vlatkovic, Claudia Quedenau, Georgi Tushev, Irina Epstein, Bernhard Schaefke, Wei Sun, Liang Fang, Guipeng Li, Yuhui Hu, Erin M Schuman, Wei Chen. (2019) __Full-length transcriptome reconstruction reveals a large diversity of RNA isoforms and open reading frames in the rat hippocampus__, _under review_. 


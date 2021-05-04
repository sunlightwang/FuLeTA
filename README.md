# FuLeTA

Toolkit for full length transcriptome annotation analysis, including isoform diversity, alternative splicing among transcripts, comparison to known transcriptome, splicing coupling, and their consequences on ORF diveristy. 

* This is part of the rat hippocampus transcriptome annotation project, where hybrid sequencing technology (PacBio + Illumina) was applied to get full-length transcriptome landscape for rat hippocampus.  
* All the executables of different functions are based on the exon block analysis, which is very useful in comparison of isoform structure annotation. 

## Content

  1. Transcript isoform diversity analysis

  2. Alternative splicing events 

  3. Comparison to known transcriptome annotation 

  4. Co-occurrence of alternative RNA processing events analysis

  5. ORF diversity analysis

## Installation
1. Operating system requirement: any OS running perl (>= 5)
2. Download the toolkit: 
```git clone https://github.com/sunlightwang/FuLeTA.git```
3. Changed to the FuLeTA directory and run the scripts: 
```cd FuLeTA```

Expected time: a couple of minutes

## Usage

* __IsoDiv.pl__   
  _Description_: Comparing each isoform to its nearest isoform in the annotation collection    
  _Usage_: IsoDiv.pl <isoform.refflat> <isoform.refflat w/ results in an additional column>  

* __IsoPairwiseCmp.pl__   
  _Description_: Pairwise comparison for isoform differences in the same gene locus  
  _Usage_: IsoPairwiseCmp.pl <isoform.refflat> <comparison results>  

* __cmpIso2known.pl__    
  _Description_: Comparing one collection of transcript isoforms to an existing (e.g. RefSeq) collection   
  _Usage_: cmpIso2known.pl <new_collection_isoform.refflat> <known_isoform.refflat> <new_collection_isoform.refflat w/ results in an additional column>   
  _Note_: The comparison results would be useful for investigating co-occurrence of alternative RNA processing events  

* __cmpORF2known.pl__    
  _Description_: Comparing one collection of open reading frames (ORFeome) to an existing (e.g. RefSeq) collection   
  _Usage_: cmpORF2known.pl <ORF.bed> <known_ORF.bed> <ORF.bed w/ results in an additional column>     


## Examples  

* __IsoDiv.pl__  
```
IsoDiv.pl example_data/FLT.rn6.chr15.refflat example_data/FLT.rn6.chr15.IsoDiv   
```

* __IsoPairwiseCmp.pl__   
```
IsoPairwiseCmp.pl example_data/FLT.rn6.chr15.refflat example_data/FLT.rn6.chr15.IsoPairwiseCmp   
```


* __cmpIso2known.pl__   
```
cmpIso2known_v2.pl example_data/FLT.rn6.chr15.refflat example_data/RefSeq.rn6.chr15.refflat example_data/FLT.rn6.chr15.cmpFLT2RefSeq   
```

* __cmpORF2known.pl__   
```
cmpORF2known_v3.pl example_data/FLT.rn6.chr15.bed example_data/RefSeq_CDS.rn6.chr15.bed example_data/FLT.rn6.chr15.cmpFLTORF2RefSeq  
```
Expected time: a couple of minutes

## Contact
```xi (dot) wang (at) dkfz (dot) de  ```   
for bug reporting or requiring additional functionality


## Citation
1.  __Xi Wang__, Xintian You, Julian D. Langer, Jingyi Hou, Fiona Rupprecht, Irena Vlatkovic, Claudia Quedenau, Georgi Tushev, Irina Epstein, Bernhard Schaefke, Wei Sun, Liang Fang, Guipeng Li, Yuhui Hu, Erin M Schuman, Wei Chen. __Full-length transcriptome reconstruction reveals a large diversity of RNA and protein isoforms in rat hippocampus__. __*Nat Commun*__ 10, 5009 (2019) [doi:10.1038/s41467-019-13037-0](https://doi.org/10.1038/s41467-019-13037-0)

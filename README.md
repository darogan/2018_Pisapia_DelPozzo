# Tristetraprolin/ZFP36 regulates the turnover of autoimmune-associated HLA-DQ mRNAs

**Pisapia Laura<sup>a,‡</sup>, Hamilton S. Russell<sup>b,‡</sup>, D’Agostino Vito<sup>c</sup>, Barba Pasquale<sup>a</sup>, Strazzullo Maria<sup>a</sup>, Provenzano Alessandro<sup>c</sup>, Gianfrani Carmen<sup>d</sup> and Del Pozzo Giovanna<sup>a,§</sup>**

<sup>
<sup>a</sup> Institute of Genetics and Biophysics “Adriano Buzzati Traverso” CNR, Via Pietro Castellino, 111, 80131, Naples, Italy <br>
<sup>b</sup> Centre for Trophoblast Research, Department of Physiology, Development and Neuroscience, University of Cambridge, Downing Site, Cambridge, CB2 3DY<br>
<sup>d</sup> Institute of Protein Biochemistry-CNR, Via Pietro Castellino, 111, 80131, Naples, Italy<br>
<sup>‡</sup> Equal contribution<br>
<sup>§</sup> Corresponding author <br>
</sup>

## Citation ##

Pisapia, L., Hamilton, R.S., D’Agostino, V, Pasqualea, B., Strazzullo, M., Provenzano, A., Gianfrani, C. & Del Pozzo, G. (2018) Tristetraprolin/ZFP36 regulates the turnover of autoimmune-associated HLA-DQ mRNAs [bioRxiv](https://doi.org/10.1101/337907)

## Abstract ##

We have previously demonstrated that the expression of HLA class II transcripts is regulated by the binding of a ribonucleoprotein complex that affects mRNA processing. We identified protein components of a complex binding transcripts encoding the HLA-DR molecule. Here we aimed to verify if the same RNA binding proteins interact with 3’UTR of messengers encoding the HLA-DQ isotype Specifically, we focused on the HLA-DQ2.5 molecule, expressed on the surface of antigen presenting cells, and representing the main susceptibility factor for celiac disease. This molecule, encoded by DQA1*05 and DQB1*02 alleles, presents the antigenic gluten peptides to CD4+ T lymphocytes, activating the autoimmune response.
Here, we identified an additional component of the RNP complex, Tristetraprolin (TTP) or ZFP36, a zinc- finger protein, widely described as a factor modulating mRNA stability. TTP shows high affinity binding to 3’UTR of CD-associated DQA1*05 and DQB1*02 alleles, in contrast to lower affinity binding to DQA1*01 and DQB1*05 non-CD associated alleles. Our in silico analysis, confirmed by molecular depletion of protein and HLA mRNA quantification, demonstrates that TTP specifically modulates the stability of the transcripts associated with celiac disease. Our work demonstrates, for the first time, that proteins of the RNP complex, affecting the processing of transcripts, interact with allele-specific transcripts.  


### Figure 2 ###

RNAVienna (Lorenz et al, 2011) dot bracket notation structures were rendered into 2D images using FORNA (Kerpedjiev et al, 2015). Additional motifs were added from FOLDALIGN (Havgaard et al, 2007) and the `AURichness.pl` script described below.

### Figure 3 ###

To calculate the AU rich motifs in the riboprobes a stand alone Perl script takes the sequences and produces frequency tables. These tables are then inputted into the R script to produce the plots as seen in Figure 3.

The following commands run from the same directory will produce frequency tables and plots from them using R.

Generate frequency tables:

    perl AURichness.pl

Example frequency table:

| Sequence  | Dinucleotide | Obs_freq           | Exp_freq           |
| ----------| ------------ | ------------------ | -------------------|
| DQB102    | AA           | 0.0144230769230769 | 0.0313408575810993 |
| DQB102    | AC           | 0.0673076923076923 | 0.0592935143426204 |
| DQB102    | AG           | 0.0576923076923077 | 0.0372702090153614 |

Example density table:

| position | DQA101_di | DQA101_tri | DQA101_quad | DQA101_penta | DQA101_hexa | DQA105_di | DQA105_tri | DQA105_quad | DQA105_penta | DQA105_hexa
|----------|----------|----------|----------|----------|----------|----------|----------|----------|----------|----------|
| 0|0|0|0|0|0|0|0|0|0|0|
| 1|0|0|0|0|0|0|0|0|1|0|
| 2|0|0|0|0|0|0|0|0|0|0|

Create plots from frequency tables:

    Rscript AURichAnalysis.R

Example AU Rich Motif Observed Vs Expected Frequency plot

![AUMotif_ObsExp](AUMotif_ObsExp.png?raw=true=100x)

Examples AU Rich Motif Position plot

![AUMotif_Positions](AUMotif_Positions.png?raw=true=100x)


## References ##

Havgaard JH, Torarinsson E, Gorodkin J. Fast pairwise structural RNA alignments by pruning of the dynamical programming matrix. PLoS Comput Biol. 2007 Oct;3(10):1896-908. doi: 10.1371/journal.pcbi.0030193. PubMed PMID: 17937495; PubMed Central PMCID: PMCPMC2014794

Kerpedjiev P, Hammer S, Hofacker IL. Forna (force-directed RNA): Simple and effective online RNA secondary structure diagrams. Bioinformatics. 2015 Oct 15;31(20):3377-9. doi: 10.1093/bioinformatics/btv372. PubMed PMID: 26099263; PubMed Central PMCID: PMCPMC4595900

Lorenz R, Bernhart SH, Honer Zu Siederdissen C, et al. ViennaRNA Package 2.0. Algorithms Mol Biol. 2011 Nov 24;6:26. doi: 10.1186/1748-7188-6-26. PubMed PMID: 22115189; PubMed Central PMCID: PMCPMC3319429

## Contact ##

Contact rsh46 -at- cam.ac.uk for bioinformatics related queries

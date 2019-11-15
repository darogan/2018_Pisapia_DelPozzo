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

HLA class II genes encode highly polymorphic heterodimeric proteins with the  function to present  antigens to T cells and to stimulate a specific immune response.  Many HLA genes are strongly associated to autoimmune diseases as they stimulate self-antigen specific CD4+  T cells that drives pathogenic response against proper tissues or organs. Recent evidences have demonstrated high expression of HLA class II risk genes associated to autoimmune disease  influencing  the strength of CD4+ T mediated autoimmune response.  The expression of HLA class II genes is regulated at both transcriptional and post-transcriptional levels. We have previously identified some proteins included in a RNP complex binding the 3’UTR and affecting mRNA processing . In this work,  we studied the regulation of HLA-DQ2.5 risk genes, the main susceptibility genetic factor for celiac disease (CD). The DQ2.5  molecule, encoded by HLA-DQA1*05 and HLA-DQB1*02 alleles, presents the antigenic gluten peptides to CD4+ T lymphocytes, activating the autoimmune response. Through molecular approaches we identified, the zinc-finger protein Tristetraprolin (TTP) or ZFP36, widely described as a factor modulating mRNA stability, among the components of the RNP complex. Since the 3’UTR of CD-associated HLA-DQA1*05 and HLA-DQB1*02 mRNA does not contain the canonical TTP binding consensus sequences, we performed a computational study focusing on mRNA secondary structure accessibility and stability. This <i>in silico</i> approach uncovered key structural differences specific to the CD-associated mRNAs allowing them to strongly interact with TTP through their 3’UTR, conferring a rapid turnover, in contrast to lower affinity binding to HLA non-CD associated mRNA.


### Figure 2 ###

Structures predicted using RNAVienna package (v2.4.12) (Lorenz et al, 2011) with the following command `RNAfold -T 37 < <FASTA FILE>`. The `-T 37` flag indicated the folding temperature should be 37&deg;C (i.e. human body temperature). Bracket notation secondary structure from RNAfold were rendered into images with selected motifs mapped onto them using FORNA (Kerpedjiev et al, 2015). Additional motifs were added from FOLDALIGN (Havgaard et al, 2007) and the `AURichness.pl` script as described below.

> Table: download sequences and structure predictions

| Allele | Sequence | RNAFold 2D Structures |
| -------| -------- | --------------------- |
| DQA101 | [[fasta](2DStructures/3DQA101.fasta)] | [[bracket](2DStructures/3DQA101.rnafold.bk)] [[stockholm](2DStructures/3DQA101.sto)] |
| DQA105 | [[fasta](2DStructures/3DQA105.fasta)] | [[bracket](2DStructures/3DQA105.rnafold.bk)] [[stockholm](2DStructures/3DQA105.sto)] |
| DQB102 | [[fasta](2DStructures/3DQB102.fasta)] | [[bracket](2DStructures/3DQB102.rnafold.bk)] [[stockholm](2DStructures/3DQB102.sto)] |
| DQB105 | [[fasta](2DStructures/3DQB105.fasta)] | [[bracket](2DStructures/3DQB105.rnafold.bk)] [[stockholm](2DStructures/3DQB105.sto)] |

> FOLDALIGN

All pairwise combinations of alleles were run through FOLDALIGN:

|  Comparison | Command | FOLDALIGN Results File |
| ----------- | ------- | ---------------------- |
| 3DQA101 vs 3DQA105 | `foldalign 3DQA101.fasta 3DQA105.fasta > 3DQA101_3DQA105.foldalign.txt` | [3DQA101_3DQA105.foldalign.txt](FOLDALIGN/3DQA101_3DQA105.foldalign.txt) |
| 3DQA101 vs 3DQB102 | `foldalign 3DQA101.fasta 3DQB102.fasta > 3DQA101_3DQB102.foldalign.txt` | [3DQA101_3DQB102.foldalign.txt](FOLDALIGN/3DQA101_3DQB102.foldalign.txt) |
| 3DQA101 vs 3DQB105 | `foldalign 3DQA101.fasta 3DQB105.fasta > 3DQA101_3DQB105.foldalign.txt` | [3DQA101_3DQB105.foldalign.txt](FOLDALIGN/3DQA101_3DQB105.foldalign.txt) |
| 3DQA105 vs 3DQB102 | `foldalign 3DQA105.fasta 3DQB102.fasta > 3DQA105_3DQB102.foldalign.txt` | [3DQA105_3DQB102.foldalign.txt](FOLDALIGN/3DQA105_3DQB102.foldalign.txt) |
| 3DQA105 vs 3DQB105 | `foldalign 3DQA105.fasta 3DQB105.fasta > 3DQA105_3DQB105.foldalign.txt` | [3DQA105_3DQB105.foldalign.txt](FOLDALIGN/3DQA105_3DQB105.foldalign.txt) |
| 3DQB102 vs 3DQB105 | `foldalign 3DQB102.fasta 3DQB105.fasta > 3DQB102_3DQB105.foldalign.txt` | [3DQB102_3DQB105.foldalign.txt](FOLDALIGN/3DQB102_3DQB105.foldalign.txt) |

### Figure 3 ###

To calculate the AU rich motifs in the riboprobes a stand alone Perl script takes the sequences and produces frequency tables. These tables are then inputted into the R script to produce the plots as seen in Figure 3.

The following commands run from the same directory will produce frequency tables and plots from them using R.

> Generate frequency tables:

    perl AURichness.pl

 > Example frequency table:

| Sequence  | Dinucleotide | Obs_freq           | Exp_freq           |
| ----------| ------------ | ------------------ | -------------------|
| DQB102    | AA           | 0.0144230769230769 | 0.0313408575810993 |
| DQB102    | AC           | 0.0673076923076923 | 0.0592935143426204 |
| DQB102    | AG           | 0.0576923076923077 | 0.0372702090153614 |

> Example density table:

| position | DQA101_di | DQA101_tri | DQA101_quad | DQA101_penta | DQA101_hexa | DQA105_di | DQA105_tri | DQA105_quad | DQA105_penta | DQA105_hexa
|----------|----------|----------|----------|----------|----------|----------|----------|----------|----------|----------|
| 0|0|0|0|0|0|0|0|0|0|0|
| 1|0|0|0|0|0|0|0|0|1|0|
| 2|0|0|0|0|0|0|0|0|0|0|

> Structure Accessibility

Calculated for each sequence and added to the AU frequency plots using RNAplfold (Lorenz et al, 2011) with a temperature of 37&deg;C, window size of 75 and mean structure score for fragments of 10 nt (`RNAplfold -T 37 -W 75 -u 10`).


> Plots

Create plots from frequency tables using the provided R script ([AURichAnalysis_V2.R](AURichAnalysis_V2.R)) (R v3.4.4):

    Rscript AURichAnalysis_V2.R

Example AU Rich Motif Observed Vs Expected Frequency plot

![AUMotif_ObsExp](AUMotif_ObsExp.png?raw=true=120x)

Examples AU Rich Motif Position plot

![AUMotif_Positions](AUMotif_Positions.png?raw=true=120x)



### Figure S1. Sfold Structure Comparison of 3’UTR of DQA1* and DQB1* ###

Sfold was run via the web server [sfold.wadsworth.org](http://sfold.wadsworth.org/cgi-bin/index.pl) (Ding et al, 2004). Additional motifs were added from FOLDALIGN (Havgaard et al, 2007)

### Figure S2. Sequence Alignment for human and mouse homologs for DQA1* and DQB1* alleles ###

Pairwise sequence alignments against the mouse sequences were performed using the Needlman-Wunch global alignment algorithm [EMBL-EBI API](http://www.ebi.ac.uk/Tools/psa/emboss_needle/) (Needleman & Wunsch, 1970; Madeira et al, 2019). Each DQ allele to mouse alignment was then manually appended and edited to optimise the alignment.

> Table: download sequences and alignments

| Alleles | Sequences | Alignments |
| ------- | --------- | ---------- |
| DQA     | [[fasta](SpeciesAlignments/DQA.human.mouse.fasta)] | [[clustal](SpeciesAlignments/DQA.human.mouse.faln)] |
| DQB     | [[fasta](SpeciesAlignments/DQB.human.mouse.fasta)] | [[clustal](SpeciesAlignments/DQB.human.mouse.faln)] |  


## References ##

Ding Y, Chan CY, Lawrence CE. Sfold web server for statistical folding and rational design of nucleic acids. Nucleic Acids Res. 2004 Jul 1;32(Web Server issue):W135-41. doi: 10.1093/nar/gkh449. PMID: 15215366; PMCID: PMC441587.

Havgaard JH, Torarinsson E, Gorodkin J. Fast pairwise structural RNA alignments by pruning of the dynamical programming matrix. PLoS Comput Biol. 2007 Oct;3(10):1896-908. doi: 10.1371/journal.pcbi.0030193. PubMed PMID: 17937495; PubMed Central PMCID: PMCPMC2014794

Kerpedjiev P, Hammer S, Hofacker IL. Forna (force-directed RNA): Simple and effective online RNA secondary structure diagrams. Bioinformatics. 2015 Oct 15;31(20):3377-9. doi: 10.1093/bioinformatics/btv372. PubMed PMID: 26099263; PubMed Central PMCID: PMCPMC4595900

Lorenz R, Bernhart SH, Honer Zu Siederdissen C, et al. ViennaRNA Package 2.0. Algorithms Mol Biol. 2011 Nov 24;6:26. doi: 10.1186/1748-7188-6-26. PubMed PMID: 22115189; PubMed Central PMCID: PMCPMC3319429

Madeira F, Park YM, Lee J, Buso N, Gur T, Madhusoodanan N, Basutkar P, Tivey ARN, Potter SC, Finn RD, Lopez R. The EMBL-EBI search and sequence analysis tools APIs in 2019. Nucleic Acids Res. 2019 Jul;47(W1) W636-W641. doi:10.1093/nar/gkz268. PMID: 30976793; PMCID: PMC6602479.

Needleman SB, Wunsch CD. A general method applicable to the search for similarities in the amino acid sequence of two proteins. J Mol Biol. 1970 Mar;48(3) 443-453. doi:10.1016/0022-2836(70)90057-4. PMID: 5420325.


## Contact ##

Contact rsh46 -at- cam.ac.uk for bioinformatics related queries

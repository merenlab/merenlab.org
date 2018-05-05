---
layout: post
title: "Assessing the completion of eukaryotic bins with anvi'o"
excerpt: "A preliminary set of 83 HMMs to assess the completion and redundancy of eukaryotic population genomes with anvi'o"
modified: 2018-05-05
categories: [anvio]
comments: true
authors: [tom]
---


{% include _toc.html %}

This post describes preliminary efforts to incorporate a collection of single copy core genes for __assessing the completion and redundancy of eukaryotic bins__ (stay tuned viruses, you are not forgotten and your time will come soon). 

##A brief background: 

Anvi'o v.4 contains two reference collections of single copy core genes for assessing the completion and redundancy of bacterial (Campbell et al.) and archaeal (Rinke et al.) bins using HMM models. __Currently missing are collections for viruses and eukaryotes__. This is a blind spot for the anvi'o users, and especially for those interested in [the genome-resolved metagenomics workflow](http://merenlab.org/2016/06/22/anvio-tutorial-v2/), where bins are determined in an interactive interface in which completion/redundancy values are estimated in real-time. Very handy to identify high-quality bacterial and archaeal genomes, not so much for the genomic content of eukaryotes and viruses...

Fortunately, the anvi'o developers have made it easy to [incorporate our own HMM models](http://merenlab.org/2016/05/21/archaeal-single-copy-genes/) to search for particular gene families of interest. Here I simply used this opportunity to test a __BUSCO__ collection of eukaryotic single copy core genes in anvi'o. This led to the creation of a promising (yet preliminary) __collection of 83 single copy core genes for eukaryotes "ready-to-be-used" from within anvi'o__.

The anvi'o compatible collection of 83 HMM models is available [here](https://www.dropbox.com/sh/z8d5kotx3tzkxnq/AADNXoNYKST1Pd5z7t_KZnOZa?dl=0).


##BUSCO: Benchmarking Universal Single-Copy Orthologs

Single copy core gene collections dedicated to different lineages of organisms have been created under the label __BUSCO__ for __B__enchmarking __U__niversal __S__ingle-__C__opy __O__rthologs. These collections are available from [here](http://busco.ezlab.org/), and those are the related articles:

BUSCO applications from quality assessments to gene prediction and phylogenomics. __doi: 10.1093/molbev/msx319__

BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs. __doi: 10.1093/bioinformatics/btv351__

In the following sections, I described the incorporation of a BUSCO collection of single copy core genes dedicated to __Protists__ into anvi'o and its subsequent testing using reference genomes and eukaryotic MAGs recovered from 3 distinct genome-resolved metagenomic projects. 

For the curious, one can define protists as all eukaryotic cells not affiliated with plants, animals or fungi. 

##The BUSCO collection of single copy core genes for Protists

The BUSCO collection of single copy core genes for Protists (named `protists_ensembl`) is available from [here](http://busco.ezlab.org/datasets/protists_ensembl.tar.gz). It contains HMM models for 215 single copy core genes.

To test the effectiveness of this collection when used from within anvi'o v.4, I generated CONTIGS databases for 4 complete eukaryotic genomes: _Bathycoccus prasinos_, _Micromonas pusilla_, _Ostreococcus tauri_ and _Pseudo-nitzschia multistriata_. 

I followed the tutorial to create our [own HMM models](http://merenlab.org/2016/05/21/archaeal-single-copy-genes/) for `protists_ensembl`, and tested it on the 4 reference eukaryotic genomes, using different e-value cut-offs. 

{:.notice} As far as I understand, BUSCO uses optimized length and e-value cut-offs for each HMM model. This flexibility is not available in anvi'o v.4, for which a single e-value cut-off can be defined for all HMM models within a collection. 

Some of the single copy core genes were detected dozens of times in all 4 reference eukaryotic genomes, even after varying e-value cut-off from E-15 to E-100. Others were detected once, as expected.

Until slight improvements are made on the side of the anvi'o developers (especially by defining optimized e-value and length cut-offs for each HMM model), I selected a __fixed e-value of E-25__. A total of 83 single copy core genes were detected __once in each of the 4 reference eukaryotic genomes__. I combined them in a new collection of single copy core genes dedicated to protists and named it __`sccgs_Protists_Preliminary`__.

As expected, here are the statistics when using this new collection:

```
Genome	Completion	Redundancy
Bathycoccus_prasinos	100.00%	0.00%
Micromonas_pusilla	100.00%	0.00%
Ostreococcus_tauri	100.00%	0.00%
Pseudo-nitzschia_multistriata	100.00%	0.00%
```
In the following section, __I tested the efficacy of this collection to estimate the completion and redundancy of eukaryotic MAGs__ characterized from different genome-resolved metagenomic projects. 

##Assessing the completion/redundancy of eukaryotic MAGs from 3 different projects

###TARA Oceans 

I recently spent some time [manually binning ~2.6 million scaffolds](http://merenlab.org/2017/05/21/thousand-genomes-from-tara/) assembled from metagenomes of the TARA Oceans project that cover the surface of four oceans and two seas. Among the [~1,000 MAGs we characterized in this project](https://www.biorxiv.org/content/early/2017/04/23/129791), 45 belong to the domain Eukarya. Until now, their completion had not been tested beyond the use of inappropriate bacterial single copy core genes. 

-__Micromonas__:

14 TARA Oceans MAGs were affiliated to _Micromonas_, and we know from cultivation that _Micromonas_ genomes are typically 20 Mbp long. This is a good case study to determine the accuracy of the new collection of single copy core genes.

Here are the related statistics:

```
Eukaryotic_MAG	Length	Genus	Completion	Redundancy
TARA_ANW_MAG_00088	19914044	Micromonas	66.27%	4.82%
TARA_ANE_MAG_00101	19009552	Micromonas	92.77%	4.82%
TARA_ASW_MAG_00046	18057117	Micromonas	91.57%	6.02%
TARA_ANW_MAG_00080	17700320	Micromonas	79.52%	1.20%
TARA_RED_MAG_00119	15250294	Micromonas	84.34%	3.61%
TARA_ASE_MAG_00032	12391309	Micromonas	66.27%	4.82%
TARA_ANW_MAG_00074	9401388	Micromonas	53.01%	1.20%
TARA_ANE_MAG_00090	5559939	Micromonas	24.10%	1.20%
TARA_ASW_MAG_00039	3501361	Micromonas	10.84%	0.00%
TARA_ASW_MAG_00034	3334468	Micromonas	22.89%	1.20%
TARA_ANW_MAG_00075	2699399	Micromonas	20.48%	1.20%
TARA_ANW_MAG_00081	2279971	Micromonas	10.84%	0.00%
TARA_ASW_MAG_00037	2179362	Micromonas	15.66%	0.00%
TARA_ANE_MAG_00099	2133633	Micromonas	8.43%	0.00%
```

__There is a strong correlation between length and completion estimates (R-scare > 0.9)__, and MAGs with a length close to 20 Mbp have high completion values (with the exception of TARA_ANW_MAG_00088 for which completion seems a little off - but this can be due to assembly or binning problems).


-__Ostreococcus__:

7 TARA Oceans MAGs were affiliated to _Ostreococcus_ (reference genomes are ~13 Mbp), and here are the related statistics:

```
Eukaryotic_MAG	Length	Genus	Completion	Redundancy
TARA_PSW_MAG_00136	11885609	Ostreococcus	93.98%	1.20%
TARA_ASE_MAG_00036	11537750	Ostreococcus	93.98%	2.41%
TARA_RED_MAG_00118	10110290	Ostreococcus	83.13%	1.20%
TARA_ANE_MAG_00093	9555746	Ostreococcus	78.31%	2.41%
TARA_PON_MAG_00082	8114168	Ostreococcus	68.67%	2.41%
TARA_PSE_MAG_00129	3144407	Ostreococcus	28.92%	0.00%
TARA_RED_MAG_00116	2042466	Ostreococcus	14.46%	0.00%
```

This time the correlation between length and completion was even stronger, with R-scare > 0.995.

 ![Alt Image Text](Figure_1.png "Figure 1")
The figure describes the linkage between genomic length and completion estimates for _Micromonas_ and _Ostreococcus_ MAGs characterized from the surface ocean, using `sccgs_Protists_Preliminary` in anvi'o v.4.
 

Overall, __average redundancy estimates for the 45 eukaryotic MAGs went from ~20% with the bacterial single copy core gene collection to 1.4% with `sccgs_Protists_Preliminary`__. This suggests the eukaryotic MAGs represent individual population genomes with no substantial amounts of contamination. As a result, it seems key information computed by anvi'o (especially differential coverage and sequence composition) is also effective for the characterization of (sometimes near-complete) eukaryotic genomes. Good.

###Pseudo-nitzschia from the Southern Ocean

In another project led by Anton Post, we have characterized using genome-resolved metagenomics a 34 Mbp MAG corresponding to _Pseudo-nitzschia_. This lineage was prevailing in the Ross Sea polynya in the coast of Antarctica during the 2013-2014 austral summer. The story is not yet published, and until now I had never estimated the completion or redundancy of this MAG. 

Based on `sccgs_Protists_Preliminary`, this MAG is 86.8% complete and 2.4% redundant. This made my day...

###Candida albicans in a gut metagenome 

[Sharon et al. (2013)](https://www.ncbi.nlm.nih.gov/pubmed/22936250) characterized a MAG corresponding to _Candida albicans_ (a fungi) from metagenomes of an infant gut, and we later recapitulated this finding using [our own metagenomic binning effort of the same dataset](https://peerj.com/articles/1319/). The _Candida albicans_ MAG we recovered is 13.8 Mbp long. Based on `sccgs_Protists_Preliminary`, this MAG is 83.1% complete and 7.2% redundant. A reference [Candida albicans genome](https://www.ncbi.nlm.nih.gov/genome/?term=Candida%20albicans) is 14.7 Mbp long, suggesting the estimates are not aberrant.  

Not to bad given that fungi are not protists. This suggests such collection could work to a certain extent beyond the realm of protists. However, much more testing is needed before a stable collection (or multiple lineage-specific collections?) of eukaryotic single copy core genes is incorporated into a new anvi'o version release.

Until then, one can download the "ready-to-use" folder required to run HMM models from the `sccgs_Protists_Preliminary` collection using [this link](https://www.dropbox.com/sh/z8d5kotx3tzkxnq/AADNXoNYKST1Pd5z7t_KZnOZa?dl=0), and store it in a ~PATH.

It is then possible to follow this simple workflow:

```
anvi-gen-contigs-database -f FASTA.file -o CONTIGS.db
anvi-run-hmms -c CONTIGS.db -H ~PATH/sccgs_Protists_Preliminary/
anvi-compute-completeness -c CONTIGS.db --completeness-source sccgs_Protists_Preliminary
```
which will print in your terminal the completion and redundancy estimates of your FASTA file of interest using the collection described in this post.

###Final notes 

It is now possible to estimate the completion and redundancy levels of genomes from the 3 domains of life with anvi'o. 

However, please keep in mind that single copy core gene collections only provide a rough estimation of the completeness and contamination level of MAGs. __Visualizing a MAG in the context of recruited reads from metagenomes for instance (the environmental signal) remains a key asset when it comes to genome-resolved metagenomics__, and this approach very often allows the detection (and subsequent removal) of suspicious contigs that contain no single copy core genes and hence are invisible from the perceptive of completion / redundancy estimates. Thus, these estimates are only one of multiple parameters one can (and maybe should) use to characterize, manually refine, and assess the biological relevance of eukaryotic MAGs.

Of course, these problems are less relevant if one is working with cultivar genomes and single cell genomes, or are they?

I leave you to this though :)

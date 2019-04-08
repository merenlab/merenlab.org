---
layout: page
title: Vocabulary
modified: 2019-04-05
excerpt: "Commonly used terms and their approximate meanings"
comments: true
---

This purpose of this document is to create a resource that defines commonly used terms in microbial 'omics.

We initially wanted to create a resource just for the terms we commonly use in anvi'o, especially those that are specific to the platform. However, we concluded that including definitions for a larger set of terms (including those that are common in 'omics studies) would be the only way to have a cohesive document.

Suggestions and edits are most welcome.

{% include _join-anvio-slack.html %}

{% include _toc.html %}

{:.notice}
{% include _fixthispage.html source="vocabulary.md" %}

## All things 'omics

This section defines terms that are commonly used in 'omics studies.

### Sequencing reads

Commonly used to define the raw output of a sequencer. These are strings of the alphabet {A, C, T, G}, representing nucleotide sequences of the input DNA or RNA.

### Reference sequence

A sequence that you know something about. In the context of metagenomics, it often refers to a sequence that one uses for read recruitment. 

### Read recruitment

A set of computational strategies to align sequencing reads to one or more reference sequences. Also known as 'mapping'.

In the context of metagenomics, read recruitment allows one to estimate the coverage of a given reference sequence in a given metagenome by identifying all short reads that match to it. Understanding this strategy, along with its power and caveats, is one of the most important steps to fully appreciate most ‘omics strategies and the ways they lend themselves to study the ecology and evolution of microbial populations.

### Coverage

Average number of sequencing reads that map to each nucleotide posotion in a reference. Also known as 'depth of coverage'.

### Contig

A contiguous segment of DNA that is often 'assembled' from short sequencing reads, but still represents only a fraction of the longer context to which it belongs.

Assembly software often takes short metagenomic reads and yields a list of contigs. Although their lengths vary depending on a myriad of factors, contigs are often orders of magnitude longer than short reads that were assembled, which makes them suitable for downstream efforts that may include gene calling, functional and/or taxonomic annotation, and/or metagenomic binning. Contigs often become reference sequences for read recruitment analyses.

### Metagenomic binning

A set of computational strategies that aims to identify and put together contigs that belong to the same population. These strategies often use differential coverage of contigs (when multiple samples are present) and/or sequence composition information (such as tetra-nucleotide frequency).

### Differential coverage

Change in coverage of a reference sequence across multiple samples. This statistic is one of the essential information for most binning algorithms.

### Tetra-nucleotide frequency

The ratio of all 4-nucleotide words in a given contig. The tetra-nucleotide frequency is largely preserved throughout microbial genomes, which enables the identification of distinct contigs that likely originate from the same population.

### Metagenome-assembled genome (MAG)

A genome bin that meets certain quality requirements and can be assumed to reprsent contigs from one bin of a metagenome, which collectively represent the DNA of (what we think is) a single population.

### Population

Although frequently used, microbiology does not have a precise and consensus definition for what is a population, and how to define boundaries of enviornmental populations. One of the operational definitions our group often uses suggests that a population is an assemblage of co-existing microbes in an environment whose genomes are similar enough to map to the context of the same reference genome.

### Metagenome

The entire DNA content of an environment. Since most environments harbor many different organisms, the metagenome includes genetic information from a large collection of genomes. High-throughput sequencing of metagenomes produce tremendous amount of sequencing reads that can be used for assembly or read recruitment.

<div class="extra-info" markdown="1">

<span class="extra-info-header">Iva Veseli on all things 'meta'</span>

In general, if you see the prefix 'meta' in front of an ‘omics term, it means that you are extending that term to apply to all the things within a given sample. For example, the metatranscriptome is the set of transcriptomes from every population in your sample. The metametabolome are all the metabolites coming from all organisms in an environment and so on. In case that's hard to remember, you can think of the 'All the Things' meme.

![](../images/allthethings.jpg){:.center-img}

Works for us. :)
</div>

### Metagenomics

The study of environmental metagenomes.

### Pangenome

From a computational standpoint, the term pangenome broadly refers to the entire collection of genes found in two or more genomes.

### Pangenomics

The family of computational strategies that determine the pangenome and make it accessible as a framework to study relationships between a set of genomes through gene clusters.

### Gene cluster

Fundamental units of pangenomes which appear in the literature also as ‘protein clusters’, ‘orthogroups’, ‘groups of orthologous genes’, or ‘operational protein families’ (and they should not be confused with biosynthetic gene clusters which describe functionally related genes that belong to the same operon in a single chromosome).

Commonly used computational strategies for pangenomics that consider entire contents of input genomes determine gene clusters typically by (1) identifying all genes among a set of genomes, (2) computing similarities between each gene using translated DNA sequences, and (3) determining which genes are homologous enough to be described in the same cluster. Hence, a gene cluster in a given pangenome corresponds to a de novo identified virtual construct that contain one or more genes from one or more genomes.

For a tutorial on how to do pangenomics with anvi’o, see [this page](http://merenlab.org/2016/11/08/pangenomics-v2/).


<div class="extra-info" markdown="1">

<span class="extra-info-header">Iva Veseli on all things 'pan'</span>

Where 'meta' indicates that we apply our analysis across all populations in a given sample, 'pan' indicates that we apply our analysis within a single group. We couldn’t find a good meme for 'pan', so please accept this analogy instead: Imagine you have a frying pan, in which you are making an omelet. To get the nutritional information for that omelet, you analyze what is in the pan. Some ingredients will contain the same nutrients, like protein or carbohydrates. Other nutrients, like specific vitamins, will be exclusive to certain ingredients. So the 'pan'-nutrients will include all the nutrients within the frying pan, regardless of which ingredient they came from (see what we did there?). 

</div>

### Single-nucleotide variant (SNV)

A nucleotide position where the identity of all bases mapping to this position varies (beyond the expected rate of sequencing error).

SNVs are characterized by (1) the position in the reference sequence where the difference occurs, and (2) a frequency vector that quantifies the frequency of nucleotide identities that mapped onto that position. See [this page](http://merenlab.org/2015/07/20/analyzing-variability/#single-nucleotide-variants) for a lengthier discussion of SNVs. Also see the definition of single-codon variant and single-amino acid variant.

## All things anvi'o

This section defines terms that are largely specific to anvi'o.

### Interactive interface

Highly customizable visualization environment in anvi'o to interact with data. It is integrated into most of the core functionalities of anvi’o, but can also be used with tabular data inputs. See [this document](http://merenlab.org/2016/02/27/the-anvio-interactive-interface/) for more details.

### Split

A fragment of a contig in your anvi'o analysis. Anvi'o, by default, breaks up long contigs you provide into multiple parts, called 'splits', purely for visualization purposes so the researcher can distinguish a very long contig from shorter once. This way, a single contig that represents 1,000,000 base pairs would not have less visual significance than 10 contigs each of which are about 10,000 base pairs. The default split size is 20,000 bp, but you can change this using the `--split-length` parameter when you generate an anvi'o contigs database. 

### Layer

Every concentric circle in anvi'o interactive interfaces (in radial display mode). Visualization of a minimal metagenome in anvi'o intearctive interface layers will include the parent, GC-content, and metagenomic samples.

### Item

Data points shown in layers. Items can be a lot of things in anvi'o: they will be splits in metagenomic mode, genes in gene mode, gene clusters in pangenome mode, or genome bins in collections mode.

### Parent layer

The special layer in anvi'o interactive interfaces that describe which splits belong to which contigs (if your contigs were long enough to be split).

### Items organization

The center piece of an anvi'o interactive display that organizes items. It could be a hierarchical clustering dendrogram based on an anvi'o clustering configuration, or a user-provided phylogenetic tree. These organizations can also utilize alphabetical orders, or additional user-provided layers.

### Clustering configuration

An advanced description of how you want anvi’o to cluster your contigs. A clustering configuration is described in a text file which specifies things like which data sources you want to use for clustering and how to normalize each data source.

### View

Statistic behind the data points shown for a given item in a given anvi'o display. A particular view considers some statistic associated with your items that anvi'o calculated automatically, or provided externally, to produce an informative (and hopefully beautiful) graph. For more information about the types of views in anvi’o, please take a look at [this post](http://merenlab.org/2017/05/08/anvio-views/) by Mike Lee.

### Contigs database

A self-contained database containing a lot of information associated with your contig (or scaffold) sequences. This includes data that isn’t dependent on which sample the contigs came from, like positions of open reading frames, k-mer frequencies, split start/end points, functional and taxonomic annotations among others. You can initialize a basic contigs database from a FASTA file with the command `anvi-gen-contigs-database`, and supplement it with additional information later in your analysis.

### Profile database

A database containing sample-specific information about your contigs; for instance, coverage information from mapping reads to the contigs in a sample. Single profiles, each of which contains data for a particular sample, can be combined into a merged profile if they link to the same contigs database. The information across samples in a merged profile can be visualized as a ‘view’ in the anvi’o interactive database.

### HMM profile

Information about known genes that can be used to search for the presence of these genes (dubbed 'hits') in contigs, using Hidden Markov Models. More formally, these profiles probabilistically describe variable and conserved regions in a set of homologous sequences. Anvi’o ships with four HMM profiles of bacterial single-copy core genes, but you can also use your own custom HMM profile if you so desire.


### Variability profile

Information on residue variants across samples. These variants, wehther they are SNVs, SCVs, or SAAVs, are usually identified for each sample when a profile database is generated for the sample. This information can later be combined across samples into a variability profile using the command `anvi-gen-variability-profile`. For an extensive tutorial on how to analyze variability profiles using anvio plese refer to [this resource](http://merenlab.org/2015/07/20/analyzing-variability/).

### External gene calls file

A text file containing gene predictions that you obtained from external software, but wish to import into anvi’o when generating your contigs database. This file must follow the format specified [here](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-gen-contigs-database).

### External functions file

A text file containing functional predictions that you obtained from external software, but wish to import into anvi’o using `anvi-import-functions`. This file must follow the format specified [here](http://merenlab.org/2016/06/18/importing-functions/).

### External taxonomy file

A text file containing taxonomy information for your genes that you obtained from external software, but wish to import into anvi’o using `anvi-import-taxonomy-for-genes`. This file must follow the format specified [here](http://merenlab.org/2016/06/18/importing-taxonomy/#simple-matrix).

### Collection

A virtual construct to store bins of items in an anvi'o profile database. Each collection contains one or more bins, and each bin contains one or more items. These items can be gene clusters, contigs, or other things depending on the display mode.

### Single-codon variant (SCV)

The lazy definition - an SCV is just like a SNV, but for a codon position. The complete definition - A codon position where the entire, 3-base codon identity (there are 64 possible codons) is different between a reference coding region and mapped reads. SCVs are characterized by 1) the position in the reference sequence and 2) the frequency of codons in that position. If you want an even longer explanation, see [this page](http://merenlab.org/2015/07/20/analyzing-variability/#single-codon-variants).

### Single-amino acid variant (SAAV)

Differences between mapped reads at the amino acid residue level. See [this page](http://merenlab.org/2015/07/20/analyzing-variability/#single-amino-acid-variants) for a more comprehensive description.


{:.notice}
{% include _fixthispage.html source="vocabulary.md" %}

{% include _join-anvio-slack.html %}

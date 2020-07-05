---
layout: post
authors: [meren,tom,evan]
title: "Studying microbial population genetics with anvi'o"
excerpt: "Exploring micro-diversity patterns using for deeper insights into ecology"
modified: 2016-07-21
tags: []
categories: [anvio]
comments: true
redirect_from: /mpg/
---

{:.notice}
This is a theoretical tutorial describing how to characterize SNVs, SCVs and SAAVs with anvi'o, and
how to interpret the output. For a more practical tutorial on the same topic, please visit
[the infant gut tutorial]({{ site.url
}}/tutorials/infant-gut/#profiling-snvs-in-a-bin) for a tutorial on profiling SNVs, SCVs, and SAAVs,
and check out
[the anvi-3dev blog post]({{ site.url }}/2018/09/04/getting-started-with-anvi-3dev/)
for a tutorial on visualizing SCVs and SAAVs directly on the protein structure they encode for. Also
check out our [reproducible workflow on oceanic SAR11]({{ site.url
}}/data/2018_Delmont_and_Kiefl_et_al_SAR11_SAAVs/), which leverages the concepts in this tutorial to
gain insight into ecologically-linked patterns of micro-diversity on a global scale.

{:.notice}
We thank [Rika Anderson](https://twitter.com/RikaEAnderson) for carefully reading this tutorial, and embarrassing us by fixing our mistakes multiple times.

{% include _toc.html %}

We use the term "**population**" to describe an assemblage of co-existing microbial genomes in an environment that are similar enough to map to the context of the same reference genome.

Cells within a population (i.e., all cells that would have classified as the same "species" or "strain" for whatever these terms mean to you) will share the vast majority of their genomes in the sequence space. Hence, a consensus contig obtained through the assembly of metagenomic reads from this environment, or a contig that originates from an isolate cultured from this environment, will recruit short reads through mapping from many other member cells in a given population.

Metagenomic short reads from a cloud mapping to the same context can have nucleotide positions that systematically differ from the reference context depending on, 

* The heterogeneity of the population,
* the fraction of the population that can be targeted by the reference genome through mapping, 
* and/or the stringency of mapping.

Since one of the mechanisms for the diversification of genomes operates at the single-nucleotide level, ability to characterize the heterogeneity of a population at the nucleotide, codon, and amino acid levels can be very critical for a more complete understanding of the environmental forces that affect communities, and adaptive strategies microbes rely on to survive in environments they reside. Anvi'o offers *de novo* strategies to characterize and visualize single nucleotide variants (SNVs), single codon variants (SCVs), and single amino acid variants (SAAVs) in a given population using mapping results, and allows researchers to delineate subtle ecological niche and ecotypes, and quantify the level of heterogeneity in a population.

---

This article will describe some key aspects of the anvi'o workflow for the recovery, profiling, and characterization of SNVs, SCVs, and SAAVs for high-resolution genomics using easy-to-use anvi'o programs. It also gives a comprehensive description of each parameter of the variability table outputs.

# An intro to single nucleotide/codon/amino acid variation

Depending on the complexity of a microbial cloud (i.e., the extent of monoclonality in it, or the lack thereof), metagenomic short reads that map to a reference context (i.e., a contig from a MAG or a cultivar genome) can generate one or more mismatches (unless the mapping software requires a 100% sequence identity for alignment).

Here we define 'variation' as the extent of disagreement between aligned sequences that map to the same context. While the source of a mismatch can be artificial (i.e., due to random sequencing or PCR errors, or non-specific mapping of local alignments), it may also represent ecologically informative variation.

## Single nucleotide variants

Simply put, single nucleotide variants (SNVs) are nucleotide positions where there is disagreement between the identity of mapped reads and the  the reference context. A SNV can be fully characterized by both _1)_ its position in the reference sequence, and _2)_, a frequency vector that quantifies the frequency of nucleotide identities that mapped onto that position. For example, a SNV where the where the mapped read composition is 70% `A` and 30% `C` has a frequency vector of `<A=0.7, C=0.3, G=0.0, T=0.0>`.

SNV positions are a minority compared to stable frequencies, where the extent of disagreement can be attributable to the rate of sequencing error or other nonbiological effects. There exists no authoritative cut off that states what level of disagreement is required for position to be considered a SNV, however below we discuss some things to bear in mind. Note that SNVs are generally dominated by two nucleotides (e.g., 60% of `A` and 40% of `T`) but can also represent three or four nucleotides in relatively high proportions. Also note that a nucleotide position may be a SNV in one metagenome, but not in another.


## Single codon variants

Single codon variants (SCVs) are just like single nucleotide variants (SNVs), however they focus on explaining the frequencies of *codons* instead of *nucleotides*.

A SCV can be uniquely and fully characterized by both (1) its unique position in the reference sequence, and (2) the _frequency vector_ that quantifies the frequency of codons that mapped onto that codon position. The frequency vector for a SCV has a length (as there are 64 codons). 

Why looking at SCVs instead of the canonical SNVs when there are 64 codons (instead of 4 nucleotides) and SCVs can only recovered from coding regions of the reference (instead of nucleotides, which can convey information about environmental variability independent of coding regions)? The primary reason for using SCVs is to be able to accurately determine whether an environmental variant is synonymous or non-synonymous. SNVs can only give reliable access to the corresponding codon variant *if and only if* number of SNVs in a given codon position does not exceed one (which is often the case for eukaryotic genomes such as the human genome, studies of which gave rise to lots of approaches for population genetics (with rather limited applicability to microbial world)). Please consider the following example:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2015-07-20-analyzing-variability/SCV_demo1.png"><img src="{{ site.url }}/images/anvio/2015-07-20-analyzing-variability/SCV_demo1.png" style="border: none; width: 100%;" /></a>
</div>

The top row is the reference sequence that contains the codon `TCC`, which encodes the amino acid serine. Each row below shows the reads that mapped at least partially to this codon. the first position in the codon is a SNV (let's call it SNV 1) and we can characterize this SNV with a frequency vector shown in red, i.e. 37% of reads are `A` and 63% are `T`. If we ignore the other positions **this SNV leads to the variants, `TTC` (serine) and `ATC` (isoleucine)**. Likewise, the second position in the codon is also a SNV (SNV 2) with the frequency vector shown in red:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2015-07-20-analyzing-variability/SCV_demo3.png"><img src="{{ site.url }}/images/anvio/2015-07-20-analyzing-variability/SCV_demo2.png" style="border: none; width: 100%;" /></a>
</div>

If we ignore again the other positions, **this SNV leads to the variants,`TCC` (serine) and `TGC` (cysteine)**.

We can not combine the independent nucleotide variants when we are working with such nucleotide pileups to resolve the true codon variants. Our implementation of SCVs addresses this issue by **working only with short reads that fully cover a given codon position**, so the linkage between nucleotides is no longer hidden:
 
 <div class="centerimg">
<a href="{{ site.url }}/images/anvio/2015-07-20-analyzing-variability/SCV_demo3.png"><img src="{{ site.url }}/images/anvio/2015-07-20-analyzing-variability/SCV_demo3.png" style="border: none; width: 100%;" /></a>
</div>

{:.notice}
Indeed this approach reduces the number of short reads we can work with. For instance, in our toy
example 4 (33%) of the reads suffer this fate. But in a hypothetical scenario where [variation in
coverage in recruitment results](http://merenlab.org/2016/12/14/coverage-variation/) is
non-existent, the expected fraction of those to be excluded should be equivalent to \$4 / (2 + L)\$
where \$L\$ is the length of your short reads.

The frequencies of remaining codons can now be read directly, which would yield the frequency vector shown in red. The resulting SCV is characterized as 50% `TCC` and 50% `AGC`, which are both codons of serine, and thus the **SCV is 100% synonymous**. The use of SNVs rather than SCVs would have overestimated non-synonymous variants.

In summary, the use of SCVs will avoid over- or under-estimation of non-synonymous variants in metagenomic read recruitment studies where the occurrence of multiple SNVs within single codons is likely.

## Single amino acid variants

Maybe you're only interested in positions that actually change the structure of the protein encoded by the gene. In such cases, single amino acid variants (SAAVs) are for you. They are just like SCVs, however they focus on explaining the frequencies of encoded amino acids instead of codons. A SAAV can be fully characterized by both _1)_ its position in the reference sequence, and _2)_, a _frequency vector_ that quantifies the frequency of amino acid identities that mapped onto that codon position. (Here when I say amino acid identities, I really referring to the amino acid encoded by the codon that mapped). The frequency vector for a SCV has a length of 20+1 because there are 20 amino acids and 1 stop codon. If you use SAAVs, do your due diligence and make sure the genes you're profiling are actually translated into proteins, otherwise you are dealing with nonsense.

## Summary bullet points

1. Use SNVs if you are interested in single nucleotide positions or non-coding regions

2. Use SCVs if you want to resolve codon variants or related concepts like synonymity

3. Use SAAVs if you are interested in variation that likely leads to structural differences in the encoded protein.


# The anvi'o way

Although the importance of it is probably almost common sense to anyone who understands the steps we follow to recover our contigs for a given environment, and although many amazing studies exist including the ones from Jill Banfield's lab that made use of SNV patterns in metagenomic settings such as [this one](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0060177), working with subtle variation in cultivar genomes or metagenomic assembly outputs hasn't been so practical to access. There are **many** challenges to address if one wants to focus on subtle things in complex environmental sequencing datasets, and in order to set the stage with reasonable expectations we would like to clarify that anvi'o will not exactly offer you a magical wand. However, anvi'o _will_ empower you to work with this aspect of your datasets by making the recovery and query of this information more or less straightforward. We hope that at the end of this post it will be clear to you how this can be done.

## A taste of how anvi'o can visualize variation

Just so you have a mental picture of how anvi'o visualizes SNVs, here is a polished screenshot from the interface:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2015-07-20-analyzing-variability/example.png"><img src="{{ site.url }}/images/anvio/2015-07-20-analyzing-variability/example.png" style="border: none; width: 100%;" /></a>
</div>

What you see is the coverage and base frequencies in reported SNVs in a single contig across two samples. Where these samples are coming from are not really relevant to the rest of this post. But usually, just like it is the case in this example, SNVs and their base frequencies do not seem to be random, and they can be extremely reproducible if the population structures are really similar in multiple samples (i.e., in biological replicates).

{:.notice}
Do you want to visualize variants directly onto protein structure? Good news, because anvi'o does that! Check out our tutorial on the anvi'o structure database for details on how ([click me](http://merenlab.org/2018/09/04/structural-biology-with-anvio)).


## How to generate variability tables

Anvi'o gives access to information about variants in two steps.

The first step is to identify variants and reporting them _for each sample separately_. This step takes place during the profiling of a given BAM file with the program `anvi-profile`. By default, SNVs are profiled during `anvi-profile`, but SAAVs and SCVs are not. This is because it takes much longer to profile. To work with SAAVs and SCVs, ensure that your `anvi-profile` command has the flag `--profile-SCVs`. There is no `--profile-SAAVs` flag because SAAVs are generated from SCVs. 

The second step is to interpret the ecological significance of sample-specific variants _across samples_ using the program `anvi-gen-variability-profile`. The first step is agnostic to the experimental design, and/or the variants in other samples: it simply does its best to find and report variation. The second step is where the user's input about the experimental design, and their questions come into the play.

The following two chapters in this article will detail these two steps, and the third chapter will explain the analysis of variability using the Infant Gut Dataset that is first appeared in [this study](http://genome.cshlp.org/content/23/1/111.long).

Everything you might need to replicate this analysis and to re-generate the third figure in the anvi'o paper is [here]({{ site.url }}/data/). We hope that we will manage to give enough details that you will be able to modify the recipe and start exploring patterns in your own datasets if you wish.

Please do not hesitate to ask questions, make suggestions, and/or start discussions regarding this topic here in the comments section down below, or at [the issues section of the anvi'o repository on Github](https://github.com/meren/anvio/issues), or by directly contacting [us]({{ site.url }}/people/).

# *De novo* characterization and reporting of SNVs

During the profiling step that should be done for each sample separately, anvi'o checks the
composition of nucleotides that map to each reference position, and keeps track of variation in the
form of 'base frequencies'.

Note that the reporting of these base frequencies is a little tricky: If you report base frequencies for each position when there is any variation at all, then the number of reported positions would be enormous. For instance, if you have a nucleotide position with 500X coverage, you can safely assume that there will always be _some_ nucleotides that do not match to the consensus nucleotide for that position--whether it is due to biology, or random sequencing errors, or other types of relevant or irrelevant sources of variation. For the sake of better management of available resources, and to be very quick, anvi'o (with its default settings) ___does not___ report variation in every single nucleotide position during profiling.

Instead, it relies on the following conservative heuristic to identify SNVs and report the variation only at these nucleotide positions:

$$ y = \left(\frac{1}{b}\right)^{x^{\frac{1}{b}} - m} + c $$

where, \$x\$ represents the coverage, and \$b\$, \$m\$, and \$c\$ represent the model parameters
equal to 3, 1.45, and 0.05, respectively. Assuming \$n_1\$ and \$n_2\$ represent the frequency of
the most frequent and the second most frequent bases in a given nucleotide position (see the table),
base frequencies are reported only if \$n_2/n_1 > y\$ criterion is satisfied for a given coverage of
\$x\$. It is this simple (and ugly, in a sense). But briefly, this approach sets a dynamic baseline
for minimum variation required for reporting _as a function of coverage depth_. According to this
heuristic, \$y\$ would be 0.29 for 20X coverage, 0.13 for 50X coverage, 0.08 for 100X coverage, and
~0.05 for very large values of coverage as \$y\$ approaches to \$c\$. The goal here is to lessen the
impact of sequencing and mapping errors in reported frequencies, and it does it's boring job.

This heuristic is affirmatively _not_ the be all end all, and it is an area or our workflow we could
stand to improve. We have chosen these parameters conservatively such that it reports SNVs that one
can really trust. However, we are definitely culling out some degree of biologically relevant
variation. If you're interested in learning more you should check out the supplemental material of
[this study](https://www.nature.com/articles/nature24287), as we think they are the people who are
thinking the best thoughts about how to identify SNVs in metagenomic contexts.

This computation- and storage-efficient strategy reports a relatively short list (could still be
millions of SNVs) of SNVs. However, the user always has the option to instruct the profiler to store
_all observed frequencies_ by declaring this intention with `--report-variability-full` flag for
more statistically appropriate downstream analyses. We are of course talking about the type of
analyses Christopher Quince-like people would like to do.

# Generating a SNV, SCV, or SAAV profile

Once SNVs are stored in a single or merged anvi'o profile database, it is time to decide how to process that information, and export a scrutinized table to make sense of them. To interpret the ecological significance of sample-specific variable positions across samples, anvi'o installs a helper program called `anvi-gen-variability-profile`.

In this section will first describe the output file structure for SNV profiles, and then describe the parameters of `anvi-gen-variability-profile`.

## The output matrix

The output generated by `anvi-gen-variability-profile` annotates each SNV with about 30 columns of information. Each line in the output file represents one SNV (reminder: a variant nucleotide position in a metagenome), and each column represents various metrics, classifications, and identifiers for that SNV. Here is the first 10 lines of an example SNV profile for the *E. facealis* bin in the Infant Gut Dataset so you have an idea of the file format, followed by the description of each column:

entry_id  |  unique_pos_identifier  |  pos    |  pos_in_contig  |  sample_id  |  corresponding_gene_call  |  in_partial_gene_call  |  in_complete_gene_call  |  base_pos_in_codon  |  codon_order_in_gene  |  codon_number  |  gene_length  |  gene_coverage       |  non_outlier_gene_coverage  |  non_outlier_gene_coverage_std  |  coverage  |  mean_normalized_coverage  |  cov_outlier_in_split  |  cov_outlier_in_contig  |  A     |  C     |  G     |  N     |  T     |  reference  |  consensus  |  competing_nts  |  departure_from_reference  |  departure_from_consensus  |  n2n1ratio             |  entropy
----------|-------------------------|---------|-----------------|-------------|---------------------------|------------------------|-------------------------|---------------------|-----------------------|----------------|---------------|----------------------|-----------------------------|---------------------------------|------------|----------------------------|------------------------|-------------------------|--------|--------|--------|--------|--------|-------------|-------------|-----------------|----------------------------|----------------------------|------------------------|---------------------
1875      |  1413                   |  654    |  654            |  DAY_17A    |  1                        |  0                     |  1                      |  3                  |  146                  |  147           |  1005         |  102.06169154228856  |  101.47290640394088         |  14.584491558359868             |  106       |  1.0385875287602855        |  0                     |  0                      |  100   |  6     |  0     |  0     |  0     |  A          |  A          |  AC             |  0.05660377358490566       |  0.05660377358490566       |  0.06                  |  0.217518571336808
1872      |  1410                   |  2562   |  2562           |  DAY_17A    |  3                        |  0                     |  1                      |  2                  |  0                    |  1             |  891          |  97.4769921436588    |  99.25173210161662          |  24.360720721889436             |  66        |  0.6770828535900152        |  1                     |  0                      |  0     |  62    |  4     |  0     |  0     |  T          |  C          |  CG             |  0.06060606060606061       |  0.06060606060606061       |  0.06451612903225806   |  0.22863187358286127
1876      |  1414                   |  2941   |  2941           |  DAY_17A    |  3                        |  0                     |  1                      |  3                  |  126                  |  127           |  891          |  97.4769921436588    |  99.25173210161662          |  24.360720721889436             |  74        |  0.7591535025100171        |  1                     |  0                      |  0     |  68    |  4     |  0     |  2     |  T          |  C          |  CG             |  0.08108108108108109       |  0.08108108108108109       |  0.058823529411764705  |  0.3330110964801871
1886      |  1420                   |  3005   |  3005           |  DAY_17A    |  3                        |  0                     |  1                      |  1                  |  148                  |  149           |  891          |  97.4769921436588    |  99.25173210161662          |  24.360720721889436             |  58        |  0.5950122046700134        |  1                     |  0                      |  2     |  2     |  54    |  0     |  0     |  C          |  G          |  AG             |  0.06896551724137931       |  0.06896551724137931       |  0.037037037037037035  |  0.2987580581893401
8560      |  5890                   |  6074   |  6074           |  DAY_16     |  8                        |  0                     |  1                      |  2                  |  21                   |  22            |  957          |  450.43782654127483  |  471.07912844036696         |  76.15909786907412              |  98        |  0.21756609730692766       |  1                     |  1                      |  0     |  6     |  0     |  0     |  92    |  G          |  T          |  CT             |  0.061224489795918366      |  0.061224489795918366      |  0.06521739130434782   |  0.23032354087587759
8561      |  5890                   |  6074   |  6074           |  DAY_22A    |  8                        |  0                     |  1                      |  2                  |  21                   |  22            |  957          |  488.79205851619645  |  513.3886925795053          |  70.26091145846179              |  121       |  0.24754903008717885       |  1                     |  1                      |  0     |  7     |  0     |  0     |  114   |  G          |  T          |  CT             |  0.05785123966942149       |  0.05785123966942149       |  0.06140350877192982   |  0.22101373435409918
1881      |  1419                   |  6076   |  6076           |  DAY_16     |  8                        |  0                     |  1                      |  1                  |  22                   |  23            |  957          |  450.43782654127483  |  471.07912844036696         |  76.15909786907412              |  98        |  0.21756609730692766       |  1                     |  1                      |  0     |  18    |  0     |  0     |  80    |  G          |  T          |  CT             |  0.1836734693877551        |  0.1836734693877551        |  0.225                 |  0.4769182703436179
1882      |  1419                   |  6076   |  6076           |  DAY_17A    |  8                        |  0                     |  1                      |  1                  |  22                   |  23            |  957          |  110.49738766980147  |  114.24320827943079         |  12.051328175500755             |  66        |  0.5972991886217918        |  1                     |  0                      |  0     |  4     |  0     |  0     |  62    |  G          |  T          |  CT             |  0.06060606060606061       |  0.06060606060606061       |  0.06451612903225806   |  0.22863187358286127
1883      |  1419                   |  6076   |  6076           |  DAY_18     |  8                        |  0                     |  1                      |  1                  |  22                   |  23            |  957          |  514.3531870428423   |  539.0582857142857          |  86.05041729742642              |  124       |  0.2410794821995953        |  1                     |  1                      |  0     |  15    |  0     |  0     |  109   |  G          |  T          |  CT             |  0.12096774193548387       |  0.12096774193548387       |  0.13761467889908258   |  0.36884872544769964
(...)     |  (...)                  |  (...)  |  (...)          |  (...)      |  (...)                    |  (...)                 |  (...)                  |  (...)              |  (...)                |  (...)         |  (...)        |  (...)               |  (...)                      |  (...)                          |  (...)     |  (...)                     |  (...)                 |  (...)                  |  (...) |  (...) |  (...) |  (...) |  (...) |  (...)      |  (...)      |  (...)          |  (...)                     |  (...)                     |  (...)                 |  (...)

1. **entry_id** refers to the unique id for the line in the output file. It is the only column that contains a unique id for each line.

2. **unique_pos_identifier** refers to the unique identifier for a given nucleotide position. Since each sample in the profile database can report variability for every nucleotide position, a **unique_pos_identifier** can appear in the file as many times as the number of the samples in the analysis. This column can be used to pull frequencies of nucleotides for a given nucleotide position from all samples.

3. **pos** refers to the nucleotide position in the split.

4. **pos_in_contig** refers to the nucleotide position in the contig (why is this called *pos_in_contig*, and the one before is not called *pos_in_split*? Well, we have been wondering about that for a long time, too circa 2015).

5. **contig_name** refers to the contig name as it appears in the contigs database. Stored if `--include-contig-names` flag is provided.

6. **split_name** refers to the split name.  Stored if `--include-split-names` flag is provided.

7. **sample_id** corresponds to the sample name a given particular line is reported from. This column allows the linking of SNVs and the sample(s) they were identified from.

8. **corresponding_gene_call** refers to a unique gene caller id (`-1`, if the position falls out of a gene call).

9. **in_partial_gene_call** indicates whether the gene call is incomplete (i.e., starts with a start codon, stops with a stop codon, etc). `1` if incomplete, `0` if both start and stop positions are detected, or if the position is not in a gene.

10. **in_complete_gene_call** indicates the gene completion status. `1` for complete, `0` if incomplete, or if the position is not in a gene.

11. **base_pos_in_codon** refers to the position of the nucleotide in a codon. `1`, `2` or `3` for codon positions, `-1` if the position is not in a detected gene.

12. **codon_order_in_gene / codon_index** refers to the order of the codon in the gene call, counting from 0 as the start position. For example, the starting methionine in a translated gene has a **codon_order_in_gene** value of 0. In anvi'o version 5 and greater, this has been renamed **codon_index**. `-1` if the position is not in a called gene. 

13. **codon_number** refers to the order of the codon in the gene call, counting from 1 as the start position. For example, the starting methionine in a translated gene has a **codon_order_in_gene** value of 1. `-1` if the position is not in a called gene. **codon_number** = **codon_index** + 1.

14. **gene_length** is the length of the gene the nucleotide position falls in, measured in base pairs. -1 if position is not in a called gene.

15. **gene_coverage** is the average coverage (number of recruited reads per nucleotide position) of the gene for which the nucleotide position lies in. -1 if position is not in a called gene.

16. **non_outlier_gene_coverage** is a measure of gene coverage that excludes nucleotide positions outlier coverage values. The criterion is that the average coverage is computed over nucleotides that are within plus or minus 1.5 times the MAD ([median absolute deviation](https://en.wikipedia.org/wiki/Median_absolute_deviation)) of the median coverage value. We took this idea from [this paper](https://www.sciencedirect.com/science/article/pii/S0022103113000668) and we like it very much because it does a good job at excluding sequence motifs that are conserved between species/populations and therefore recruit unrealistically high coverage values. -1 if position is not in a called gene.

15. **non_outlier_gene_coverage_std** is the standard deviation of coverage values for all nucleotide positions in the gene that passed the MAD filter described in the definition of **non_outlier_gene_coverage**. -1 if position is not in a called gene.

17. **coverage** refers to the coverage, the number of recruited reads mapping to this position.

18. **mean_normalized_coverage** refers to the **coverage** divided by **gene_coverage**.

19. **cov_outlier_in_split** indicates whether the coverage of this position is marked as an outlier compared to all other positions in the split. `1` if outlier, and `0` if not.

20. **cov_outlier_in_contig** has the same purpose with the one above, except it is at the contig-level.

21. **A** refers to the number of mapped reads covering this position with a `A`.

22. **C** refers to the number of mapped reads covering this position with a `C`.

23. **G** refers to the number of mapped reads covering this position with a `G`.

24. **N** refers to the number of mapped reads covering this position with an ambiguous base.

25. **T** refers to the number of mapped reads covering this position with a `T`.

26. **reference** refers to the reference nucleotide in the mapping context.

27. **consensus** refers to the most frequent nucleotide mapping to this position. This column is used to define the **departure_from_consensus** ratio.

28. **competing_nts** refers to the two most represented nucleotides. Note that if competing nucleotides are stored as `CT`, it does not necessarily mean that `C` occurs more than `T`, since they are ordered alphabetically. For positions that do not have any variation, competing nucleotides will contain two identical nucleotides. I hope here you are asking yourself here "*why would a position without any variation would appear in this table?*". The answer is `--quince-mode`. See below.

29. **departure_from_reference** refers to the ratio of nucleotides in a given position that diverge from the reference nucleotide. If a position with a coverage of 100X has the nucleotide `A` in the reference, the departure from reference would be `0.2` if the frequency of mapping nucleotides are as follows: `A: 80X`, `T: 0X`, `C: 12X`, and `G: 8X`. Note that the departure from reference can dramatically change across samples, revealing subtle differences at the single nucleotide level.

30. **departure_from_consensus** refers to the ratio of nucleotides in a given position that diverge from the most frequent nucleotide. If the nucleotide frequencies for a given position are, `A: 10X`, `T: 20X`, `C: 30X`, and `G: 70X`, then the departure from reference would be `0.46` (from `(10X + 20X + 30X) / (10X + 20X + 30X + 70X)`). Compared to the **departure from reference**, this column relies less on the reference sequence and thus might better reflect the variation in the environment regardless of the context you are mapping against.

31. **n2n1ratio** refers to the ratio of the second most frequent nucleotide to the consensus nucleotide. If the frequency of nucleotides mapped to this position are `A: 10X`, `T: 20X`, `C: 30X`, and `G: 70X`, then the n2n1ratio would be `0.42`, ... and you would be as close as two orders of magnitude to the answer to the ultimate question of life, the universe, and everything (which, in an ideal world, should include an answer to your research question, too, you lucky you (but then is there really an answer to your research question? BAM! (If you survived this, you can survive the rest of this tutorial. Please carry on))).

32. **entropy** is the site entropy of a nucleotide position. Read about that [here](https://en.wikipedia.org/wiki/Entropy_(information_theory)). People always say, "entropy is a measure of information" as if you're then supposed to go, "OH! Now it makes sense". To me, entropy quantifies the degree of disagreement ("OH! Now it makes sense"), much like departure from consensus. If the entropy is 0, it means every read at that position was identical. For example, if you read my spiel above about frequency vectors, the frequency vector of `<A=1, C=0, G=0, T=0>` has an entropy of 0 because there is no variation. On the flipside, the entropy is maximized when all frequencies are equal, i.e. when `<A=0.25, C=0.25, G=0.25, T=0.25>`.

33. **kullback_leibler_divergence_raw** is a weighted version of the site entropy that measures how different a metagenome is compared to the average across all metagenomes. Consider for example, a SNV with a frequncy vector `<A=0.25, C=0.25, G=0.25, T=0.25>`. This SNV has maximal **entropy**, however if all other metagenomes at this site also have the same frequency vector, the **kullback_leibler_diverence_raw** will be 0. On the otherhand, if the frequency vectors of all other metagenomes are `<A=1, C=0, G=0, T=0>`, this SNV stands out like a sore thumb and will therefore have a maximal **kullback_leibler_diverence_raw**. For more info on Kullback-Leibler Divergence, check out the [wiki](https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence). This output is only present when `--quince-mode` flag is provided.

34. **kullback_leibler_divergence_normalized** is just like **kullback_leibler_divergence_raw**, however with a small and important difference. **kullback_leibler_divergence_raw** compares the metagenome of interest to the "average" metagenome in your dataset by calculating a frequency vector for the "average" metagenome. This average is calculated by summing up the counts of `A`, `C`, `G`, and `T` across all metagenomes and then dividing each of the 4 counts by the sum of coverage values across all metagenomes. The consequence of this is that metagenomes with very large coverage values bias the average frequency vector towards their own frequency vector. To equally weight each metagenome, **kullback_leibler_divergence_normalized** calculates the frequency vector for each metagenome separately, and then averages those frequency vectors so that each metagenome gives an equal weight to the average frequency vector. This output is only present when `--quince-mode` flag is provided.

## Matrix output (those unique to SCVs)

As mentioned above, a prerequisite for reporting SCVs is having the flag `--profile-SCVs` when your single or merged anvi'o profile databases were created. After that, it is as simple as running `anvi-gen-variability-profile` with the parameter `--engine CDN`. Besides obvious changes such as the columns `A`, `C`, `G`, `T` being replaced with `AAA`, `AAC`, etc., the matrix contains some extra information. Here is a sample output:

entry_id  |  unique_pos_identifier  |  sample_id   |  corresponding_gene_call  |  codon_order_in_gene  |  codon_number  |  gene_length  |  gene_coverage       |  non_outlier_gene_coverage  |  non_outlier_gene_coverage_std  |  coverage  |  mean_normalized_coverage  |  AAA   |  AAC  | (...) |  TTC  |  TTG   |  TTT  |  reference  |  consensus  |  competing_codons  |  departure_from_reference  |  departure_from_consensus  |  n2n1ratio             |  entropy              |  synonymity
----------|-------------------------|--------------|---------------------------|-----------------------|----------------|---------------|----------------------|-----------------------------|---------------------------------|------------|----------------------------|--------|-------|-------|-------|--------|-------|-------------|-------------|--------------------|----------------------------|----------------------------|------------------------|-----------------------|--------------------
30        |  25                     |  SRR5405970  |  6                        |  21                   |  22            |  966          |  25.25879917184265   |  24.774764150943398         |  3.428167856163931              |  30        |  1.187704918032787         |  0     |  0    | (...) |  0    |  0     |  0    |  CAG        |  CAC        |  CACCAC            |  1.0                       |  0.0                       |  0.0                   |  0.0                  |  1.0
31        |  26                     |  SRR5405970  |  6                        |  44                   |  45            |  966          |  25.25879917184265   |  24.774764150943398         |  3.428167856163931              |  26        |  1.029344262295082         |  0     |  0    | (...) |  0    |  0     |  0    |  CGC        |  CGT        |  CGCCGT            |  0.8076923076923077        |  0.19230769230769232       |  0.23809523809523808   |  0.4895518936229288   |  1.0
43        |  38                     |  SRR5405970  |  6                        |  283                  |  284           |  966          |  25.25879917184265   |  24.774764150943398         |  3.428167856163931              |  17        |  0.6730327868852459        |  0     |  0    | (...) |  0    |  0     |  0    |  ATT        |  ATC        |  ATCATT            |  0.8823529411764706        |  0.11764705882352941       |  0.13333333333333333   |  0.362210557135449    |  1.0
44        |  39                     |  SRR5405970  |  11                       |  6                    |  7             |  1236         |  23.100323624595468  |  22.929310344827588         |  4.19710526766056               |  32        |  1.3852619781451387        |  0     |  0    | (...) |  0    |  0     |  0    |  ACC        |  ACT        |  ACCACT            |  0.9375                    |  0.0625                    |  0.06666666666666667   |  0.2337916587064593   |  1.0
45        |  40                     |  SRR5405970  |  11                       |  29                   |  30            |  1236         |  23.100323624595468  |  22.929310344827588         |  4.19710526766056               |  25        |  1.0822359204258896        |  25    |  0    | (...) |  0    |  0     |  0    |  GAA        |  AAA        |  AAAAAA            |  1.0                       |  0.0                       |  0.0                   |  0.0                  |  1.0
56        |  51                     |  SRR5405970  |  11                       |  321                  |  322           |  1236         |  23.100323624595468  |  22.929310344827588         |  4.19710526766056               |  21        |  0.9090781731577473        |  0     |  0    | (...) |  0    |  0     |  0    |  GAT        |  AAT        |  AATGAT            |  0.8095238095238095        |  0.19047619047619047       |  0.23529411764705882   |  0.48691270946460224  |  0.8176943699731903
57        |  52                     |  SRR5405970  |  11                       |  331                  |  332           |  1236         |  23.100323624595468  |  22.929310344827588         |  4.19710526766056               |  21        |  0.9090781731577473        |  0     |  0    | (...) |  0    |  0     |  0    |  ACA        |  ACG        |  ACGACG            |  1.0                       |  0.0                       |  0.0                   |  0.0                  |  1.0
0         |  0                      |  SRR5405970  |  36                       |  38                   |  39            |  2136         |  23.000936329588015  |  22.270707070707072         |  5.975408401160299              |  33        |  1.4347242010991248        |  0     |  0    | (...) |  4    |  0     |  29   |  TTC        |  TTT        |  TTCTTT            |  0.8787878787878788        |  0.12121212121212122       |  0.13793103448275862   |  0.3693331216456527   |  1.0
58        |  53                     |  SRR5405970  |  71                       |  12                   |  13            |  639          |  21.77151799687011   |  21.132246376811594         |  6.405860373728439              |  28        |  1.2860839562967223        |  0     |  0    | (...) |  0    |  0     |  0    |  GGC        |  GGC        |  GGCGGT            |  0.07142857142857142       |  0.07142857142857142       |  0.07692307692307693   |  0.25731864054383163  |  1.0
(...)     |  (...)                  |  (...)       |  (...)                    |  (...)                |  (...)         |  (...)        |  (...)               |  (...)                      |  (...)                          |  (...)     |  (...)                     |  (...) | (...) | (...) | (...) |  (...) | (...) |  (...)      |  (...)      |  (...)             |  (...)                     |  (...)                     |  (...)                 |  (...)                |  (...)

The additional columns are described here:

1. **competing_codons** is just like **competing_nts**, except it refers to the two most represented codons.

2. **synonymity** describes the degree of synonymous vs non-synonymous changes in a SCV. It is calculated by taking each read that mapped to a SCV and comparing it to every other read, marking them as either synonymous or non-synonymous. Synonymity is simply the fraction of synonymous pairs. For example, if all reads are `AAA`, then the synonymity is 1 because `AAA` to `AAA` is a synonymous mutation. If half the reads are `AAA` and half are `ACA`, the synonymity is 0.5. Wait, what? Not 0? That's because every `AAA` will pair with `AAA` half of the time, and `ACA` half the time. Likewise, every `ACA` will pair with `ACA` half of the time, and `AAA` half the time. We are currently thinking of better metrics to quantify synonymity, so please share your idea if you have a metric that you think could work.


## Matrix output (those unique to SAAVs)

Just like for SCVs, a prerequisite for reporting SCVs is having the flag `--profile-SCVs` (that's not a typo, SAAVs are computed from SCVs) when your single or merged anvi'o profile databases were created. After that, it is as simple as running `anvi-gen-variability-profile` with the parameter `--engine AA`. Besides obvious changes such as the columns `A`, `C`, `G`, `T` being replaced with `Ala`, `Arg`, etc., the matrix contains some extra information. Here is a sample output:

entry_id  |  unique_pos_identifier  |  sample_id   |  corresponding_gene_call  |  codon_order_in_gene  |  codon_number  |  gene_length  |  gene_coverage       |  non_outlier_gene_coverage  |  non_outlier_gene_coverage_std  |  coverage  |  mean_normalized_coverage  |  Ala   |  Arg   |  (...) |  Tyr   |  Val   |  reference  |  consensus  |  competing_aas  |  departure_from_reference  |  departure_from_consensus  |  n2n1ratio             |  entropy              |  BLOSUM62  |  BLOSUM62_weighted    |  BLOSUM90  |  BLOSUM90_weighted
----------|-------------------------|--------------|---------------------------|-----------------------|----------------|---------------|----------------------|-----------------------------|---------------------------------|------------|----------------------------|--------|--------|--------|--------|--------|-------------|-------------|-----------------|----------------------------|----------------------------|------------------------|-----------------------|------------|-----------------------|------------|---------------------
30        |  25                     |  SRR5405970  |  6                        |  21                   |  22            |  966          |  25.25879917184265   |  24.774764150943398         |  3.428167856163931              |  30        |  1.187704918032787         |  0     |  0     |  (...) |  0     |  0     |  Gln        |  His        |  HisHis         |  1.0                       |  0.0                       |  0.0                   |  0.0                  |  8.0       |                       |  8.0       |
31        |  26                     |  SRR5405970  |  6                        |  44                   |  45            |  966          |  25.25879917184265   |  24.774764150943398         |  3.428167856163931              |  26        |  1.029344262295082         |  0     |  26    |  (...) |  0     |  0     |  Arg        |  Arg        |  ArgArg         |  0.0                       |  0.0                       |  0.0                   |  0.0                  |  5.0       |                       |  6.0       |
43        |  38                     |  SRR5405970  |  6                        |  283                  |  284           |  966          |  25.25879917184265   |  24.774764150943398         |  3.428167856163931              |  17        |  0.6730327868852459        |  0     |  0     |  (...) |  0     |  0     |  Ile        |  Ile        |  IleIle         |  0.0                       |  0.0                       |  0.0                   |  0.0                  |  4.0       |                       |  5.0       |
44        |  39                     |  SRR5405970  |  11                       |  6                    |  7             |  1236         |  23.100323624595468  |  22.929310344827588         |  4.19710526766056               |  32        |  1.3852619781451387        |  0     |  0     |  (...) |  0     |  0     |  Thr        |  Thr        |  ThrThr         |  0.0                       |  0.0                       |  0.0                   |  0.0                  |  5.0       |                       |  6.0       |
45        |  40                     |  SRR5405970  |  11                       |  29                   |  30            |  1236         |  23.100323624595468  |  22.929310344827588         |  4.19710526766056               |  25        |  1.0822359204258896        |  0     |  0     |  (...) |  0     |  0     |  Glu        |  Lys        |  LysLys         |  1.0                       |  0.0                       |  0.0                   |  0.0                  |  5.0       |                       |  6.0       |
56        |  51                     |  SRR5405970  |  11                       |  321                  |  322           |  1236         |  23.100323624595468  |  22.929310344827588         |  4.19710526766056               |  21        |  0.9090781731577473        |  0     |  0     |  (...) |  0     |  0     |  Asp        |  Asn        |  AsnAsp         |  0.8095238095238095        |  0.19047619047619047       |  0.23529411764705882   |  0.48691270946460224  |  1.0       |  0.9999999999999999   |  1.0       |  0.9999999999999999
57        |  52                     |  SRR5405970  |  11                       |  331                  |  332           |  1236         |  23.100323624595468  |  22.929310344827588         |  4.19710526766056               |  21        |  0.9090781731577473        |  0     |  0     |  (...) |  0     |  0     |  Thr        |  Thr        |  ThrThr         |  0.0                       |  0.0                       |  0.0                   |  0.0                  |  5.0       |                       |  6.0       |
0         |  0                      |  SRR5405970  |  36                       |  38                   |  39            |  2136         |  23.000936329588015  |  22.270707070707072         |  5.975408401160299              |  33        |  1.4347242010991248        |  0     |  0     |  (...) |  0     |  0     |  Phe        |  Phe        |  PhePhe         |  0.0                       |  0.0                       |  0.0                   |  0.0                  |  6.0       |                       |  7.0       |
58        |  53                     |  SRR5405970  |  71                       |  12                   |  13            |  639          |  21.77151799687011   |  21.132246376811594         |  6.405860373728439              |  28        |  1.2860839562967223        |  0     |  0     |  (...) |  0     |  0     |  Gly        |  Gly        |  GlyGly         |  0.0                       |  0.0                       |  0.0                   |  0.0                  |  6.0       |                       |  6.0       |
(...)     |  (...)                  |  (...)       |  (...)                    |  (...)                |  (...)         |  (...)        |  (...)               |  (...)                      |  (...)                          |  (...)     |  (...)                     |  (...) |  (...) |  (...) |  (...) |  (...) |  (...)      |  (...)      |  (...)          |  (...)                     |  (...)                     |  (...)                 |  (...)                |  (...)     |  (...)                |  (...)     |

The additional columns are described here:

1. **competing_aas** is just like **competing_nts**, except it refers to the two most represented amino acids.

2. **BLOSUM62 / BLOSUM90** are both specific types of BLOSUM (BLOcks SUbstitution Matrix) matrices. BLOSUM matrices describe the interchangeability of two amino acids based on massive datasets of multiple sequence alignments. For example, valine-isoleucine is found in these multiple sequence alignments to be extremely common, and so the valine-isoleucine entry has a score of 3 in the BLOSUM62 matrix and also in the BLOSUM90. In contrast, tryptophan-glycine has an entry of -2 and -4, respectively. Thus, the more positive, the more common, and the more negative, the less common. The value of each entry in these columns corresponds to the BLOSUM entry between the **competing_aas**. This isn't the place to get into the details of BLOSUM matrices, but the [original paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC50453/) does an excellent job with an easy to follow example if you want details. What is important to distinguish is the difference between BLOSUM62 and BLOSUM90. 62 and 90 refer to the percent identity cutoffs for each multiple sequence alignment. What this means is that BLOSUM62 generates its substitution scores from a more divergent collection of species, whereas BLOSUM90 has stricter cutoffs. This means that relative to BLOSUM90, BLOSUM62 will overestimate the likelihood of substitution between chemically dissimilar amino acids, which is what we observe in the case of tryptophan and glycine: BLOSUM62 has a score of -2 and BLOSUM90 has a score of -4. Since SAAVs belong to highly genetically similar populations, BLOSUM90 is probably a more accurate representation of substitution rates.

3. **BLOSUM62_weighted / BLOSUM90_weighted** are weighted versions of **BLOSUM62** and **BLOSUM90**. Rather than taking the 2 most common amino acids (i.e. the **competing_aas**), the weighted BLOSUM scores determine the mean BLOSUM scores when comparing all reads against all reads. By design, we chose to exclude BLOSUM scores of a read with itself, so that only different amino acids were compared against one another. If there is no variability for an entry, the weighted BLOSUM scores are therefore undefined.


## Matrix output (those unique to structure database integration)

So in addition to the above parameters, your variability table entries can be further extended with protein structure annotation if you have a structure database and pass it to `anvi-gen-variability-profile` with the parameter `-s <your_structure_database.db>`. You can read about creating your own structure database here: [structural biology with anvi'o](http://merenlab.org/2018/09/04/structural-biology-with-anvio/). These are the columns that can extend your variability output with:


1. **sec_struct** is secondary structure of the corresponding protein structure in your structure database. The program responsible for these annotations is [DSSP](http://merenlab.org/2016/06/18/installing-third-party-software/#dssp)). The 8 classes of secondary structure are 3-10 helix (G), alpha helix (H), pi-helix (I), beta bridge (B), beta bulge (E), turns (T), high curvature (S), and generic loop (C). These 8 classes fall into 3 broader categories you may be more familiar with: helix (G, H, and I), strand (E, and B), and loop (S, T, and C) (taken from [wikipedia](https://en.wikipedia.org/wiki/DSSP_(hydrogen_bond_estimation_algorithm))).

2. **rel_solvent_acc** is a measure of how accessible the reference residue is to water. To normalize for the fact that some residues have more surface area then others, **rel_solvent_acc** is the fraction of the residues total surface area that is accessible to water. Solvent accessibility is an extremely important determinant in the permissibility of amino acid-impacting mutations.

3. **phi / psi** are two bond angles that describe the orientation of a residue. They are the x- and y-coordinates of a [Ramachandran plot](https://en.wikipedia.org/wiki/Ramachandran_plot) and are useful to determine the overall structure of protein. For example, a helix-rich protein will be enriched in phi and psi angles of around -60 and -40 degrees, respectively.

4. **contact_numbers** is a list of neighboring residues determined to be in physical contact with the residue's side chain. For example, if **contact_numbers** is equal to "11,12,56,59", then residues with **codon_numbers** 11, 12, 56, and 59 are all in contact with the residue. The criterion for being "in contact" is an adopted dogma that the alpha carbons from each of the residues are within 6 angstroms of one another.


## Parameters to refine the output

This program allows the researcher to specify _filters_ and refine reported variable nucleotide positions with a rich array of parameters:

- With `anvi-gen-variability-profile` you can request all variable positions identified in a genome bin stored in a collection using the `--collection-name` and `--bin-id` parameters. Or alternatively you can focus on an arbitrary list of splits using `--splits-of-interest` parameter. You can also focus on a specific gene calls using the `--genes-of-interest` parameter which is a file of gene caller ids. If you just care about a couple of genes, you can provide them with the `--gene-caller-ids` parameter (e.g. `--gene-caller-ids 23,235,867`).

- As already mentioned, you can profile SNVs, SCVs, or SAAVs by using `--engine NT`, `--engine AA`, or `--engine CDN`, respectively.

- By default, the program will report variants from every sample found in the profile database. You can use the `--samples-of-interest` parameter to define which samples should be included in the analysis.

- The output file will not have columns to display in which split or contig a given SNV was identified. But you can change that behavior using the `--include-contig-names` and/or `--include-split-names` flags.

- You can use the `--num-positions-from-each-split` to define a maximum number of SNVs to be reported from each split. If a given split has more SNVs than this number, the program will randomly select a matching number of them.

- You can set a minimum coverage value for each SNV using the parameter `--min-coverage-in-each-sample`. Any SNV that is covered less than this value would be omitted from the output.

- You can filter nucleotide positions based upon frequencies using parameters `--min-departure-from-reference`, `--max-departure-from-reference`, `--min-departure-from-consensus`, or `--max-departure-from-consensus`.

- You can remove nucleotide positions from the final report that show variation in less than any number of samples using the `--min-occurrence` parameter.

- You can include structural information in your output by providing a structure database with `-s <your_structure.db>`. Entries without corresponding structures will not be annotated, and can be removed from the output with `--only-if-structure`.

- One of the very important flags is the flag `--quince-mode`. If you read the previous section you know that the base frequencies are reported for a given nucleotide position in a sample *only* if that position has variation. But it doesn't work well when you want to be able to say "give me back SNVs only if the nucleotide position they occur into is covered more than X in every sample". But not every sample will be in the variability profile, hence, the program will not have access to the coverage scores for nucleotide positions for samples with 0 variation at those positions. When you use the `--quince-mode` flag, the program goes back to the auxiliary files generated during profiling, and recovers coverage values for nucleotide positions in sampels for each **unique_pos_identifier**. In other words, if there are 10 samples in a dataset, and a given position has been reported as a variable site during the profiling in only one of those samples, there will be no information stored in the database for the remaining 9. When this flag is used, base frequencies will be recovered. It will take considerably longer to report when this flag is on, and the use of it will increase the file size dramatically, however, it may be necessary to use for some statistical approaches (as well as for some beautiful visualizations).

{:.notice}
**Fact**: The answer is "*yes*", if you are asking yourself "*is this 'quince' in `--quince-mode` the 'Quince' in '[Christopher Quince](http://www2.warwick.ac.uk/fac/med/staff/cquince/)'?*". We added this flag because he complained about the fact that this functionality was not available in anvi'o at the time.

{:.notice}
**Fact**: The answer is "*very very likely*", if you are now asking yourself "*would I get my own flag if I complain about a lacking functionality in anvi'o?*". We, the anvi'o developers, are often a civil bunch of people.

---

We certainly hope the variant characterization framework in anvi'o can be useful for you. Please don't hesitate get in touch with us with questions, or suggestions.

If you are interested, you can find an easy-to-follow tutorial here, demonstrating a variant analysis on the Infant Gut Dataset: [http://merenlab.org/tutorials/infant-gut/#profiling-snvs-in-a-bin]({{ site.url }}/tutorials/infant-gut/#profiling-snvs-in-a-bin). 

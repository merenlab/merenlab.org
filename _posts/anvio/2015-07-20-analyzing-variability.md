---
layout: post
authors: [meren,tom]
title: "Analyzing single nucleotide variations (SNVs) with anvi'o"
excerpt: "Exploring micro-diversity patterns using for deeper insights into ecology"
modified: 2016-07-21
tags: []
categories: [anvio]
comments: true
---

{:.notice}
This is more of a theoretical tutorial. For a more practical one on the same topic, please visit [http://merenlab.org/tutorials/infant-gut/#profiling-snvs-in-a-bin]({{ site.url }}/tutorials/infant-gut/#profiling-snvs-in-a-bin)

{% include _toc.html %}

We use the term "**microbial clouds**" to describe an assemblage of co-existing microbial genomes in an environment that are similar enough to map to the context of the same reference genome (this concept is widely defined as "microbial populations", although the definition of the term [population](https://en.wikipedia.org/wiki/Population) makes it somewhat irrelevant to microbiology).

Cells within a microbial cloud (i.e., all cells that would have classified as the same "species" or "strain" for whatever these terms mean to you) will share the vast majority of their genomes in the sequence space. Hence, a consensus contig obtained through the assembly of metagenomic reads from this environment, or a contig that originates from an isolate cultured from this environment, will recruit short reads through mapping from many other member cells in a given cloud.

Metagenomic short reads from a cloud mapping to the same context can have nucleotide positions that systematically differ from the reference context depending on, 

* The heterogeneity of the microbial cloud,
* the fraction of the cloud that can be targeted by the reference genome through mapping, 
* and/or the stringency of mapping.

Since one of the mechanisms for the diversification of genomes operates at the single-nucleotide level, ability to characterize the heterogeneity of a microbial cloud nucleotide by nucleotide can be very critical for a more complete understanding of the environmental forces that affect communities, and adaptive strategies microbes rely on to survive in environments they reside. Anvi'o can *de novo* characterize single nucleotide variations (SNVs) in a given microbial cloud using mapping results, and allow researchers to delineate subtle ecological niche and ecotypes, and quantify the level of heterogeneity in a microbial cloud.

This article will describe some key aspects of the anvi'o workflow for the recovery, profiling, and characterization of SNVs for high-resolution genomics using easy-to-use anvi'o programs.

# Defining "*single nucleotide variation*"

Depending on the complexity of a microbial cloud (i.e., the extent of monoclonality in it, or the lack thereof), metagenomic short reads that map to a reference context (i.e., a contig from a MAG or a cultivar genome) can generate one or more mismatches (unless the mapping software does not require a 100% sequence identity for alignment).

Here we define 'variation' as the extent of disagreement between the aligned nucleotides that map to the same context. While the source of a mismatch can be artificial (i.e., due to random sequencing or PCR errors, or non-specific mapping of local alignments), it may also represent ecologically informative variation. The context is a single nucleotide position from the reference sequence. This disagreement is used to differentiate SNVs (generally in a minority of positions) from stable nucleotide positions upon 'base frequencies'. Note that SNVs are generally dominated by two nucleotides (e.g., 60% of `A` and 40% of `T`) but can also represent three or four nucleotides in relatively high proportions.

Although the importance of it is probably almost common sense to anyone who understands the steps we follow to recover our contigs for a given environment, and although many amazing studies exist including the ones from Jill Banfield's lab that made use of SNV patterns in metagenomic settings such as [this one](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0060177), working with subtle nucleotide variation in cultivar genomes or metagenomic assembly outputs hasn't been so practical to access. There are **many** challenges to address if one wants to focus on subtle things in complex environmental sequencing datasets, and in order to set the stage with reasonable expectations we would like to clarify that anvi'o will not exactly offer you a magical wand. However, anvi'o _will_ empower you to work with this aspect of your datasets by making the recovery and query of this information more or less straightforward. We hope that at the end of this post it will be clear to you how this can be done.

---

Just so you have a mental picture of how anvi'o visualizes SNVs, here is a polished screenshot from the interface:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2015-07-20-analyzing-variability/example.png"><img src="{{ site.url }}/images/anvio/2015-07-20-analyzing-variability/example.png" style="border: none; width: 100%;" /></a>
</div>

What you see is the coverage and base frequencies in reported SNVs in a single contig across two samples. Where these samples are coming from are not really relevant to the rest of this post. But usually, just like it is the case in this example, SNVs and their base frequencies do not seem to be random, and they can be extremely reproducible if the population structures are really similar in multiple samples (i.e., in biological replicates).

# The anvi'o way

Anvi'o gives access to information about single nucleotide variants in two steps.

The first step is to identify SNVs and reporting them _for each sample separately_. This step takes place during the profiling of a given BAM file with the program `anvi-profile`. 

The second step is to interpret the ecological significance of sample-specific SNVs _across samples_ using the program `anvi-gen-variability-profile`. The first step is agnostic to the experimental design, and/or the SNVs in other samples: it simply does its best to find and report variation. The second step is where the user's input about the experimental design, and their questions come into the play.

The following two chapters in this article will detail these two steps, and the third chapter will explain the analysis of variability using the Infant Gut Dataset that is first appeared in [this study](http://genome.cshlp.org/content/23/1/111.long).

Everything you might need to replicate this analysis and to re-generate the third figure in the anvi'o paper is [here]({{ site.url }}/data/). We hope that we will manage to give enough details that you will be able to modify the recipe and start exploring patterns in your own datasets if you wish.

Please do not hesitate to ask questions, make suggestions, and/or start discussions regarding this topic here in the comments section down below, or at [the issues section of the anvi'o repository on Github](https://github.com/meren/anvio/issues), or by directly contacting [us]({{ site.url }}/people/).

# *De novo* characterization and reporting of SNVs

During the profiling step that should be done for each sample separately, anvi'o checks the composition of nucleotides that map to each reference position, and keeps track of variation in the form of 'base frequencies'.

Note that the reporting of these base frequencies is a little tricky: If you report base frequencies for each position when there is any variation at all, then the number of reported positions would be enormous. For instance, if you have a nucleotide position with 500X coverage, you can safely assume that there will always be _some_ nucleotides that do not match to the consensus nucleotide for that position --whether it is due to biology, or random sequencing errors, or other types of relevant or irrelevant sources of variation. For the sake of better management of available resources, and to be very quick, anvi'o (with its default settings) ___does not___ report variation in every single nucleotide position during profiling.

Instead, it relies on the following conservative heuristic to identify SNVs and report the variation only at these nucleotide positions:
 
<div class='centerimg'>
<a href='https://www.desmos.com/calculator/qwocua4zi5'><img src='{{ site.url }}/images/anvio/2015-07-20-analyzing-variability/function.png' style='border: none; width: 200px;' /></a>
</div>

where, `x` represents the coverage, and `b`, `m`, and `c` represent the model parameters equal to 3, 1.45, and 0.05, respectively. Assuming `n1` and `n2` represent the frequency of the most frequent and the second most frequent bases in a given nucleotide position (see the table), base frequencies are reported only if `n2/n1 > y` criterion is satisfied for a given coverage of `x`. It is this simple (and ugly, in a sense). But briefly, this approach sets a dynamic baseline for minimum variation required for reporting _as a function of coverage depth_. According to this heuristic, `y` would be 0.29 for 20X coverage, 0.13 for 50X coverage, 0.08 for 100X coverage, and ~0.05 for very large values of coverage as `y` approaches to `c`. The goal here is to lessen the impact of sequencing and mapping errors in reported frequencies, and it does it's boring job.

This computation- and storage-efficient strategy reports a relatively short list of sample-specific SNVs that one can really trust. However, the user has always the option to instruct the profiler to store _all observed frequencies_ by declaring this intention with `--report-variability-full` flag for more statistically appropriate downstream analyses. We are of course talking about the type of analyses Christopher Quince-like people would like to do.


# Generating a SNV profile

Once SNVs are stored in a single or merged anvi'o profile database, it is time to decide how to process that information, and export a scrutinized table to make sense of them. To interpret the ecological significance of sample-specific variable positions across samples, anvi'o installs a helper program called `anvi-gen-variability-profile`.

In this section will first describe the output file structure for SNV profiles, and then describe the parameters of `anvi-gen-variability-profile`.

## The output matrix

The output generated by `anvi-gen-variability-profile` contains 24 to 26 columns per nucleotide position (2 columns are optional). Each line in the output file represents one nucleotide position. Here is the first 10 lines of an example SNV profile for the *E. facealis* bin in the Infant Gut Dataset so you have an idea of the file format, followed by the description of each column:

|entry_id|unique_pos_identifier|sample_id|pos|pos_in_contig|corresponding_gene_call|in_partial_gene_call|in_complete_gene_call|base_pos_in_codon|codon_order_in_gene|coverage|cov_outlier_in_split|cov_outlier_in_contig|departure_from_reference|competing_nts|reference|A|T|C|G|N|consensus|departure_from_consensus|n2n1ratio|
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
|365|0|DAY_17A|1600|1600|-1|0|0|0|-1|42|0|0|0.428571428571|CT|T|0|24|14|4|0|T|0.428571428571|0.583333333333|
|366|1|DAY_17A|1602|1602|-1|0|0|0|-1|38|0|0|0.210526315789|CT|T|0|30|8|0|0|T|0.210526315789|0.266666666667|
|367|2|DAY_17A|1603|1603|-1|0|0|0|-1|36|0|0|0.277777777778|AG|G|10|0|0|26|0|G|0.277777777778|0.384615384615|
|368|3|DAY_17A|1606|1606|-1|0|0|0|-1|28|0|0|0.142857142857|AT|A|24|4|0|0|0|A|0.142857142857|0.166666666667|
|369|4|DAY_17A|1596|1596|-1|0|0|0|-1|52|0|0|0.153846153846|AT|A|44|4|0|4|0|A|0.153846153846|0.0909090909091|
|370|5|DAY_17A|624|624|11188|0|1|2|88|324|0|0|0.104938271605|AG|A|290|0|2|32|0|A|0.104938271605|0.110344827586|
|371|6|DAY_17A|627|627|11188|0|1|2|89|326|0|0|0.539877300613|AG|G|176|0|0|150|0|A|0.460122699387|0.852272727273|
|372|7|DAY_17A|821|821|11188|0|1|1|154|104|0|0|0.0961538461538|CT|C|0|10|94|0|0|C|0.0961538461538|0.106382978723|
|373|8|DAY_17A|470|470|11188|0|1|1|37|184|0|0|0.184782608696|AG|G|34|0|0|150|0|G|0.184782608696|0.226666666667|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|

1. **entry_id** refers to the unique id for the line in the output file. It is the only column that contains a unique id for each line.

2. **unique_pos_identifier** refers to the unique identifier for a given nucleotide position. Since each sample in the profile database can report variability for every nucleotide position, a **unique_pos_identifier** can appear in the file as many times as the number of the samples in the analysis. This column can be used to pull frequencies of nucleotides for a given nucleotide position from all samples.

3. **sample_id** corresponds to the sample name a given particular line is reported from. This column allows the linking of SNVs and the sample(s) they were identified from.

4. **pos** refers to the nucleotide position in the split.

5. **pos_in_contig** refers to the nucleotide position in the contig (why is this called *pos_in_contig*, and the one before is not called *pos_in_split*? Well, we have been wondering about that for a long time, too).

6. **corresponding_gene_call** refers to a unique gene caller id (`-1`, if the position falls out of a gene call).

7. **in_partial_gene_call** indicates whether the gene call is incomplete (i.e., starts with a start codon, stops with a stop codon, etc). `1` if incomplete, `0` if both start and stop positions are detected, or if the position is not in a gene.

8. **in_complete_gene_call** indicates the gene completion status. `1` for complete, `0` if incomplete, or if the position is not in a gene.

9. **base_pos_in_codon** refers to the position of the nucleotide in a codon. `1`, `2` or `3` for codon positions, `-1` if the position is not in a detected gene.

10. **codon_order_in_gene** refers to the order of the codon in the gene call, starting from the start position. `-1` if the position is not in a called gene.

11. **coverage** refers to the coverage, the number of recruited reads mapping to this position.

12. **cov_outlier_in_split** indicates whether the coverage of this position is marked as an outlier compared to all other positions in the split. `1` if outlier, and `0` if not.

13. **cov_outlier_in_contig** has the same purpose with the one above, except it is at the contig-level.

14. **departure_from_reference** refers to the ratio of nucleotides in a given position that diverge from the reference nucleotide. If a position with a coverage of 100X has the nucleotide `A` in the reference, the departure from reference would be `0.2` if the frequency of mapping nucleotides are as follows: `A: 80X`, `T: 0X`, `C: 12X`, and `G: 8X`. Note that the departure from reference can dramatically change across samples, revealing subtle differences at the single nucleotide level.

15. **competing_nts** refers to the two most represented nucleotides. Note that if competing nucleotides are stored as `CT`, it does not necessarily mean that `C` occurs more than `T`, since they are ordered alphabetically. For positions positions that does not have any variation, competing nucleotides will contain two identical nucleotides. I hope here you are asking yourself here "*why would a position without any variation would appear in this table?*". The answer is `--quince-mode`. See below.

16. **reference** refers to the reference nucleotide in the mapping context.

17. **A** refers to the number of mapped reads covering this position with a `A`.

18. **T** refers to the number of mapped reads covering this position with a `T`.

19. **C** refers to the number of mapped reads covering this position with a `C`.

20. **G** refers to the number of mapped reads covering this position with a `G`.

21. **N** refers to the number of mapped reads covering this position with an ambiguous base.

22. **consensus** refers to the most frequent nucleotide mapping to this position. This column is used to define the **departure_from_consensus** ratio.

23. **departure_from_consensus** refers to the ratio of nucleotides in a given position that diverge from the most frequent nucleotide. If the nucleotide frequencies for a given position are, `A: 10X`, `T: 20X`, `C: 30X`, and `G: 70X`, then the departure from reference would be `0.46` (from `(10X + 20X + 30X) / (10X + 20X + 30X + 70X)`) if the frequency of mapping nucleotides for this position is . Compared to the **departure from reference**, this column relies less on the reference sequence and thus might better reflect the variation in the environment regardless of the context you are mapping against.

24. **n2n1ratio** refers to the ratio of the second most frequent nucleotide to the consensus nucleotide. If the frequency of nucleotides mapped to this position are `A: 10X`, `T: 20X`, `C: 30X`, and `G: 70X`, then the n2n1ratio would be `0.42`, .. and you would be as close as two orders of magnitude to the answer to the ultimate question of life, the universe, and everything (which, in an ideal world, should include an answer to your research question, too, you lucky you (but then is there really an answer to your research question? BAM! (If you survived this, you can survive the rest of this tutorial. Please carry on))).

25. **contig_name** refers to the contig name as it appears in the contigs database.

26. **split_name** refers to the split name.

So this is the output format. The following section will detail filters to scrutinize it.

## Parameters to filter the output

This program allows the researcher to specify _filters_ and refine reported variable nucleotide positions with a rich array of parameters:

- With `anvi-gen-variability-profile` you can request all variable positions identified in a genome bin stored in a collection using the `--collection-name` and `--bin-id` parameters. Or alternatively you can focus on an arbitrary list of splits using `--splits-of-interest` parameter. You can also focus on a specific gene calls using the `--genes-of-interest` parameter.

- By default, the program will report SNVs from every sample found in the profile database. You can use the `--samples-of-interest` parameter to define which samples should be reported.

- The output file will not have columns to display in which split or contig a given SNV was identified. But you can change that behavior using the `--include-contig-names` and/or `--include-split-names` flags.

- You can use the `--num-positions-from-each-split` to define a maximum number of SNVs to be reported from each split. If a given split has more SNVs than this number, the program will randomly select a matching number of them.

- Using the `--min-scatter` parameter you can eliminate some SNV positions based on how they partition samples. This one is a bit tricky, but Meren wants to keep it in the code base. If you skip this you will not lose anything, but for the nerd kind, this is how it goes: If you have `N` samples in your dataset, a given variable position `x` in one of your splits can split your `N` samples into `t` groups based on the identity of the variation they harbor at position `x`. For instance, `t` would have been 1, if all samples had the same type of variation at position `x` (which would not be very interesting, because in this case position `x `would have zero contribution to a deeper understanding of how these samples differ based on variability. When `t` > 1, it would mean that identities at position `x` across samples *do* differ. But how much scattering occurs based on position `x` when `t > 1`? If `t=2`, how many samples would end up in each group? Obviously, the even distribution of samples across groups may tell us something different than uneven distribution of samples across groups. So, this parameter filters out any `x` if 'the number of samples in the second largest group' (=scatter) is less than the value you choose. Here is an example: lets assume you have 7 samples. While 5 of those have `AG`, 2 of them have `TC` at position x. This would mean the scatter of `x` is 2. If you set `-m` to 2, this position would not be reported in your output matrix. The default value for `-m` is 0, which means every `x` found in the database and survived previous filtering criteria will be reported. Naturally, `-m` can not be more than half of the number of samples.

- You can set a minimum coverage value for each SNV using the parameter `--min-coverage-in-each-sample`. Any SNV that is covered less than this value would be omitted from the output.

- You can filter nucleotide positions based upon frequencies using parameters `--min-departure-from-reference`, `--max-departure-from-reference`, `--min-departure-from-consensus`, or `--max-departure-from-consensus`.

- You can remove nucleotide positions from the final report that show variation in less than any number of samples using the `--min-occurrence` parameter.

- One of the very important flags is the flag `--quince-mode`. If you read the previous section you know that the base frequencies are reported for a given nucleotide position in a sample *only* if that position has variation. But it doesn't work well when you want to be able to say "give me back SNVs only if the nucleotide position they occur into is covered more than X in every sample". But not every sample will be in the variability profile, hence, the program will not have access to the coverage scores for nucleotide positions for samples with 0 variation at those positions. When you use the `--quince-mode` flag, the program goes back to the auxiliary files generated during profiling, and recovers coverage values for nucleotide positions in sampels for each **unique_pos_identifier**. In other words, if there are 10 samples in a dataset, and a given position has been reported as a variable site during the profiling in only one of those samples, there will be no information stored in the database for the remaining 9. When this flag is used, base frequencies will be recovered. It will take considerably longer to report when this flag is on, and the use of it will increase the file size dramatically, however, it may be necessary to use for some statistical approaches (as well as for some beautiful visualizations).

{:.notice}
**Fact**: The answer is "*yes*", if you are asking yourself "*is this 'quince' in `--quince-mode` the 'Quince' in '[Christopher Quince](http://www2.warwick.ac.uk/fac/med/staff/cquince/)'?*". We added this flag because he complained about the fact that this functionality was not available in anvi'o at the time.

{:.notice}
**Fact**: The answer is "*very very likely*" if you are now asking yourself "*would I get my own flag if I complain about lacking functionality in anvi'o?*". We the anvi'o developers are often a civil bunch of people.
 Note: if you are asking yourself "*is this 'quince' in `--quince-mode` the 'Quince' in '[Christopher Quince](http://www2.warwick.ac.uk/fac/med/staff/cquince/)'?*", the answer is "yes". We added this flag because he complained about the fact that that functionality was not available at the time. If you are now asking yourself "*would I get my own flag if I complain about lacking functionality?*", the answer to that is "yes", too.
---

We certainly hope SNV characterization framework in anvi'o would be useful to you. Please don't hesitate get in touch with us with questions, or suggestions.

If you are interested, you can find an easy-to-follow tutorial here, demonstrating a SNV analysis on the Infant Gut Dataset: [http://merenlab.org/tutorials/infant-gut/#profiling-snvs-in-a-bin]({{ site.url }}/tutorials/infant-gut/#profiling-snvs-in-a-bin).

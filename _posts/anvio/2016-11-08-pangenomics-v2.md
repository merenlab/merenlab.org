---
layout: post
title: "An anvi'o workflow for microbial pangenomics"
excerpt: "The user-friendly interface anvi'o provides to work with pangenomes."
modified: 2016-11-08
tags: []
categories: [anvio]
comments: true
redirect_from: /2015/11/14/pangenomics/
---

{% include _toc.html %}

{% include _project-anvio-version.html %}

{: .notice}
This pangenomic workflow is for `2.1.0` and later versions of anvi'o. You can learn which version you have on your computer by typing `anvi-profile --version` in your terminal. If your anvi'o installation is `v2.0.2` or earlier, you can follow the pangenomic workflow described [here]({% post_url anvio/2015-11-14-pangenomics-v1 %}) (and not read this one to not be upset about what you are missing).

{:.notice}
{% include _fixthispage.html source="anvio/2016-11-08-pangenomics-v2.md" %}

{% capture images %}{{site.url}}/images/anvio/2016-11-08-pan-genomics{% endcapture %}

With the [anvi'o](https://peerj.com/articles/1319/) pangenomic workflow you can,

* **Identify protein clusters**,
* **Visualize** the distribution of protein clusters across your genomes,
* **Organize** genomes based on shared protein clusters,
* Interactively (or programmatically) **bin protein clusters** into collections, and **summarize your bins**,
* **Annotate** your genes, and **inspect** amino acid alignments within your protein clusters,
* Extend your pangenome with **contextual information** about your genomes or protein clusters,
* **Combine metagenome-assembled genomes** from your anvi'o projects directly **with cultivar genomes** from other sources.

{:.notice}
You can use anvi'o for pangenomic workflow even if you haven't done any metagenomic work with anvi'o. All you need is an anvi'o installation, and a FASTA file for each of your genomes.

## Introduction

The anvi'o pangenomic workflow is composed of three main steps:

* Generating an anvi'o genomes storage using the program `anvi-gen-genomes-storage`.
* Generating an anvi'o pan database using the program `anvi-pan-genome` (which requires a genomes storage)
* Displaying results using the program `anvi-display-pan` (which requires a genomes storage and an anvi'o pan database).

You can then use the interactive interface to bin your protein clusters into collections, or use the program `anvi-import-collection` to import bins for your protein clusters, and finally you can use the program `anvi-summarize` to create a static HTML summary of your results. Easy peasy!

Following sections will detail each step, and finally an example run will follow, but let's first make sure you have all the dependencies installed and your installation is good to go.

### Dependencies

Even if you have a complete installation of anvi'o, the pangenomic workflow uses some additional software that may not be installed on your system.

* DIAMOND or NCBI's blastp for search.
* [MCL]({% post_url anvio/2016-06-18-installing-third-party-software %}#mcl) for clustering.
* [muscle]({% post_url anvio/2016-06-18-installing-third-party-software %}#muscle) for alignment. Which is optional. If you don't have it, amino acid sequences within your protein clusters will not be aligned, and will look ugly.

If your system is properly setup, this command should run without any errors:

{% highlight bash %}
$ anvi-self-test --suite pangenomics
{% endhighlight %}

## Generating an anvi'o genomes storage

An anvi'o *genomes storage* is a special database that stores information about genomes. A genomes storage can be generated only from *external* genomes, only from *internal* genomes, or it can contain both types. Before we go any further, here are some definitions to clarify things:

* **An external genome** is anything you have in a FASTA file format (i.e., a genome you downloaded from NCBI, or obtained through any other way). 

* **An internal genome** is any *genome bin* you stored in an anvi'o collection at the end of your metagenomic analysis (if you are not familiar with the anvi'o metagenomic workflow please take a look at [this post]({% post_url anvio/2016-06-22-anvio-tutorial-v2 %})). 

{:.notice}
**Converting FASTA files into anvi'o contigs databases**: Working with *internal* genomes is quite straightforward since you already have an anvi'o contigs and an anvi'o profile database for them. But if all you have is a bunch of FASTA files, this workflow will require you to convert each of them into an anvi'o contigs database. There is a lot of information about how to create an anvi'o contigs database [here]({% post_url anvio/2016-06-22-anvio-tutorial-v2 %}/#anvi-gen-contigs-database), but if you feel lazy, you can also use the script `anvi-script-FASTA-to-contigs-db`, which takes a single parameter: the FASTA file path. Power users should consider taking a [look at the code](https://github.com/meren/anvio/blob/master/sandbox/anvi-script-FASTA-to-contigs-db#L50), and create their own batch script with improvements on those lines based on their needs (for instance, increasing the number of threads when running HMMs, etc). Also, you may want to run `anvi-run-ncbi-cogs` on your contigs database to annotate your genes.

You can create a new anvi'o genomes storage using the program `anvi-gen-genomes-storage`, which will require you to provide descriptions of genomes to be included in this storage. File formats for external genome and internal genome descriptions differ slightly. For instance, this is an example `--external-genomes` file:

|name|contigs_db_path|
|:--|:--|
|Name_01|/path/to/contigs-01.db|
|Name_02|/path/to/contigs-02.db|
|Name_03|/path/to/contigs-03.db|
|(...)|(...)|

and this is an example file for `--internal-genomes`:

|name|bin_id|collection_id|profile_db_path|contigs_db_path|
|:--|:--:|:--:|:--|:--|
|Name_01|Bin_id_01|Collection_A|/path/to/profile.db|/path/to/contigs.db|
|Name_02|Bin_id_02|Collection_A|/path/to/profile.db|/path/to/contigs.db|
|Name_03|Bin_id_03|Collection_B|/path/to/another_profile.db|/path/to/another/contigs.db|
|(...)|(...)|(...)|(...)|(...)|

{:.notice}
For names in the first column please use only letters, digits, and the underscore character.

Thanks to these two files, genome bins in anvi'o collections and cultivar genomes could be combined and analyzed together seamlessly. The most essential need for the coherence of the genomes storage is to make sure each internal and external genome is generated identically with respect to how genes were called, how functions were assigned, etc. Anvi'o will check for some things, but it can't stop you from doing terrible mistakes. For instance, if the version of Prodigal is not identical across all contigs databases, the genes described in genomes storage will not necessarily be comparable. If you are not sure about something, send us an e-mail, and we will be happy to try to clarify.

## Running a pangenome analysis

Once you have your genomes storage ready, you can use the program `anvi-pan-genome` to run the actual pangenomic analysis. This is the simplest form of this command:

{% highlight bash %}
$ anvi-pan-genome -g MY-GENOMES.h5 -n PROJECT_NAME
{% endhighlight %}

When you run it, `anvi-pan-genome`,

* Will use all genomes in the genomes storage. If you would like to focus on a subset, you can use the parameter `--genome-names`.

* Will use only a single core by default. Depending on the number of genomes you are analyzing, this process can be very time consuming, hence you should consider increasing the number of threads to be used via the parameter `--num-threads`.

* Will use [DIAMOND](http://ab.inf.uni-tuebingen.de/software/diamond/) ([Buchnfink et al., 2015](http://www.nature.com/nmeth/journal/v12/n1/abs/nmeth.3176.html)) in 'fast' mode (you can ask DIAMOND to be 'sensitive' by using the flag `--sensitive`) by default to calculate similarities of each protein in every genome against every other protein (which clearly requires you to have DIAMOND installed). Alternatively you could use the flag `--use-ncbi-blast` to use NCBI's `blastp` for protein search.

{:.notice}
***A note from Meren***: I strongly suggest you to do your analysis with the `--use-ncbi-blast` flag. Yes, DIAMOND is very fast, and it may take 8-9 hours to analyze 20-30 genomes with `blastp`. But dramatic increases in speed *rarely* comes without major trade-offs in sensitivity and accuracy, and some of my observations tell me that DIAMOND is not one of those *rare* instances. This clearly deserves a more elaborate discussion, and maybe I will have a chance to write it later, but for now take this as a friendly reminder.

* Will use every gene call, whether they are complete or not. Although this is not a big concern for complete genomes, metagenome-assembled genomes (MAGs) will have many incomplete gene calls at the end and at the beginning of contigs. Our experiments so far suggest that they do not cause major issues, but if you want to exclude them, you can use the `--exclude-partial-gene-calls` flag.

* Will use the *minbit heuristic* that was originally implemented in ITEP ([Benedict et al, 2014](http://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-8)) to eliminate weak matches between two protein sequences. You see, the pangenomic workflow first identifies proteins that are somewhat similar by doing similarity searches, and then resolves clusters based on those similarities. In this scenario, weak similarities can connect protein clusters that should not be connected. Although the network partitioning algorithm can recover from these weak connections, it is always better to eliminate as many noise as possible at every step. So the minbit heuristic provides a mean to set a to eliminate weak matches between two protein sequences. We learned it from ITEP (Benedict MN et al, doi:10.1186/1471-2164-15-8), which is another comprehensive analysis workflow for pangenomes, and decided to use it because it makes a lot of sense. Briefly, If you have two protein sequences, `A` and `B`, the minbit is defined as `BITSCORE(A, B) / MIN(BITSCORE(A, A), BITSCORE(B, B))`. So the minbit score between two sequences goes to `1.0` if they are very similar over the entire length of the 'shorter' protein sequence, and goes to `0.0` if (1) they match over a very short stretch compared even to the length of the shorter protein sequence or (2) the match betwen sequence identity is low. The default minbit is `0.5`, but you can change it using the parameter `--minbit`.

{:.notice}
This parameter has wrongly named as `maxbit` in anvi'o `v2.4.0` and earlier versions. [Please take a look at the relevant issue](https://github.com/merenlab/anvio/issues/581) if you are using an anvi'o version that is impacted by this typo.

* Will use the [MCL](http://micans.org/mcl/) algorithm ([van Dongen and Abreu-Goodger, 2012](http://www.ncbi.nlm.nih.gov/pubmed/22144159)) to identify clusters in protein similarity search results. We use `2` as the *MCL inflation parameter* by default. This parameter defines the sensitivity of the algorithm during the identification of the protein clusters. More sensitivity means more clusters, but of course more clusters does not mean better inference of evolutionary relationships. More information on this parameter and it's effect on cluster granularity is here [http://micans.org/mcl/man/mclfaq.html#faq7.2](http://micans.org/mcl/man/mclfaq.html#faq7.2), but clearly, we, the metagenomics people will need to talk much more about this. So far in the Meren Lab we have been using `2` if we are comparing many distantly related genomes (i.e., genomes classify into different families or farther), and `10` if we are comparing very closely related genomes (i.e., 'strains' of the same 'species' (whatever they mean to you)). You can change it using the parameter `--mcl-inflation`. Please experiment yourself, and consider reporting!

* Will utilize every protein cluster, even if they occur in only one genome in your analyses. Of course, the importance of singletons or doubletons will depend on the number of genomes in your analysis, or the question you have in mind. However, if you would like to define a cut-off, you can use the parameter `--min-occurrence`, which is 1, by default. Increasing this cut-off will improve the clustering speed and make the visualization much more manageable, but again, this parameter should be considered in the context of each study.

* Will use `euclidean` distance and `ward` linkage to organize protein clusters and genomes. You can change those using `--distance` and `--linkage` parameters.

* Will try to utilize previous search results if there is already a directory. This way you can play with the `--minbit`, `--mcl-inflation`, or `--min-occurrence` parameters without having to re-do the protein search.  However, if you have changed something, either you need to remove the output directory, or use the `--overwrite-output-destinations` flag to redo the search.

{:.notice}
You need another parameter? Well, of course you do! Let us know, and let's have a discussion. We love parameters.

Once you are done, a new directory with your analysis results will appear. Among other files, it will contain an anvi'o pan database, and an anvi'o [samples database]({% post_url anvio/2015-11-10-samples-db %}). If you would like to update the samples database, you can edit the `*-samples-information.txt` or `*-samples-order.txt` files and re-create an updated samples database using the program `anvi-gen-samples-database` as explained [here]({% post_url anvio/2015-11-10-samples-db %}/#the-command-line-to-generate-the-database). 

## Displaying and exploring the pan genome

Once your analysis is done, you will use the program `anvi-display-pan` to display your results.

This is the simplest form of this command:

{% highlight bash %}
$ anvi-display-pan -p PROJECT-PAN.db -s PROJECT-PAN-SAMPLES.db -g PROJECT-PAN-GENOMES.h5
{% endhighlight %}

The program `anvi-display-pan` is very similar to the program `anvi-interactive`, and the interface that will welcome you is nothing but the standard [anvi'o interactive interface]({% post_url anvio/2016-02-27-the-anvio-interactive-interface %}/#using-the-anvio-interactive-interface) with slight adjustments for pangenomic analyses. Of course `anvi-display-pan` will allow you to set the IP address and port number to serve, add additional view data, additional layers, and/or additional trees, and more. Please familiarize yourself with it by running `anvi-display-pan -h` in your terminal.

Here is an example to show how the pangenomic analysis of 31 *Prochlorococcus* genomes we just recently finished looks like when we run `anvi-display-pan` the very first time on it:

[![31 Prochlorococcus raw]({{images}}/prochlorococcus-pangenomics-raw.png)]({{images}}/prochlorococcus-pangenomics-raw.png){:.center-img .width-60}

Looks ugly. But do not despair. For instance, to improve things a little, you may like to organize your genomes based on protein clustering results by selecting the 'PC frequencies' tree from the Samples Tab > Sample Order menu:

[![31 Prochlorococcus samples]({{images}}/prochlorococcus-pangenomics-samples-tab.png)]({{images}}/prochlorococcus-pangenomics-samples-tab.png){:.center-img .width-80}

This is what happens when you draw it again (note the tree that appears on the right):

[![31 Prochlorococcus ordered]({{images}}/prochlorococcus-pangenomics-ordered.png)]({{images}}/prochlorococcus-pangenomics-ordered.png){:.center-img .width-60}

Looks more meaningful .. but still ugly.

Well, this is exactly where you need to start using the interface more efficiently. For instance, this is the same pangenome after we are done with it:

[![31 Prochlorococcus final]({{images}}/prochlorococcus-pangenomics-final.png)]({{images}}/prochlorococcus-pangenomics-final.png){:.center-img .width-60}

No excuses for bad looking pangenomes.

Every protein cluster in your analysis will contain one or more amino acid sequences that originate from one or more genomes. While there will likely be a 'core' section, in which all protein cluster will appear in every genome, it is also common to find protein clusters that contain more than one gene call from a single genome (i.e., all multi-copy genes in a given genome will end up in the same protein cluster). Sooner or later you will start getting curious about some of the protein clusters, and want to learn more about them. Luckily you can right-click on to any protein cluster, and you would see this menu (or maybe even more depending on when you are reading this article):

[![31 Prochlorococcus final]({{images}}/pc-right-click.png)]({{images}}/pc-right-click.png){:.center-img .width-80}

For instance, if you click 'Inspect PC', you will see all the amino acid sequences from each genome that went into that protein cluster (with the same order and background colors of genomes as they are arranged in the main display):

[![31 Prochlorococcus final]({{images}}/pc-inspect.png)]({{images}}/pc-inspect.png){:.center-img .width-60}

It is not only fun but also very educational to go through protein clusters one by one. Fine. But what do you do if you want to make sense of large selections?

As you already know, the anvi'o interactive interface allows you to make selections from the tree. So you can select groups of protein clusters into bins (don't mind the numbers on the left panel, there clearly is a bug, and will be fixed in your version):

[![31 Prochlorococcus selection]({{images}}/pc-selection.gif)]({{images}}/pc-selection.gif){:.center-img .width-80}

You can create multiple bins with multiple selections, and even give them meaningful names if you fancy:

[![31 Prochlorococcus collection]({{images}}/pc-collection.png)]({{images}}/pc-collection.png){:.center-img .width-60}

Then you can save your selection as a "collection", and then run the `anvi-summarize` to generate a summary. Which would be as simple as this:

{% highlight bash %}
$ anvi-summarize -p PROJECT-PAN.db -g PROJECT-PAN-GENOMES.h5 -C COLLECTION_NAME -o PROJECT-SUMMARY
{% endhighlight %}

If you open the `index.html` file in the summary directory, you will see an output with some essential information about the analysis:

[![31 summary]({{images}}/summary.png)]({{images}}/summary.png){:.center-img .width-60}

As well as a TAB-delimited file that describes every protein cluster:

[![31 summary]({{images}}/summary-file.png)]({{images}}/summary-file.png){:.center-img .width-60}

The structure of this file will look like this, and will give you an opportunity to investigate your protein clusters in a more detailed manner (you may need to scroll towards right to see more of the table):

|unique_id|protein_cluster_id|bin_name|genome_name|gene_callers_id|COG_FUNCTION_ACC|COG_FUNCTION|COG_CATEGORY_ACC|COG_CATEGORY|aa_sequence|
|:--:|:--|:--|:--|:--|:--|:--|:--|:--|:--|
|1|PC_00001990||MIT9303|30298|COG1199|Rad3-related DNA helicase|L|Replication, recombination and repair|MLEARSHQQLKHLLLQNSSP(...)|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|
|91|PC_00001434|Core_HL|MIT9322|42504|COG3769|Predicted mannosyl-3-phosphoglycerate phosphatase|G|Carbohydrate transport and metabolism|MIENSSIWVVSDVDGTLMDH(...)|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|
|257|PC_00000645|Core_all|MIT9322|42217|COG1185|Polyribonucleotide nucleotidyltransferase|J|Translation, ribosomal structure and biogenesis|MEGQNKSITFDGREIRLTTG(...)|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|
|2019|PC_00001754|Core_LL|NATL2A|49129|COG0127|Inosine/xanthosine triphosphate pyrophosphatase|F|Nucleotide transport and metabolism|MDNVPLVIASGNKGKIGEFK(...)|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|
|5046|PC_00001653||PAC1|52600|COG1087|UDP-glucose 4-epimerase|M|Cell wall/membrane/envelope biogenesis|MRVLLTGGAGFIGSHIALLL(...)|
|5047|PC_00001653||LG|7488|COG1087|UDP-glucose 4-epimerase|M|Cell wall/membrane/envelope biogenesis|MNRILVTGGAGFIGSHTCIT(...)|
|5048|PC_00001653||SS35|56661|COG1087|UDP-glucose 4-epimerase|M|Cell wall/membrane/envelope biogenesis|MNRILVTGGAGFIGSHTCIT(...)|
|5049|PC_00001653||NATL2A|49604|COG1087|UDP-glucose 4-epimerase|M|Cell wall/membrane/envelope biogenesis|MRVLLTGGSGFIGSHVALLL(...)|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|

The summary directory is a self-contained, static HTML output you can send to your colleagues, or add to your publications as a supplementary dataset. 


## Final words

We are realizing that there is a lot to explore in this front, and excited to work with the community. Please let us know if you have any suggestions using our [discussion group]({% post_url anvio/2015-10-14-anvio-discussion-group %}), the comments section below. You can also keep an eye on our [public code repository](http://github.com/meren/anvio) for new releases.

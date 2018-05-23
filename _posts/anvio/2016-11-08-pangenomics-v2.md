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
With changed all instances of the term '**protein cluster**' in our pangenomic workflow with '**gene cluster**' in anvi'o `v4`. The reason behind this change is detailed [here](https://github.com/merenlab/anvio/issues/644).

{:.notice}
{% include _fixthispage.html source="anvio/2016-11-08-pangenomics-v2.md" %}

{% capture images %}{{site.url}}/images/anvio/2016-11-08-pan-genomics{% endcapture %}

With the [anvi'o](https://peerj.com/articles/1319/) pangenomic workflow you can,

* **Identify gene clusters** among your genomes,
* **Visualize** the distribution of gene clusters across your genomes,
* **Estimate relationships** between your genomes based on shared gene clusters,
* Interactively (or programmatically) **bin gene clusters** into collections, and **summarize your bins**,
* **Annotate** your genes, and **inspect** amino acid alignments within your gene clusters,
* Extend your pangenome with **contextual information** about your genomes or gene clusters,
* **Combine metagenome-assembled genomes** from your anvi'o projects directly **with cultivar genomes** from other sources.

{:.notice}
You can use anvi'o for pangenomic workflow even if you haven't done any metagenomic work with anvi'o. All you need is an anvi'o installation, and a FASTA file for each of your genomes.

## Introduction

The anvi'o pangenomic workflow is composed of three main steps:

* Generating an anvi'o genomes storage using the program `anvi-gen-genomes-storage`.
* Generating an anvi'o pan database using the program `anvi-pan-genome` (which requires a genomes storage)
* Displaying results using the program `anvi-display-pan` (which requires a genomes storage and an anvi'o pan database).

You can then use the interactive interface to bin your gene clusters into collections, or use the program `anvi-import-collection` to import bins for your gene clusters, and finally you can use the program `anvi-summarize` to create a static HTML summary of your results. Easy peasy!

Following sections will detail each step, and finally an example run will follow, but let's first make sure you have all the dependencies installed and your installation is good to go.

### Dependencies

Even if you have a complete installation of anvi'o, the pangenomic workflow uses some additional software that may not be installed on your system.

* DIAMOND or NCBI's blastp for search.
* [MCL]({% post_url anvio/2016-06-18-installing-third-party-software %}#mcl) for clustering.
* [muscle]({% post_url anvio/2016-06-18-installing-third-party-software %}#muscle) for alignment. Which is optional. If you don't have it, amino acid sequences within your gene clusters will not be aligned, and will look ugly.

If your system is properly setup, this command should run without any errors:

``` bash
$ anvi-self-test --suite pangenomics
```

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

<div class="extra-info" markdown="1">

<span class="extra-info-header">A real example for hardcore tutorial readers: A *Prochlorococcus* pangenome</span>

For the sake of reproducibility, the rest of the tutorial will follow a real example.

We will simply create a pangenome of 31 [Prochlorococcus isolate genomes that we used in this study](https://peerj.com/articles/4320/).

If you wish to follow the tutorial on your computer, you can download the Prochlorococcus data pack ([doi:10.6084/m9.figshare.6318833](https://doi.org/10.6084/m9.figshare.6318833)) which contains anvi'o contigs databases for these isolate genomes on your computer:

``` bash
wget https://ndownloader.figshare.com/files/11568539 -O Prochlorococcus_31_genomes.tar.gz
tar -zxvf Prochlorococcus_31_genomes.tar.gz
cd Prochlorococcus_31_genomes
```

The directory contains anvi'o contigs databases, an external genomes file, and a layer additional data file. You can generate a genomes storage as described in this section the following way:

``` bash
anvi-gen-genomes-storage -e external-genomes.txt \
                         -o PROCHLORO-GENOMES.db
```

</div>

## Running a pangenome analysis

Once you have your genomes storage ready, you can use the program `anvi-pan-genome` to run the actual pangenomic analysis. This is the simplest form of this command:

``` bash
$ anvi-pan-genome -g MY-GENOMES.db -n PROJECT_NAME
```

When you run it, `anvi-pan-genome`,

* Will use all genomes in the genomes storage. If you would like to focus on a subset, you can use the parameter `--genome-names`.

* Will use only a single core by default. Depending on the number of genomes you are analyzing, this process can be very time consuming, hence you should consider increasing the number of threads to be used via the parameter `--num-threads`.

* Will use [DIAMOND](http://ab.inf.uni-tuebingen.de/software/diamond/) ([Buchnfink et al., 2015](http://www.nature.com/nmeth/journal/v12/n1/abs/nmeth.3176.html)) in 'fast' mode (you can ask DIAMOND to be 'sensitive' by using the flag `--sensitive`) by default to calculate similarities of each protein in every genome against every other protein (which clearly requires you to have DIAMOND installed). Alternatively you could use the flag `--use-ncbi-blast` to use NCBI's `blastp` for protein search.

{:.notice}
***A note from Meren***: I strongly suggest you to do your analysis with the `--use-ncbi-blast` flag. Yes, DIAMOND is very fast, and it may take 8-9 hours to analyze 20-30 genomes with `blastp` (using one core on a standard laptop). But dramatic increases in speed *rarely* comes without major trade-offs in sensitivity and accuracy, and some of my observations tell me that DIAMOND is not one of those *rare* instances. This clearly deserves a more elaborate discussion, and maybe I will have a chance to write it later, but for now take this as a friendly reminder.

* Will use every gene call, whether they are complete or not. Although this is not a big concern for complete genomes, metagenome-assembled genomes (MAGs) will have many incomplete gene calls at the end and at the beginning of contigs. Our experiments so far suggest that they do not cause major issues, but if you want to exclude them, you can use the `--exclude-partial-gene-calls` flag.

* Will use the *minbit heuristic* that was originally implemented in ITEP ([Benedict et al, 2014](http://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-8)) to eliminate weak matches between two protein sequences. You see, the pangenomic workflow first identifies proteins that are somewhat similar by doing similarity searches, and then resolves clusters based on those similarities. In this scenario, weak similarities can connect gene clusters that should not be connected. Although the network partitioning algorithm can recover from these weak connections, it is always better to eliminate as many noise as possible at every step. So the minbit heuristic provides a mean to set a to eliminate weak matches between two protein sequences. We learned it from ITEP (Benedict MN et al, doi:10.1186/1471-2164-15-8), which is another comprehensive analysis workflow for pangenomes, and decided to use it because it makes a lot of sense. Briefly, If you have two protein sequences, `A` and `B`, the minbit is defined as `BITSCORE(A, B) / MIN(BITSCORE(A, A), BITSCORE(B, B))`. So the minbit score between two sequences goes to `1.0` if they are very similar over the entire length of the 'shorter' protein sequence, and goes to `0.0` if (1) they match over a very short stretch compared even to the length of the shorter protein sequence or (2) the match betwen sequence identity is low. The default minbit is `0.5`, but you can change it using the parameter `--minbit`.

{:.notice}
This parameter has wrongly named as `maxbit` in anvi'o `v2.4.0` and earlier versions. [Please take a look at the relevant issue](https://github.com/merenlab/anvio/issues/581) if you are using an anvi'o version that is impacted by this typo.

* Will use the [MCL](http://micans.org/mcl/) algorithm ([van Dongen and Abreu-Goodger, 2012](http://www.ncbi.nlm.nih.gov/pubmed/22144159)) to identify clusters in protein similarity search results. We use `2` as the *MCL inflation parameter* by default. This parameter defines the sensitivity of the algorithm during the identification of the gene clusters. More sensitivity means more clusters, but of course more clusters does not mean better inference of evolutionary relationships. More information on this parameter and it's effect on cluster granularity is here [http://micans.org/mcl/man/mclfaq.html#faq7.2](http://micans.org/mcl/man/mclfaq.html#faq7.2), but clearly, we, the metagenomics people will need to talk much more about this. So far in the Meren Lab we have been using `2` if we are comparing many distantly related genomes (i.e., genomes classify into different families or farther), and `10` if we are comparing very closely related genomes (i.e., 'strains' of the same 'species' (whatever they mean to you)). You can change it using the parameter `--mcl-inflation`. Please experiment yourself, and consider reporting!

* Will utilize every gene cluster, even if they occur in only one genome in your analyses. Of course, the importance of singletons or doubletons will depend on the number of genomes in your analysis, or the question you have in mind. However, if you would like to define a cut-off, you can use the parameter `--min-occurrence`, which is 1, by default. Increasing this cut-off will improve the clustering speed and make the visualization much more manageable, but again, this parameter should be considered in the context of each study.

* Will use `euclidean` distance and `ward` linkage to organize gene clusters and genomes. You can change those using `--distance` and `--linkage` parameters.

* Will try to utilize previous search results if there is already a directory. This way you can play with the `--minbit`, `--mcl-inflation`, or `--min-occurrence` parameters without having to re-do the protein search.  However, if you have changed something, either you need to remove the output directory, or use the `--overwrite-output-destinations` flag to redo the search.

{:.notice}
You need another parameter? Well, of course you do! Let us know, and let's have a discussion. We love parameters.

Once you are done, a new directory with your analysis results will appear. You can add or remove additional data items into your pan profile database using anvi'o [additional data tables subsystem]({% post_url anvio/2017-12-11-additional-data-tables %}).

<div class="extra-info" markdown="1">

<span class="extra-info-header">A *Prochlorococcus* pangenome [*continued*]</span>

Let's run the pangenome using the genomes storage we created using the 31 Prochlorococcus isolates:

``` bash
anvi-pan-genome -g PROCHLORO-GENOMES.db \
                --project-name "Prochlorococcus_Pan" \
                --output-dir PROCHLORO \
                --num-threads 6 \
                --minbit 0.5 \
                --mcl-inflation 10 \
                --use-ncbi-blast
```

Every parameter after the `--project-name` is optional (yet matches to the way we run the pangenome for our publication).

The directory you have downloaded also contains a file called "layer-additional-data.txt", which summarizes what clade these genomes belong to. Once the pangenome is computed, we can add it into the pan database:

``` bash
anvi-import-misc-data layer-additional-data.txt \
                      -p PROCHLORO/Prochlorococcus_Pan-PAN.db \
                      --target-data-table layers

New layers additional data...
===============================================
Data key "clade" .............................: Predicted type: str
Data key "light" .............................: Predicted type: str

New data added to the db for your layers .....: clade, light.
```

You can learn more about the items and layers additional data tables [here](http://merenlab.org/2017/12/11/additional-data-tables/).

</div>

## Displaying the pan genome

Once your analysis is done, you will use the program `anvi-display-pan` to display your results.

This is the simplest form of this command:

``` bash
$ anvi-display-pan -p PROJECT-PAN.db -g PROJECT-PAN-GENOMES.db
```

The program `anvi-display-pan` is very similar to the program `anvi-interactive`, and the interface that will welcome you is nothing but the standard [anvi'o interactive interface]({% post_url anvio/2016-02-27-the-anvio-interactive-interface %}/#using-the-anvio-interactive-interface) with slight adjustments for pangenomic analyses. Of course `anvi-display-pan` will allow you to set the IP address and port number to serve, add additional view data, additional layers, and/or additional trees, and more. Please familiarize yourself with it by running `anvi-display-pan -h` in your terminal.

Here is the pangenome for the 31 *Prochlorococcus* isolate genomes we have created in the previous sections of this tutorial:

``` bash
anvi-display-pan -g PROCHLORO-GENOMES.db \
                 -p PROCHLORO/Prochlorococcus_Pan-PAN.db
```

[![31 Prochlorococcus raw]({{images}}/prochlorococcus-pangenomics-raw.png)]({{images}}/prochlorococcus-pangenomics-raw.png){:.center-img .width-60}

Looks ugly. But do not despair. For instance, to improve things a little, you may like to organize your genomes based on gene clustering results by selecting the 'gene_cluster frequencies' tree from the Samples Tab > Sample Order menu:

[![31 Prochlorococcus samples]({{images}}/prochlorococcus-pangenomics-samples-tab.png)]({{images}}/prochlorococcus-pangenomics-samples-tab.png){:.center-img .width-50}

This is what happens when you draw it again (note the tree that appears on the right):

[![31 Prochlorococcus ordered]({{images}}/prochlorococcus-pangenomics-ordered.png)]({{images}}/prochlorococcus-pangenomics-ordered.png){:.center-img .width-60}

Looks more meaningful .. but still ugly.

Well, this is exactly where you need to start using the interface more efficiently. For instance, this is the same pangenome after some changes using the additional settings items in the settings panel of the interactive interface:

[![31 Prochlorococcus final]({{images}}/prochlorococcus-pangenomics-final.png)]({{images}}/prochlorococcus-pangenomics-final.png){:.center-img .width-60}

The anvi'o state file for this display is also in your work directory if you have downloaded the Prochlorococcus data pack, and you can import it this way:

``` bash
anvi-import-state -p PROCHLORO/Prochlorococcus_Pan-PAN.db \
                  --state pan-state.json \
                  --name default
```

No excuses for bad looking pangenomes.

## Inspecting gene clusters

Every gene cluster in your analysis will contain one or more amino acid sequence that originate from one or more genomes. While there will likely be a 'core' section, in which all gene cluster will appear in every genome, it is also common to find gene clusters that contain more than one gene call from a single genome (i.e., all multi-copy genes in a given genome will end up in the same gene cluster). Sooner or later you will start getting curious about some of the gene clusters, and want to learn more about them. Luckily you can right-click on to any gene cluster, and you would see this menu (or maybe even more depending on when you are reading this article):

[![31 Prochlorococcus final]({{images}}/pc-right-click.png)]({{images}}/pc-right-click.png){:.center-img .width-80}

For instance, if you click 'Inspect gene cluster', you will see all the amino acid sequences from each genome that went into that gene cluster (with the same order and background colors of genomes as they are arranged in the main display):

[![31 Prochlorococcus final]({{images}}/pc-inspect.png)]({{images}}/pc-inspect.png){:.center-img .width-60}

{:.notice}
You can learn more about the amino acid color coding algorithm [here](http://merenlab.org/2018/02/13/color-coding-aa-alignments/).

It is not only fun but also very educational to go through gene clusters one by one. Fine. But what do you do if you want to make sense of large selections?

As you already know, the anvi'o interactive interface allows you to make selections from the tree. So you can select groups of gene clusters into bins (don't mind the numbers on the left panel, there clearly is a bug, and will be fixed in your version):

[![31 Prochlorococcus selection]({{images}}/pc-selection.gif)]({{images}}/pc-selection.gif){:.center-img .width-80}

You can create multiple bins with multiple selections, and even give them meaningful names if you fancy:

[![31 Prochlorococcus collection]({{images}}/pc-collection.png)]({{images}}/pc-collection.png){:.center-img .width-60}

While selecting gene clusters manually using the dendrogram is an option, it is also possible to identify them using the search interface that allows you to define very specific search criteria:

[![31 Prochlorococcus collection]({{images}}/search-gene-clusters.png)]({{images}}/search-gene-clusters.png){:.center-img .width-60}

You can highlight these selections in the interface, or you can add them to a collection to summarize them later.

In addition, you can search gene clusters also based on functions:

[![31 Prochlorococcus collection]({{images}}/search-functions.png)]({{images}}/search-functions.png){:.center-img .width-60}

Similarly, you can add these gene clusters into collections with whatever name you like, and summarize those collections later.

{:.warning}
Advanced access to gene clusters is also possible through the command line through the program `anvi-get-sequences-for-gene-clusters`. For more information see [this issue](https://github.com/merenlab/anvio/issues/668#issuecomment-354195886) or this [vignette](http://merenlab.org/software/anvio/vignette/#anvi-get-sequences-for-gene-clusters).


## Summarizing an anvi'o pan genome

When you store your selections of gene clusters as a collection, anvi'o will allow you to summarize these results.

{:.notice}
Even if you want to simply summarize everything in the pan gneome without making any selections in the interface, you still will need a collection in the pan database. But luckily, you can use the program `anvi-script-add-default-collection` to add a default collection that contains every gene cluster.

The summary step gives you two important things: a static HTML web page you can compress and share with your colleagues or add into your publication as a supplementary data file, and a copmrehensive TAB-delimited file in the output directory that describes every gene cluster.

You can summarize a collection using the program `anvi-summarize`, and a generic form of this command looks like this:

``` bash
$ anvi-summarize -p PROJECT-PAN.db \
                 -g PROJECT-PAN-GENOMES.db \
                 -C COLLECTION_NAME \
                 -o PROJECT-SUMMARY
```

If you open the `index.html` file in the summary directory, you will see an output with some essential information about the analysis:

[![31 summary]({{images}}/summary.png)]({{images}}/summary.png){:.center-img .width-60}

As well as the TAB-delimited file for gene clusters:

[![31 summary]({{images}}/summary-file.png)]({{images}}/summary-file.png){:.center-img .width-60}

The structure of this file will look like this, and will give you an opportunity to investigate your gene clusters in a more detailed manner (you may need to scroll towards right to see more of the table):

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

This file will link each gene from each genome with every selection you've made and named through the interface or through the program `anvi-import-collection`, and will also give you access to the amino acid seqeunce and function of each gene.


## Final words

We are realizing that there is a lot to explore in this front, and excited to work with the community. Please let us know if you have any suggestions using our [discussion group]({% post_url anvio/2015-10-14-anvio-discussion-group %}), the comments section below. You can also keep an eye on our [public code repository](http://github.com/meren/anvio) for new releases.

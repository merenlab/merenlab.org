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
* **Combine metagenome-assembled genomes** from your anvi'o projects directly **with cultivar genomes** from other sources.
* **Visualize** the distribution of gene clusters across your genomes,
* **Estimate relationships** between your genomes based on shared gene clusters,
* Interactively (or programmatically) **bin gene clusters** into collections, and **summarize your bins**,
* **Annotate** your genes, and **inspect** amino acid alignments within your gene clusters,
* Extend your pangenome with **contextual information** about your genomes or gene clusters,
* Perform **functional enrichment analyses** on groups of genomes in your pangenome,
* Compute and visualize **average nucleotide identity** scores between you genomes.

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
anvi-migrate-db *.db
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

<div class="extra-info" markdown="1">

<span class="extra-info-header">A real example for hardcore tutorial readers: A *Prochlorococcus* pangenome [*continued*]</span>

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

We all are looking for ways to enrich our pangenomic displays, and anvi'o's dditional data tables are excellent ways to do that. Please take a moment to learn more about them [here](http://merenlab.org/2017/12/11/additional-data-tables/).

</div>

When you run `anvi-pan-genome`, the program will,

* Use all genomes in the genomes storage. If you would like to focus on a subset, you can use the parameter `--genome-names`.

* Use only a single core by default. Depending on the number of genomes you are analyzing, this process can be very time consuming, hence you should consider increasing the number of threads to be used via the parameter `--num-threads`.

* Use [DIAMOND](http://ab.inf.uni-tuebingen.de/software/diamond/) ([Buchnfink et al., 2015](http://www.nature.com/nmeth/journal/v12/n1/abs/nmeth.3176.html)) in 'fast' mode (you can ask DIAMOND to be 'sensitive' by using the flag `--sensitive`) by default to calculate similarities of each protein in every genome against every other protein (which clearly requires you to have DIAMOND installed). Alternatively you could use the flag `--use-ncbi-blast` to use NCBI's `blastp` for protein search.

{:.notice}
***A note from Meren***: I strongly suggest you to do your analysis with the `--use-ncbi-blast` flag. Yes, DIAMOND is very fast, and it may take 8-9 hours to analyze 20-30 genomes with `blastp` (using one core on a standard laptop). But dramatic increases in speed *rarely* comes without major trade-offs in sensitivity and accuracy, and some of my observations tell me that DIAMOND is not one of those *rare* instances. This clearly deserves a more elaborate discussion, and maybe I will have a chance to write it later, but for now take this as a friendly reminder.

* Use every gene call, whether they are complete or not. Although this is not a big concern for complete genomes, metagenome-assembled genomes (MAGs) will have many incomplete gene calls at the end and at the beginning of contigs. Our experiments so far suggest that they do not cause major issues, but if you want to exclude them, you can use the `--exclude-partial-gene-calls` flag.

* Use the *minbit heuristic* that was originally implemented in ITEP ([Benedict et al, 2014](http://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-8)) to eliminate weak matches between two protein sequences. You see, the pangenomic workflow first identifies proteins that are somewhat similar by doing similarity searches, and then resolves clusters based on those similarities. In this scenario, weak similarities can connect gene clusters that should not be connected. Although the network partitioning algorithm can recover from these weak connections, it is always better to eliminate as many noise as possible at every step. So the minbit heuristic provides a mean to set a to eliminate weak matches between two protein sequences. We learned it from ITEP (Benedict MN et al, doi:10.1186/1471-2164-15-8), which is another comprehensive analysis workflow for pangenomes, and decided to use it because it makes a lot of sense. Briefly, If you have two protein sequences, `A` and `B`, the minbit is defined as `BITSCORE(A, B) / MIN(BITSCORE(A, A), BITSCORE(B, B))`. So the minbit score between two sequences goes to `1.0` if they are very similar over the entire length of the 'shorter' protein sequence, and goes to `0.0` if (1) they match over a very short stretch compared even to the length of the shorter protein sequence or (2) the match betwen sequence identity is low. The default minbit is `0.5`, but you can change it using the parameter `--minbit`.

{:.notice}
This parameter has wrongly named as `maxbit` in anvi'o `v2.4.0` and earlier versions. [Please take a look at the relevant issue](https://github.com/merenlab/anvio/issues/581) if you are using an anvi'o version that is impacted by this typo.

* Use the [MCL](http://micans.org/mcl/) algorithm ([van Dongen and Abreu-Goodger, 2012](http://www.ncbi.nlm.nih.gov/pubmed/22144159)) to identify clusters in protein similarity search results. We use `2` as the *MCL inflation parameter* by default. This parameter defines the sensitivity of the algorithm during the identification of the gene clusters. More sensitivity means more clusters, but of course more clusters does not mean better inference of evolutionary relationships. More information on this parameter and it's effect on cluster granularity is here [http://micans.org/mcl/man/mclfaq.html#faq7.2](http://micans.org/mcl/man/mclfaq.html#faq7.2), but clearly, we, the metagenomics people will need to talk much more about this. So far in the Meren Lab we have been using `2` if we are comparing many distantly related genomes (i.e., genomes classify into different families or farther), and `10` if we are comparing very closely related genomes (i.e., 'strains' of the same 'species' (whatever they mean to you)). You can change it using the parameter `--mcl-inflation`. Please experiment yourself, and consider reporting!

* Utilize every gene cluster, even if they occur in only one genome in your analyses. Of course, the importance of singletons or doubletons will depend on the number of genomes in your analysis, or the question you have in mind. However, if you would like to define a cut-off, you can use the parameter `--min-occurrence`, which is 1, by default. Increasing this cut-off will improve the clustering speed and make the visualization much more manageable, but again, this parameter should be considered in the context of each study.

* Use `euclidean` distance and `ward` linkage to organize gene clusters and genomes. You can change those using `--distance` and `--linkage` parameters.

* Try to utilize previous search results if there is already a directory. This way you can play with the `--minbit`, `--mcl-inflation`, or `--min-occurrence` parameters without having to re-do the protein search.  However, if you have changed something, either you need to remove the output directory, or use the `--overwrite-output-destinations` flag to redo the search.

{:.notice}
You need another parameter? Well, of course you do! Let us know, and let's have a discussion. We love parameters.

Once you are done, a new directory with your analysis results will appear. You can add or remove additional data items into your pan profile database using anvi'o [additional data tables subsystem]({% post_url anvio/2017-12-11-additional-data-tables %}).

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


## Making some sense of functions in your pangenome

Once we have our pangenome, one of the critical things that we usually want to do is look at the functions associated with our gene clusters. This is a crucial yet complicated challenge with multiple ways to approach. Here, we will describe how you can identify functions that are enriched for some of the clades or sub clades that are included in your pangenome. In addition, we will discuss how you can find the functional core of the pangenome. This is done with our new program `anvi-get-enriched-functions-per-pan-group`.

This program utilizes information in the layers additional data table of your pan database to identify 'groups' within your genomes, and find functions that are enriched in those groups, i.e. functions that are characteristic of these genomes, and predominantly absent from genomes from outside this group. To use this feature you must have at least one categorical additional layer information (which can easily be done via `anvi-import-misc-data`), and at least one functional annotation source for your genomes storage (which will automatically be the case if every contigs database that were used when you run `anvi-gen-genomes-storage` was annotated with the same functional source).

<div class="extra-info" markdown="1">

<span class="extra-info-header">Alon's short 'behind the scenes' lecture on anvi'o functional enrichment analysis for the curious</span>

For those of you who like to dive into the details, here is some information about what goes on behind the scenes.

In order for this analysis to be compatible with everything else relating to the pangenome, we decided that it should be focused on gene clusters, since they are the heart of our pangenomes. This means that the first step of the functional analysis is to try to associate every gene cluster with a function. Because gene clusters do not have functional annotation by default as functional annotation is done in anvi'o at the level of individual genes.

In an ideal world, all genes in a single gene cluster would be annotated with the same identical function. While in reality this is usually the case, it is not always true. In cases in which there are multiple functions associated with a gene cluster, we chose the most frequent function, i.e. the one that the largest number of genes in the gene cluster matches (if there is a tie of multiple functions, then we simply choose one arbitrarily). If none of the genes in the gene cluster were annotated, then the gene cluster has no function associated to it.

Naturally, when we associate each gene cluster with a single function, we could end up with multiple gene clusters in a pangenome with the same function. From our experience, most functions are associated with a single gene cluster, but there are still plenty of functions that associate with multiple gene clusters, and this will be more common in pangenomes that contain distantly related genomes. In these cases, in order to find the occurrence of a given function among genomes, we simply 'merge' the occurrences of all gene clusters associated with that same function (for you computational readers, we simply take the "or" product of the presence/absence vectors of the gene clusters across genomes).

The careful readers will notice that we distinguish between 'functional annotation' and 'functional association' in the following text. When we mention 'functional annotation', we refer to the annotation of a single gene with a function by the functional annotation source (i.e. COGs, EggNOG, etc.), whereas 'functional association' of a gene cluster is the association of gene clusters with a single function as described above.

Ok, so now we have a presence/absence table of functions in genomes, and we can calculate different scores for each function, and also visualize it. See the details below!
</div>


Let's use the *Prochlorococcus* example to learn what we can do with this.

First we will compare the low-light vs. the high-light genomes in order to see if there are any functions that are unique to either group using the 'light' categorical additional layer data (if this doesn't make sense to you please go back to one of the pangenome figures above and see the layer additional data 'light'):

```bash
anvi-get-enriched-functions-per-pan-group -p PROCHLORO/Prochlorococcus_Pan-PAN.db \
                                          -g PROCHLORO-GENOMES.db \
                                          --category light \
                                          --annotation-source COG_FUNCTION \
                                          -o PROCHLORO-PAN-enriched-functions-light.txt
```

Here is the structure of the output file *PROCHLORO-PAN-enriched-functions-light.txt* (there are more columns, scroll towards right to see them):

|category|COG_FUNCTION|enrichment_score|weighted_enrichment_score|portion_occurrence_in_group|portion_occurrence_outside_of_group|occurrence_in_group|occurrence_outside_of_group|gene_clusters_ids|core_in_group|core|wilcoxon_p_value|wilcoxon_statistic|wilcoxon_corrected_p_value|
|:--:|:--|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
|LL|Ser/Thr protein kinase RdoA involved in Cpx stress response, MazF antagonist|1.00|29.09|1.00|0.00|11.00|0.00|GC_00002783, GC_00003936, GC_00004631, GC_00005468|True|False|0.00|4.54|0.00|
|LL|3-polyprenyl-4-hydroxybenzoate decarboxylase|1.00|29.09|1.00|0.00|11.00|0.00|GC_00001766, GC_00001810|True|False|0.00|4.54|0.00|
|LL|Exonuclease VII, large subunit|1.00|29.09|1.00|0.00|11.00|0.00|GC_00001726|True|False|0.00|4.54|0.00|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|
|LL | Periplasmic protein TonB, links inner and   outer membranes | 0.58 | 16.79 | 0.73 | 0.15 | 8 | 3 | GC_00002213, GC_00003430, GC_00005286,   GC_00006014, GC_00006187, GC_00006496, GC_00006928 | FALSE | FALSE | 0.01 | 2.62 | 0.17|
LL | Gamma-glutamyltranspeptidase | 0.55 | 15.87 | 0.55 | 0 | 6 | 0 | GC_00002294 | FALSE | FALSE | 0.01 | 2.48 | 0.23|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|
|LL | Predicted NAD/FAD-dependent oxidoreductase | -1 | -29.09 | 0 | 1 | 0 | 20 | GC_00000970 | FALSE | FALSE | 0 | -4.54 | 0|
|LL | ATP-dependent DNA ligase | -1 | -29.09 | 0 | 1 | 0 | 20 | GC_00000992, GC_00001635 | FALSE | FALSE | 0 | -4.54 | 0|
|HL | Redox-sensitive bicupin YhaK, pirin superfamily | 1 | 29.09 | 1 | 0 | 20 | 0 | GC_00001240 | TRUE | FALSE | 0 | 4.54 | 0|
|HL | N-acetylglutamate synthase or related acetyltransferase, GNAT family | 1 | 29.09 | 1 | 0 | 20 | 0 | GC_00001298 | TRUE | FALSE | 0 | 4.54 | 0|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|
|HL | 1,6-Anhydro-N-acetylmuramate kinase | -1 | -29.09 | 0 | 1 | 0 | 11 | GC_00001728 | FALSE | FALSE | 0 | -4.54 | 0|
|HL | Proteasome lid subunit RPN8/RPN11, contains Jab1/MPN domain metalloenzyme   (JAMM) motif | -1 | -29.09 | 0 | 1 | 0 | 11 | GC_00002219, GC_00003850, GC_00004483 | FALSE | FALSE | 0 | -4.54 | 0|

The following describes each column:

1. **category** is the column you chose from your [layers additional data table](http://merenlab.org/2017/12/11/additional-data-tables/#layers-additional-data-table), in order to identify groups within your genomes. In the Prochlorococcus case we have two light groups: low light (LL), and high light (HL). The output table is ordered according to this column first, and within each group, the table is ordered according to the enrichment score.

2. **COG_FUNCTION** this column has the name of the specific function for which enrichment was calculated. In this example we chose to use `COG_FUNCTION` for functional annotation, and hence the column title is `COG_FUNCTION`. You can specify whichever functional annotation source you have in your PAN database using the `--annotation-source`, and then the analysis would be done according to that annotation source. Even if you have multiple functional annotation sources in your genome storage, only one source could be used for a single run of this program. If you wish, you could run it multiple times and each time use a different annotation source. If you don't remember which annotation sources are available in your genomes storage, you can use `--list-annotation-sources`.

3. **enrichment_score** is a score to measure how much is this function unique to the genomes that belong to a specific group vs. all other genomes in your pangenome. It is simply the output of,<br /><br />
$$
portion\_occurrence\_in\_group -  portion\_occurrence\_outside\_of\_group
$$
<br /><br />In the comparison above, each genome belongs to one of the two groups (HL, LL), but if the column you chose from your layers additional data table has more than two groups, then when comparing a function for each group, the occurrence is compared between the group members, and the rest of the genomes, i.e. the comparison is not pair-wise between groups (you can see the example below of comparisonbetween clades of Prochlorococcus for more details). When the occurrence of a function in the group is lower than outside the group, then we get negative values (makes sense, right?).

4. **weighted_enrichment_score** is meant to help you compare between enrichment scores for different groups even if these groups have different sizes. The ides is that if you get a high enrichment score for a group of only two members it is not as impressive as getting high enrichment for a larger group of genomes. In addition, if the comparison is done between groups similar in size, then it is more reliable, i.e. if the number of genomes in a group is similar to the number of genomes outside that group, then the comparison is more valid. To account for these, the weighted enrichment is the output of the following equation:<br /><br />
$$
weighted\_enrichment\_score = enrihment \times N_{g} \times (p_{group} \times log_{2} (p_{group}) + (1-p_{group}) \times log_{2} (1 - p_{group}))
$$
<br /><br />Where $$N_{g}$$ is the total number of genomes in the pangenome, and $$p_{group}$$ is the portion of group members out of all genomes (i.e. the number of genomes that belong to the group divided by the total number of genomes in the pangenome). In other words the weighted enrichment score is the product of the enrichment score with the total number of genomes in the pangenome, multiplied by the entropy of group membership. If the number of group members is equal to the number of non-members, then the entropy is maximized (and equals 1), and if the groups are not balanced then the entropy is lower.

5. **portion_occurrence_in_group** is the number of genomes in the group that were associated with the function, divided by the total number of genomes in the group

6. **portion_occurrence_outside_of_group** is the number of genomes not in the group that were associated with the function, divided by the total number of genomes not in the group.

7. **occurrence_in_group** is the number of genomes in the group that were associated with the function.

8. **occurrence_outside_of_group** is the number of genomes outside the group that were associated with the function.

9. **gene_clusters_ids** are the gene clusters that were associated with this function. Notice that each gene cluster would be associated with a single function, but a function could be associated with multiple gene clusters.

10. **core_in_group** is "true" if the function occurs in all members of the group.

11. **core** is "true" if the function occurs in all genomes in the pangenome.

12. **wilcoxon_p_value** is the p value for a [Wilcoxon rank-sum test](https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test) to compare the occurrence of the function in the group vs. outside the group (we use [scipy.stats.ranksums](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ranksums.html)).

13. **wilcoxon_statistic** the statistic for the Wilcoxon rank sum test.

14. **wilcoxon_corrected_p_value** a correction for multiple tests using the [Benjamini-Hochberg false discovery rate method](https://en.wikipedia.org/wiki/False_discovery_rate#Benjamini%E2%80%93Hochberg_procedure). This is done using [statsmodels.stats.multitest.multipletests](http://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html#statsmodels.sandbox.stats.multicomp.multipletests).

Now let's search for the top function in the table "Ser/Thr protein kinase RdoA involved in Cpx stress response, MazF antagonist", which is enriched for the members of the LL group, and we can see in the table that it matches four gene clusters.

[![layers]({{images}}/ser_thr_kinase_gene_clusters.png)]({{images}}/ser_thr_kinase_gene_clusters.png){:.center-img .width-60}

In fact, if we look carefuly, then we find that this function matches a gene cluster that is unique for each one of the four low light clades, and is single copy core gene for the low light genomes in our pangenome. Cool!

Let's look at another function from the table: `Exonuclease VII, large subunit` (third line). When we search this function, we should be careful since the name of the function contains a comma, and the function search option in the interactive interface treats the comma as if it is separating multiple functions that are to be searched at the same time. Hence I just searched for `Exonuclease VII`, and here are the results:


[![layers]({{images}}/Exonuclease-VII.png)]({{images}}/Exonuclease-VII.png){:.center-img .width-40}

We can see that the search matched hits for both Exonuclease VII, large and small subunits, with a total of 22 hits.

[![layers]({{images}}/Exonuclease-VII-2.png)]({{images}}/Exonuclease-VII-2.png){:.center-img .width-60}

The large subunit matches a single gene cluster which is in the CORE LL, and the small subunit matches a gene cluster in each one of the clade specific cores (similar to what we saw above for the Ser/Thr protein kinase). Both of these genes are also part of the single copy core unique to low light members and absent from high light members.

### Creating a quick pangenome with functions

Next, we will explore whether there are any functions enriched for any of the sub clades. In addition, we will introduce another feature `--functional-occurrence-table-output`. This optional output is a TAB-delimited file with the presence/absence information for functions in genomes.

```bash
anvi-get-enriched-functions-per-pan-group -p PROCHLORO/Prochlorococcus_Pan-PAN.db \
                                          -g PROCHLORO-GENOMES.db \
                                          --category clade\
                                          --annotation-source COG_FUNCTION \
                                          -o PROCHLORO-PAN-enriched-functions-clade.txt \
                                          --functional-occurrence-table-output PROCHLORO-functions-occurrence.txt
```

Let's look at some results. ***This is how PROCHLORO-PAN-enriched-functions-clade.txt*** looks like:

| category | COG_FUNCTION                                                  | enrichment_score | weighted_enrichment_score | portion_occurrence_in_group | portion_occurrence_outside_of_group | occurrence_in_group | occurrence_outside_of_group | gene_clusters_ids        | core_in_group | core  | wilcoxon_p_value | wilcoxon_statistic | wilcoxon_corrected_p_value |
|----------|---------------------------------------------------------------|------------|---------------------|----------------------------|------------------------------------|--------------------|----------------------------|--------------------------|---------------|-------|------------------|--------------------|----------------------------|
| LL_II    | Archaeal DNA polymerase II, large subunit                     | 1          | 19.76               | 1                          | 0                                  | 5                  | 0                          | GC_00002540              | TRUE          | FALSE | 0                | 3.49               | 0.06                       |
| LL_II    | Outer membrane protein assembly factor BamD, BamD/ComL family | 1          | 19.76               | 1                          | 0                                  | 5                  | 0                          | GC_00002403              | TRUE          | FALSE | 0                | 3.49               | 0.06                       |
| LL_II    | tRNA G37 N-methylase Trm5                                     | 0.96       | 19                  | 1                          | 0.04                               | 5                  | 1                          | GC_00002586, GC_00007078 | TRUE          | FALSE | 0                | 3.36               | 0.06                       |
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|

And this is how ***PROCHLORO-functions-occurrence.txt*** looks like:

|| AS9601                                                                                                      | CCMP1375 | EQPAC1 | GP2 | LG | MED4 | MIT9107 | MIT9116 | MIT9123 | MIT9201 | MIT9202 | MIT9211 | MIT9215 | MIT9301 | MIT9302 | MIT9303 | MIT9311 | MIT9312 | MIT9313 | MIT9314 | MIT9321 | MIT9322 | MIT9401 | MIT9515 | NATL1A | NATL2A | PAC1 | SB | SS2 | SS35 | SS51 |
|-------------------------------------------------------------------------------------------------------------|----------|--------|-----|----|------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|--------|--------|------|----|-----|------|------|---|
| 3-deoxy-D-manno-octulosonate 8-phosphate phosphatase KdsC and related HAD superfamily phosphatases          | 0        | 0      | 0   | 0  | 0    | 0       | 0       | 0       | 0       | 0       | 0       | 0       | 0       | 0       | 0       | 1       | 0       | 0       | 1       | 0       | 0       | 0       | 0       | 0      | 0      | 0    | 0  | 0   | 0    | 0    | 0 |
| Creatinine amidohydrolase/Fe(II)-dependent formamide hydrolase involved in riboflavin and F420 biosynthesis | 1        | 1      | 1   | 1  | 1    | 1       | 1       | 1       | 1       | 1       | 1       | 1       | 1       | 1       | 1       | 1       | 1       | 1       | 1       | 1       | 1       | 1       | 1       | 1      | 1      | 1    | 1  | 1   | 1    | 1    | 1 |
| RNA recognition motif (RRM) domain                                                                          | 1        | 1      | 1   | 1  | 1    | 1       | 1       | 1       | 1       | 1       | 1       | 1       | 1       | 1       | 1       | 1       | 1       | 1       | 1       | 1       | 1       | 1       | 1       | 1      | 1      | 1    | 1  | 1   | 1    | 1    | 1 |
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|

Let's use the functional occurrence table for visualization. First we fix the names of functions to get rid of things like commas etc. Basically, we only keep alphanumeric characters, and replace any sequence of non-alphanumeric characters by a single `_`.

```bash
sed "s/[^[:alnum:]	_]/_/g" PROCHLORO-functions-occurrence.txt |\
	tr -s \_ _ | sed 's/^	/name	/' \
	> PROCHLORO-functions-occurence-fixed-a-little.txt
```

Unfortunately there could be some functions with very similar names, let's check if there are any duplicate names now:

```
cut -f 1 PROCHLORO-functions-occurence-fixed-a-little.txt | sort | uniq -d
```

We can see that `Fatty_acid_desaturase`, and `Protein_tyrosine_phosphatase` are duplicated now. This is because, for example, there are two COG functions `Protein-tyrosine phosphatase` (accession `COG2453`) and `Protein-tyrosine-phosphatase` (accession `COG0394`), which match different gene clusters in our pangenome. Because we cannot create a tree with duplicated nodes, and since we can't really say that these are different functions, we will use a script to merge these duplicated occurences (i.e. merge their occurences with a `logical or`). Your working directory includes the script `fix_functional_occurence_table.py`:

```
./fix_functional_occurence_table.py PROCHLORO-functions-occurence-fixed-a-little.txt PROCHLORO-functions-occurence-fixed.txt
```

Then we create trees for the interactive interface:

```bash
anvi-matrix-to-newick PROCHLORO-functions-occurrence-fixed.txt \
                      -o PROCHLORO-functions-tree.txt

anvi-matrix-to-newick PROCHLORO-functions-occurrence-fixed.txt \
                      -o PROCHLORO-functions-layers-tree.txt \
                      --transpose
```

We run a quick dry run to create a manual mode profile database.

```bash
anvi-interactive -p PROCHLORO-functions-manual-profile.db \
                 --tree PROCHLORO-functions-tree.txt \
                 --data PROCHLORO-functions-occurrence-fixed.txt \
                 --manual \
                 --dry-run
```

We import the tree for the layers order:

```bash
echo -e "item_name\tdata_type\tdata_value" > PROCHLORO-functions-layers-order.txt
echo -e "PROCHLORO_functions_tree\tnewick\t`cat PROCHLORO-functions-layers-tree.txt`" \
                             >> PROCHLORO-functions-layers-order.txt

anvi-import-misc-data PROCHLORO-functions-layers-order.txt \
                      -p PROCHLORO-functions-manual-profile.db \
                      -t layer_orders \
                      --just-do-it
```

We can get some information about the genomes from the PAN database:

```bash
anvi-export-misc-data -p PROCHLORO/Prochlorococcus_Pan-PAN.db \
                      -t layers \
                      -o PROCHLORO-layer-additional-data.txt
```

And then import that into our manual mode profile database:

```bash
anvi-import-misc-data PROCHLORO-layer-additional-data.txt \
                      -p PROCHLORO-functions-manual-profile.db \
                      -t layers
```

The work directory you downloaded includes a nice state that we created for this profile database along with a collection, let's import these:

```bash
anvi-import-state -p PROCHLORO-functions-manual-profile.db \
                  -s PROCHLORO-manual-default-state.json \
                  -n default
    
anvi-import-collection -p PROCHLORO-functions-manual-profile.db \
                       -C default \
                       PROCHLORO-functions-collection.txt
```

And now we can take a look:

```bash
anvi-interactive -p PROCHLORO-functions-manual-profile.db \
                 -t PROCHLORO-functions-tree.txt \
                 -d PROCHLORO-functions-occurrence-fixed.txt \
                 --title "Prochlorococcus Pan - functional occurrence" \
                 --manual
```

[![layers]({{images}}/Functional-occurrence.png)]({{images}}/Functional-occurrence.png){:.center-img .width-60}

A collection of 869 core functions emmerges.

We can also see that the occurrence of functions is recapitulating all four LL clades. In contrast, the two HL clades seem to be mixed together.

A little BASH one-liner can give us all gene clusters associated with core functions:

``` bash
grep LL PROCHLORO-PAN-enriched-functions-light.txt |\
  awk  -F $'\t' '$11 == "True" { print $9 }' |\
    tr ',' '\n' |\
      sed 's/ //g' > core_functions_gene_clusters.txt
```

And we get a file with 2613 gene clusters. Now we can compare this to the collection of gene clusters we have [above](#displaying-the-pan-genome). We can find, for example, how many of the gene clusters that are in the CORE_LL are in the functional core of all 31 Prochlorococcus genomes:

``` bash
for gc in `grep CORE_LL PROCHLORO-PAN-default-collection.txt`
do
    grep $gc core_functions_gene_clusters.txt
done > CORE_LL_included_in_functional_core.txt
```

We find that 103 of the 144 gene clusters are part of the functional core. And for HL it is 294 of 499 that are found to be part of the functional core. An important thing to remember is that gene clusters that have no function associated with them are not included in this analysis.

We can find all the gene clusters that are associated with a function:

``` bash
grep LL PROCHLORO-PAN-enriched-functions-light.txt |
  awk  -F $'\t' '{ print $9 }' |
    tr ',' '\n' |
      sed 's/ //g' > all_gene_clusters_with_functions.txt
```

There are 3629 gene clusters with functions. How many gene clusters that belong to the CORE_HL have functions?

```
for gc in `grep CORE_HL PROCHLORO-PAN-default-collection.txt`
do
   grep $gc all_gene_clusters_with_functions.txt;
done > HL_CORE_GC_with_functions.txt

wc -l all_gene_clusters_with_functions.txt
321
```

There are 321 (of a total of 499) gene clusters in CORE_HL that have functions. Hence many of the gene clusters in CORE_HL that were not found to be part of the functional core, simply don't have any functional annotation.

We hope you find great uses for the functional enrichment analysis framework in anvi'o pangenomic workflow! Feel free to send us any questions you may have.


## Computing the average nucleotide identity for genomes

Anvi'o also contains a program, [`anvi-compute-ani`](http://merenlab.org/software/anvio/vignette/#anvi-compute-ani), which uses [PyANI](https://doi.org/10.1039/C5AY02550H) to compute average nucleotide identity across your genomes. It expects you to provide an external genomes file to declare which genomes you are interested in for an ANI analysis, but it also optionally accepts a pan database to add its findings into it as additional layer data.

Here is an example with our *Prochlorococcus* Pan genome:

``` bash
anvi-compute-ani --external-genomes external-genomes.txt \
                 --output-dir ANI \
                 --num-threads 6 \
                 --pan-db PROCHLORO/Prochlorococcus_Pan-PAN.db
```

Once it is complete, we can visualize the pan genome again to see what is new:

``` bash
anvi-display-pan -g PROCHLORO-GENOMES.db \
                 -p PROCHLORO/Prochlorococcus_Pan-PAN.db
```

When you first look at it you will not see anything out of ordinary. But if you go to the 'Layers' tab, you will see the following addition:


[![layers]({{images}}/layer-groups-ani.png)]({{images}}/layer-groups-ani.png){:.center-img .width-60}

If you click the checkbox for, say, `ANI_percentage_identity`, you will see a new set of entries will be added into the list of layer data entries. Then, if you click the button 'Draw' or 'Redraw layer data', you should see the ANI added to your display:

[![prochlorococcus-ani]({{images}}/prochlorococcus-ani.png)]({{images}}/prochlorococcus-ani.png){:.center-img .width-60}

Yes. [Magic](https://github.com/merenlab/anvio/pull/822).

{:.notice}
You may need to change the `min` value from the interface for a better representation of ANI across your genomes in your own pangenome.

{:.notice}
***A note from Meren***: We have `anvi-compute-ani`, because someone asked for it on GitHub. We thank our proactive users, like [Mike Lee](https://astrobiomike.github.io/), on behalf of the community.


<div class="extra-info" markdown="1">

<span class="extra-info-header">A little note from Meren on how to order genomes</span>

It has always been a question mark for me how to order genomes with respect to gene clusers. Should one use the gene cluster presence/absence patterns to organize genomes (where one ignores paralogs and works with a binary table to cluster genomes)? Or should one rely on gene cluster frequency data to order things (where one does consider paralogs, so their table is not binary anymore, but contains the frequencies for number of genes from each genome in each gene cluster)? 

Thanks to Özcan’s [new addition](https://github.com/merenlab/anvio/commit/aa007cf902dea2de4bd63524cd49f0566cf2511d) to the codebase while I was working on this section of the tutorial --so the ANI matrix is always ordered based on the order of genomes in the pan display, I made an observation that shed some light on the question of how to order genomes.

This is how ANI matrix looks like when I order *Prochlorococcus* genomes based on gene cluster frequency data: 

[![ani]({{images}}/ani-gc-freqs.png)]({{images}}/ani-gc-freqs.png){:.center-img .width-40}

In contrast, this is how ANI matrix looks like when I order genomes based on gene cluster presence/absence data:

[![ani]({{images}}/ani-gc-preabs.png)]({{images}}/ani-gc-preabs.png){:.center-img .width-40}

There you have it.

If you are working with isolate genomes, maybe you should consider ordering them based on gene cluster frequencies since you *would* expect a proper organization of genomes based on gene clusters should often match with their overall similarity. Here are some small notes:

* In this case it is easy to see that gene cluster presence/absence data is doing a poor job since we know the *Prochlorococcus* clades to see that clearly (i.e., gene cluster presence/absence data is mixing up high-light group II (purple) and high-light group I (poopcolor) genomes. But when you don't know to which clades your genomes belong, gene cluster frequencies can be more accurate to organize them nicely. 

* It is sad to know that MAGs will often lack multi-copy genes, and they will have less paralogs (since assembly of short reads create very complex de Bruijn graphs which obliterate all redundancy in a given genome), hence their proper organization will not truly benefit from gene cluster frequency data, at least not as much as isolate genomes, and the errors we see with gene cluster presence/absence data-based organization will have a larger impact on our pangenomes.

* See the emergence of those sub-clades in high-light group II? If you are interested, please take a look at [our paper on *Prochlorococcus* metapangenome](https://peerj.com/articles/4320/) (especially the section "**The metapangenome reveals closely related isolates with different levels of fitness**"), in which we divided HL-II into multiple subclades and showed that despite small differences between these genomes, those sub-clades displayed remarkably different environmental distribution patterns. Sub-clades we could define in that study based on gene cluster frequencies and environmental distribution patterns match quite well to the groups ANI reveals. This goes back to some of the most fundamental questions in microbiology regarding how to define ecologically relevant units of microbial life. I don't know if we are any close to making any definition, but I can tell you that those 'phylogenetic markers' we are using ... I am not sure if they are really working that well (smiley).

</div>


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

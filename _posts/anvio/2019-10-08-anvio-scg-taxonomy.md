---
layout: post
authors: [meren]
title: "Real-time genome taxonomy estimation with anvi'o"
excerpt: "An attempt in alchemy combining the magic of GTDB and SCGs in anvi'o."
modified: 2019-10-08
tags: []
categories: [anvio]
comments: true
redirect_from:
  - /scg-taxonomy/
---

{% include _toc.html %}

{% include _project-anvio-version.html %}

{:.warning}
The shortened URL for this page is [{{ site.url }}/scg-taxonomy]({{ site.url }}/scg-taxonomy)

{: .notice}
This workflow is for `v6` and later versions of anvi'o. You can identify which version you have on your computer by typing `anvi-self-test --version` in your terminal.

{% capture images %}{{site.url}}/images/anvio/2019-10-08-anvio-scg-taxonomy{% endcapture %}

With the [anvi'o](https://peerj.com/articles/1319/) scg-taxonomy workflow you can learn,

* Taxonomic affiliation of a single genome
* Taxonomic profile of an entire assembled metagenome
* Relative abundances of taxa in a given assembled metagenome
* Taxonomic profile of an anvi'o collection (i.e., binning results)
* Relative abundances of taxa found in a colleciton across metagenomes

{:.notice}
You can use anvi'o scg-taxonomy workflow even if you haven't done any metagenomic work with anvi'o. All you need is an anvi'o installation, and a FASTA file for each of your genomes.

## Thanks

The new anvi'o scg-taxonomy module described on this page stands on more than 1,500 lines of code that is fully integrated to the rest of the anvi'o platform. This effort is partially a result of a pilot award from the [FACCTS](https://fcc.uchicago.edu/) (*France and Chicago Collaborating in The Sciences*) program, which supported the travel expenses of **Quentin Clayssen** from France to Chicago, who put the initial effort to investigate the feasibility of this idea. I would also like to thank [**Alon Shaiber**](https://twitter.com/alon_shaiber) for his excellent suggestions that influenced the speed of this algorithm, and [**Özcan Esen**](https://twitter.com/ozcanesen) for helping with the integration of this module to the anvi'o interface.

## Introduction

Often when we deal with genomes we need some sort of insight into taxonomy. This becomes a critical need especially when you reconstruct genomes from metagenomes as you want to know where do these new genomes fit in the tree of life rapidly. So far anvi'o did not offer a rapid solution for that, but starting with the `v6` we have a solution for that. Everything is new, so probably there are going to be hickups before we manage to stabilize things.

The new anvi'o scg-taxonomy worklow uses the taxonomy determined by the [The Genome Taxonomy Database](https://gtdb.ecogenomic.org/) (GTDB, [@ace_gtdb](https://twitter.com/ace_gtdb)), an effort led by scientists who are primarily affiliated with the [Australian Centre for Ecogenomics](https://ecogenomic.org/).

{:.warning}
If you use anvi'o scg-taxonomy module, please don't forget to cite the initial work that introduced this resource: [doi:10.1038/nbt.4229](http://doi.org/10.1038/nbt.4229). There is also a tool from the developers of the GTDB, [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk), to assign taxonomy to genomes.

In a nutshell, anvi'o (1) downloads both single-copy core genes (SCGs) and the taxonomy of the genomes -as defined by the GTDB- from which these genes are coming from, (2) bulids search databases for a subset of these SCGs (depending on whether they are in the Bacterial and Archaeal collection of SCGs anvi'o recognizes and based on other heuristics), (3) searches SCGs from every new genome in the anvi'o ecosystem in these databases and creates consensus taxonomy for each genome.

There following sections will cover the three major steps for anvi'o scg-taxonomy: (1) setting things up, (2) population contigs database with taxonomy, and (3) running taxonomy estimations. The text will cover how to run anvi'o scg-taxonomy both from the terminal environment and through the graphical user interface.

For the sake of reproducibility, here I will use the infant gut dataset that [we use to demonstrate many things anvi'o](/tutorials/infant-gut/). If you would like to follow this tutorial not with your own datasets but with the infant gut dataset, you can simply [click here first](/tutorials/infant-gut/#downloading-the-pre-packaged-infant-gut-dataset), follow the download instructions, and come back here once the text tells you to **go back**.


## Setting up anvi'o SCG taxonomy

{:.notice}
This is something you will do only once after installing anvi'o `v6`.

A fresh installation of anvi'o will have all the FASTA files that are going to be used for scg-taxonomy, but not the search databases that must be generated from them. We are not shipping the databases because they need to be built with the DIAMOND installed on your computer to make sure we are not running into various incompatibilty problems.

To setup your SCG taxonomy databases, run the following program, which will not take longer than a minute:

```
anvi-setup-scg-databases
```

{:.notice}
Please see the help menu for every program even if the text does not ask you to do it. There are many parameters for every anvi'o program so you can tweak things depending on your system and needs.

When I run this on my computer, this is the output I get:

```
(anvio-master) (base) meren ~ $ anvi-setup-scg-databases

WARNING
===============================================
Please remember that the data anvi'o uses for SCG taxonomy is a courtesy of The
Genome Taxonomy Database (GTDB), an initiative to establish a standardised
microbial taxonomy based on genome phylogeny, primarly funded by tax payers in
Australia. Please don't forget to cite the original work, doi:10.1038/nbt.4229
by Parks et al to explicitly mention the source of databases anvi'o relies upon
to estimate genome level taxonomy. If you are not sure how it should look like
in your methods sections, anvi'o developers will be happy to help you if you
can't find any published example to get inspiration.


WARNING
===============================================
Anvi'o found your FASTA files in place, but not the databases. Now it will
generate all the search databases using the existing FASTA files.

* Every FASTA is now turned into a fancy search database. It means you are now
allowed to run `anvi-run-scg-taxonomy` on anvi'o contigs databases.
```

If you don't have permission to write to the part of the disk anvi'o is installed on your system, you should ask your system administrator to run this command for you.

The command above will only create the databases without downloading original files. In some rare cases you may have to re-download everything and start with a clean slate. If you think this is what you need, try this command (but you should run this only if you have a problem):

```
anvi-setup-scg-databases --reset
```

## Populating contigs db with SCG taxonomy

{:.notice}
This is something you will do once for every contigs database you wish to work with.

For anvi'o scg-taxonomy to work, you need to run the following program on your contigs database:

```
anvi-run-scg-taxonomy -c CONTIGS.db
```

{:.warning}
If you have multiple CPUs and memory you may want to use `--num-parallel-processes` and `--num-threads` parameters to speed things up. Please refer to the help menu to find out what they actually mean. Just ot give you an idea, runing this program on the infant gut dataset without any parallelization takes 50 seconds. When I run the same command with `--num-parallel-processes 3` and `--num-threads 3`, it takes 10 seconds. Just so you know.

Here is a crude list of what is going on behind the scenes when you run this program:

- Anvi'o finds all [the 22 single-copy core genes](https://github.com/merenlab/anvio/tree/master/anvio/data/misc/SCG_TAXONOMY/GTDB/SCG_SEARCH_DATABASES) it uses for scg-taxonomy in a given contigs database, 
- Searches them against the corresponding databases from the GTDB, and,
- Chooses the best hit among all the hits, and puts stores them in a table in the contigs database.

{:.notice}
If you run the program with `--debug` flag you can see all the underlying taxonomy strings that lead to the consensus. Beware, though, your terminal will be overwhelmed by some text :)

## Estimating taxonomy in the terminal

There are many ways to estimate taxonomy, but all of them are going to use the program `anvi-estimate-scg-taxonomy`. Following sections will demonstrate multiple uses of this tool.

### Contigs db of a single genome

If you have an anvi'o contigs for a single genome, you can simply run this command:

```
anvi-estimate-scg-taxonomy -c CONTIGS.db
```

If your contigs database is a metagenome, rather than a single genome, this command line will yield an error that complains about the fact that the redundancy of the single-copy core genes is too much for this to be a single genome. In this case you can either force anvi'o to give you an answer ANYWAY by using the flag `--just-do-it` (which will definitely not going to be a useful answer) or tell anvi'o to treat this contigs database as a metagenome.

### Contigs db of a metagenome

If you have generated your contigs database from a metagenomic assembly, then you should be running the same command in metageome mode:

```
anvi-estimate-scg-taxonomy -c CONTIGS.db \
                           --metagenome-mode
```

Here is a simplified list of what is going on behind the scenes when you run this command:

- Anvi'o counts the occurrence of each of the 22 single-copy core genes in this contigs database,
- Chooses one of them to use for taxonomy estimation, and,
- Generates an output (which can also be stored into a tab delimited file with `--output` parameter).

Since the contigs database of the infant gut dataset is generated from the assembly of all short reads, this is in fact quite an appropriate way to study it with `anvi-estimate-scg-taxonomy`. When you run the command, anvi'o will first tell you the frequencies of SCGs in your contigs database. Here is an example output from my run on the infant gut dataset:

```
* A total of 171 single-copy core genes with taxonomic affiliations were
successfuly initialized from the contigs database 🎉 Following shows the
frequency of these SCGs: Ribosomal_S6 (10), Ribosomal_S8 (10), Ribosomal_L27A
(10), Ribosomal_S3_C (9), Ribosomal_S9 (9), Ribosomal_L6 (9), Ribosomal_L9_C
(9), Ribosomal_L16 (9), ribosomal_L24 (9), Ribosomal_L1 (8), Ribosomal_L2 (8),
Ribosomal_L4 (8), Ribosomal_S20p (7), Ribosomal_L13 (7), Ribosomal_L17 (7),
Ribosomal_L22 (7), Ribosomal_S2 (6), Ribosomal_S11 (6), Ribosomal_L3 (6),
Ribosomal_L20 (6), Ribosomal_L21p (6), Ribosomal_S7 (5).
```

Because there are 10 `Ribosomal_S6` genes, anvi'o will automatically use that one to offer a taxonomy table. For this dataset, that output looks like this:

[![scg-taxonomy]({{images}}/igd-metagenome-mode.png)]({{images}}/igd-metagenome-mode.png){:.center-img .width-90}

Not bad, eh?

{:.warning}
One of the most significant shortcomings of this workflow is that you will not find taxonomic signal for those populations that did not assemble. 

I could also set the SCG myself. For instance I can see in the output that the SCG `Ribosomal_S9` occurs 9 times instead of 10. If I try that SCG instad,

```
anvi-estimate-scg-taxonomy -c CONTIGS.db \
                           --metagenome-mode \
                           --scg-name Ribosomal_S9
```

I get a similar picture of the same data, but of course I am missing one of the entries:

[![scg-taxonomy]({{images}}/igd-metagenome-mode-S9.png)]({{images}}/igd-metagenome-mode-S9.png){:.center-img .width-90}

These are not conclusive taxonomic insights, but still very useful to have a rapid idea about what you have.

### Contigs db + profile db

If you have a contigs database for your genomes or metagenomes AS WELL AS a single or merged profile database that you have generated from the read recruitment analyses of your metageomes given those contigs (i.e., the standard metagenomic workflow), in fact you can also estimate the coverages of those SCGs you are assigning taxonomy to. To do that, I can run the following command on the infant gut dataset:

```
anvi-estimate-scg-taxonomy -c CONTIGS.db \
                           -p PROFILE.db \
                           --metagenome-mode \
                           --compute-scg-coverages
``` 

Tadaa:

[![scg-taxonomy]({{images}}/igd-metagenome-mode-coverages.png)]({{images}}/igd-metagenome-mode-coverages.png){:.center-img .width-90}

Almost instantenously, you learn about the coverages of a given SCG across samples.

{:.notice}
You can always use `--output-file` parameter to store the full output into a TAB-delmited file.

### Contigs db + profile db + collection

When we have a metagenomic assembly, in most cases we either perform manual binning and store our bins in a collection, or use an external binning algorithm to import a collection using the program `anvi-import-collection` into our profile database. In the infant gut dataset there is already a collection, so you can import that for testing:

```
anvi-import-collection additional-files/collections/merens.txt \
                       -p PROFILE.db \
                       -c CONTIGS.db \
                       -C merens
```

And you can always check collections in your profile databases easily. For instance this is what I get after importing the collection `merens`:

```
anvi-show-collections-and-bins -p PROFILE.db

Collection: "CONCOCT"
===============================================
Collection ID ................................: CONCOCT
Number of bins ...............................: 34
Number of splits described ...................: 4,784
Bin names ....................................: Bin_1, Bin_10, Bin_11, Bin_12, Bin_13, Bin_14, Bin_15, Bin_16, Bin_17, Bin_18, Bin_19, Bin_2, Bin_20, Bin_21, Bin_22, Bin_23, Bin_24, Bin_25, Bin_26, Bin_27, Bin_28, Bin_29, Bin_3, Bin_30, Bin_31, Bin_32, Bin_33, Bin_34, Bin_4, Bin_5, Bin_6, Bin_7, Bin_8, Bin_9

Collection: "merens"
===============================================
Collection ID ................................: merens
Number of bins ...............................: 13
Number of splits described ...................: 4,451
Bin names ....................................: Aneorococcus_sp, C_albicans, E_facealis, F_magna, L_citreum, P_acnes, P_avidum, P_rhinitidis, S_aureus, S_epidermidis, S_hominis, S_lugdunensis, Streptococcus
```

The bin names in the collection meren are taxonomic affiliations for these bins as they appeared in the original publicatoin of this dataset by [Sharon et al.](https://www.ncbi.nlm.nih.gov/pubmed/22936250), so they can serve as ground truths for this little experiment. Now we can run the scg-taxonomy estimation on this collection:

```
anvi-estimate-scg-taxonomy -c CONTIGS.db \
                           -p PROFILE.db \
                           -C merens
```

which returns this:

[![scg-taxonomy]({{images}}/igd-collection.png)]({{images}}/igd-collection.png){:.center-img .width-90}

Everything seems to be working well. But it is important to remember two things. One, this is a simple human gut metagenome and for your environmental metageomes things will not be as smooth (but those of us who study environmental metagenomes are already aware of that). Two, please note the last two lines. Why are there blank? If you were to look at the completion estimates of the bins in this collection,

```
anvi-estimate-genome-completeness -p PROFILE.db -c CONTIGS.db -C merens
```

You would see that `P_acnes` has no completion estimate:

[![scg-taxonomy]({{images}}/igd-collection-completion.png)]({{images}}/igd-collection-completion.png){:.center-img .width-90}

It doesn't have a completion estimate, becasue it doesn't have enough number of SCGs according to anvi'o to even make a suggestion regarding how complete it is. In contrast, `C_albicans` has a fairly well completion. But even though BUSCO HMMs [Tom Delmont had curated]({% post_url anvio/2018-05-05-eukaryotic-single-copy-core-genes %}) for anvi'o predicts its completion well, this bin is a eukaryotic bin. Sadly scg-taxonomy will not work for eukaryotic genomes (until GTDB starts considering them).


## Estimating taxonomy in the interface

In other words, the fun stuff.

We embarked upon this journey because we needed something extremely fast and relatively accurate while working with metageomes. Identifying a genome bin in anvi'o is extremely easy, but figuring out the taxonomy of these bins in the interface has not been possible. The scg-taxonomy subsystem bridges that gap now. For instance, if I initiate the interactive interface for the infant gut metagenome,

```
anvi-interactive -p PROFILE.db \
                 -c CONTIGS.db
```

and click the draw button, switch the Bin panel, I will see this display:

[![scg-taxonomy]({{images}}/igd-vanilla.png)]({{images}}/igd-vanilla.png){:.center-img .width-90}

Look what happens if I enable the real-time taxonomy estimation for my bins:

[![scg-taxonomy]({{images}}/igd-selection.gif)]({{images}}/igd-selection.gif){:.center-img .width-90}

Well. I don't know how to put this, but **this is real-time**. This is how fast a new set of contigs are affiliated with taxonomy in the interface.

Here is what happens when I load the colleciton `merens` to this interface:

[![scg-taxonomy]({{images}}/igd-collection.gif)]({{images}}/igd-collection.gif){:.center-img .width-90}

It is happening so rapidly, it looks as if anvi'o knew the taxonomy for these bins. But in reality the taxonomy shown underneath bin names are comoputed on the fly.

One thing to remember is that the taxonomy information shown on the bins menu is a consensus taxonomy calculated from the taxonomic affiliations of all relevant SCGs in a given bin. It is possible to see the underlying data for this through the taxonomy table interface you can open by clicking "Recalculate / Show Taxonomy for Bins" button in the Bins panel. When you click that, you see the following table that shows the consensus taxonomy, and the SCGs from which it was computed:

[![scg-taxonomy]({{images}}/igd-taxonomy-table.gif)]({{images}}/igd-taxonomy-table.gif){:.center-img .width-90}

## Final words and some warnings

You will not get perfect taxonomy with scg-taxonomy, but you will some taxonomy *very* rapidly. Anvi'o scg-taxonomy module will not estimate taxonomy for bins that are very low completion. It will not estimate taxonomy for eukaryotic bins. In `--metagenome-mode` you will not see any evidence of the populations from which the assembler did not recover any of the ribosomal proteins anvi'o uses for taxonomy estimation.

Yes. You will not get perfect taxonomy, but the reality is nothing will give you perfect taxonomy. There is no such thing as perfect taxonomy. Taxonomy is a mess. Thanks to heroic efforts by those who try to find an order in this mess (such as those who gave us SILVA or GTDB), we, as a community, try to find the best place to attach our genomes in the tree of life. But from building trees to determine single-copy core genes to identify them in genomes to define and assign taxonomic names, everything slowly falls apart as we struggle between sensitivity, feasibility, and accuracy of our approaches. For those of us who need taxonomy, the struggle is real. Our best option is to make sure we confirm the taxonomic affiliations of our genomes and bins. Which, despite the popular belief, may require manual work in some occasions using pangenomics, phylogenomics, studing ANI, BLAST'ing things around, etc. But I still hope the anvi'o scg-taxonomy will help you to develop some initial insights into *your* taxonomic mess 😇

I would also suggest you to take look at [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk). It will almost certainly do a better job to estimate taxonomy (we are hoping to publish a small note with proper benchmarks on accuracy). The GitHub page of GTDB-Tk suggests that it requires ~100Gb of memory, ~27Gb of storage, and takes about ~1 hour per 1,000 genomes when using 64 CPUs. Anvi'o scg-taxonomy takes about 100Mb space on a laptop computer wiht 16Gb memory and 4 CPUs, and takes about 20 seconds to assign taxonomy to 1,000 genomes. Of course this time doesn't include standard items in anvi'o workflow that takes place before anvi'o scg-taxonomy but is required for it to run, such as the time it takes to create an anvi'o contigs database, running HMMs on it, etc. But regardless, you should always be extra careful with things that are fast.

Both for `anvi-run-scg-taxonomy` and `anvi-estimate-scg-taxonomy` programs accept `--debug` flags. If things don't make sense, you should re-run these programs with the `--debug` flag to see additional information.

This is an extremely new feature and we hope it will contribute to your journey in microbial 'omics and your anvi'o experience. Please keep in mind just like everyting in anvi'o, the evolution of this feature will be heavily influenced by your feedback and suggestions. Please let us know if something doesn't work for you!

{% include _join-anvio-slack.html %}

{:.notice}
{% include _fixthispage.html source="_posts/anvio/2019-10-07-anvio-scg-taxonomy.md" %}







































































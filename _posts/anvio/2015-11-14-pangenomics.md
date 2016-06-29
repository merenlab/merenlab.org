---
layout: post
title: "An anvi'o workflow for microbial pangenomics"
excerpt: "The user-friendly interface anvi'o provides to work with pangenomes."
modified: 2015-10-14
tags: []
categories: [anvio]
comments: true
authors: [meren, tom]
---

{% include _toc.html %}
Cultivation of closely related microorganisms, and the subsequent recovery of their genomes, revealed that not all isolates share the same functional traits for a given population. On the other hand, individual isolates do not fully echo the functional complexity of naturally occurring microbial populations either. Large genomic collections is the only way to get closer to a more realistic picture of the pool of functions in a given clade of bacterial tree.

Overlapping and differing functions among the genomes of closely related organisms led to the introduction of a new concept, "pangenomics" [(Tettelin et al., 2005)](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1216834/), where genes across genomes are segregated into distinct groups: (1) the core-genome, (2) genes detected in multiple yet not in all genomes, and finally (3) isolate-specific genes that are detected only in a single genome. Pangenomic investigations are now widely used to dissect the functional traits of microorganisms, and to uncover environmentally- and clinically-important clusters of genes or functions.Although [multiple bioinformatics software](http://omictools.com/pangenomics-c1590-p1.html) are available to generate and/or visualize pangenomes, these solutions do not necessary offer flexible work environments, and hence limit the user's ability to interact with their data. 

In this post we will demonstrate the [anvi'o](https://peerj.com/articles/1319/) workflow for pangenomics using a simple dataset. Although we are still [working](http://github.com/meren/anvio) on the pan-genomic workflow, using the [v2 version of anvi'o]({% post_url anvio/2016-06-26-installation-v2 %}) you can,

* **Identify protein clusters**,
* **Visualize** the distribution of protein clusters across your genomes,
* **Organize** genomes based on shared protein clusters,
* Interactively identify **core**, and **core-like** protein clusters,
* Extend analysis results with **contextual information** about your genomes,
* **Compare MAGs** from your projects **with cultivar genomes** from other sources.

{:.notice}
You can use anvi'o for pangenomic analysis of your genomes even if you haven't done any metagenomic work with anvi'o. All you need is an anvi'o installation and a FASTA file for each of your genomes.

## Introduction

The program `anvi-pan-genome` is the main entry to the pangenomic workflow. Next chapters will demonstrate how it is used, but first, a brief understanding of what it does.

You can run `anvi-pan-genome` on *external genomes*, *internal genomes*, or using a mixture of the two. The difference between external and internal genomes in the context of anvi'o workflow is simple:

* **External genomes**: Anything you have in a FASTA file format (i.e., a genome you downloaded from NCBI, or any other resource).

* **Internal genomes**: Any genome bin you stored in an anvi'o collection (you can create an anvi'o collection by importing automatic binning results, or by using anvi'o to do human guided binning; if you would like to lear more about this please take a look at [the metagenomic workflow]({% post_url anvio/2016-07-22-anvio-tutorial-v2 %})). 

File formats for external genome and internal genome descriptions differ slightly. This is an example `--external-genomes` file:

|name|contigs_db_path|
|:--|:--|
|Name_01|/path/to/contigs-01.db|
|Name_02|/path/to/contigs-02.db|
|Name_03|/path/to/contigs-03.db|
|(...)|(...)|

The only thing you need is to define a `name` for your genome, and point the `contigs_db_path` for it (see the later chapters to see how can you get from you FASTA files to here).

The file format to describe internal genomes of interest is slightly more complex, but allows you to do quite a lot. Here is an example file for `--internal-genomes`:

|name|bin_id|collection_id|profile_db_path|contigs_db_path|
|:--|:--:|:--:|:--|:--|
|Name_01|Bin_id_01|Collection_A|/path/to/profile.db|/path/to/contigs.db|
|Name_02|Bin_id_02|Collection_A|/path/to/profile.db|/path/to/contigs.db|
|Name_03|Bin_id_03|Collection_B|/path/to/another_profile.db|/path/to/another/contigs.db|
|(...)|(...)|(...)|(...)|(...)|

Any genome bin in any anvi'o collection can be included in the analysis (and you can combine these genome bins with cultivar genomes you have downloaded from elsewhere). Let's talk about the paramters, and output files before we start with examples.

## Parameters

When you run `anvi-pan-genome` with `--external-genomes`, or `--internal-genomes`, or both paramters, the program will

* Use [DIAMOND](http://ab.inf.uni-tuebingen.de/software/diamond/) ([Buchnfink et al., 2015](http://www.nature.com/nmeth/journal/v12/n1/abs/nmeth.3176.html)) in 'fast' mode (you can ask DIAMOND to be 'sensitive' by using the flag `--sensitive`) by default to calculate similarities of each protein in every genome against every other protein (which clearly requires you to have DIAMOND installed). Alternatively you could use the flag `--use-ncbi-blast` to use NCBI's `blastp` for protein search.

{:.notice}
***A note from Meren***: I strongly suggest you to do your analysis with the `--use-ncbi-blast` flag. Yes, DIAMOND is very fast, and it may take 8-9 hours to analyze 20-30 genomes with `blastp`. But dramatic increases in speed *rarely* comes without major trade-offs in sensitivity and accuracy, and some of my observations tell me that DIAMOND is not one of those *rare* instances. This clearly deserves a more elaborate discussion, and maybe I will have a chance to write it later, but for now take this as a friendly reminder.

* Use every gene call, whether they are complete or not. Although this is not a big concern for complete genomes, metagenome-assembled genomes (MAGs) will have many incomplete gene calls at the end and at the beginning of contigs. They shouldn't cause much issues, but if you want to exclude them, you can use the `--exclude-partial-gene-calls` flag.

* Use the *maxbit heuristic* that was originally implemented in ITEP ([Benedict et al, 2014](http://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-8)) to eliminate weak matches between two protein seqeunces (you see, the pangenomic workflow first identifies proteins that are somewhat similar by doing similarity searches, and then resolves clusters based on those similarities. In this scenario, weak similarities can connect protein clusters that should not be connected). Here is a definition of maxbit: If you have two protein sequences `A` and `B`, the maxbit is defined as `BITSCORE(A, B) / MIN(BITSCORE(A, A), BITSCORE(B, B))`. Therefore, the maxbit score of two sequences goes to `1` if they are very similar over their entire length, and it goes to `0` if they match over a very short stretch compared to their entire length (regardless of how similar they are at the aligned region), or their sequence identity is low. The default maxbit is `0.5`, but you can change it using the parameter `--maxbit`.

* Use the [MCL](http://micans.org/mcl/) algorithm ([van Dongen and Abreu-Goodger, 2012](http://www.ncbi.nlm.nih.gov/pubmed/22144159)) to identify clusters in protein similarity search results. We use `2` as the *MCL inflation parameter* by default. This parameter defines the sensitivity of the algorithm during the identification of the protein clusters. More sensitivity means more clusters, but of course more clusters does not mean better inference of evolutionary relationships. More information on this parameter and it's effect on cluster granularity is here [http://micans.org/mcl/man/mclfaq.html#faq7.2](http://micans.org/mcl/man/mclfaq.html#faq7.2), but clearly, we, the metagenomics people will need to talk much more about this.

* Utilize every protein cluster, even if they occur in only one genome in your analyses. Of course, the importance of singletons or doubletons will depend on the number of genomes in your analysis, or the qeustion you have in mind. However, if you would like to define a cut-off, you can use the parameter `--min-occurrence`, which is 1, by default. Increasing this cut-off will improve the clustering speed and make the visualization much more manageable, but again, this parameter should be considered in the context of each study.

* Use only one CPU. If you have multiple cores, you should consider increasing that number using the `--num-threads` parameter. The more threads, the merrier the time required for computation.

* Utilze previous search results if there is already a directory. This way you can play with the `--maxbit`, `--mcl-inflation`, or `--min-occurrence` parameters without having to re-do the protein search. However, if you have changed the input protein sequences, either you need to remove the output directory, or use the `--overwrite-output-destinations` flag to redo the search.

{:.notice}
You need another parameter? Well, of course you do! Let us know, and let's have a discussion. We love paramters.

Once `anvi-pan-genome` is done, it will generate an output directory that will be called `pan-output` by default (and you should change it to a more meaningful one using the `--output-dir` paramter). This directory will include essential files such as search results and protein clusters, and another special directory with all the necessary files to invoke an anvi'o interactive interface using this command:

{% highlight bash %}
anvi-interactive -s samples.db -p profile.db -t tree.txt -d view_data.txt -A additional_view_data.txt --manual --title "My Pangenome"
{% endhighlight %}

Let's see some examples.
## Case I: *All I have is a bunch of FASTA files*

Assume you have three FASTA files for three *E. coli* strains you want to analyze with `anvi-pan-genome` (you can download these genomes and play with them if you click those links):

* [E_coli_BL21.fa.gz](https://github.com/meren/anvio/raw/master/tests/sandbox/anvi_pangenome_files/E_coli_BL21.fa.gz)
* [E_coli_B_REL606.fa.gz](https://github.com/meren/anvio/blob/master/tests/sandbox/anvi_pangenome_files/E_coli_B_REL606.fa.gz)
* [E_coli_O111.fa.gz](https://github.com/meren/anvio/blob/master/tests/sandbox/anvi_pangenome_files/E_coli_O111.fa.gz)

Let's assume you downloaded each one of them, and now your work directory looks like this:

{% highlight bash %}
$ ls -l 
-rw-r--r--   1 meren  staff  1362351 Jun 29 14:42 E_coli_BL21.fa.gz
-rw-r--r--   1 meren  staff  1383249 Jun 29 14:43 E_coli_B_REL606.fa.gz
-rw-r--r--   1 meren  staff  1604922 Jun 29 14:43 E_coli_O111.fa.gz
{% endhighlight %}

The first thing is to generate an anvi'o contigs database for each one of them, and run HMMs on them to identify bacterial single-copy core genes. Here is an example for one genome:

{% highlight bash %}
$ gzip -d E_coli_BL21.fa.gz
$ anvi-gen-contigs-database -i E_coli_BL21.fa -o E_coli_BL21.db --split-length -1
$ anvi-run-hmms -c E_coli_BL21.db
{% endhighlight %}

{:.notice}
If you haven't created an anvi'o contigs database before, please read [this section from the metagenomics tutorial]({% post_url anvio/2016-07-22-anvio-tutorial-v2 %}/#creating-an-anvio-contigs-database) to learn the details of this step.

Doing the same thing above for all genomes in a serial manner would have looked like this, in case you are a work-in-progress BASH scripting guru:

{% highlight bash %}
for g in E_coli_BL21 E_coli_B_REL606 E_coli_O111
do
    gzip -d $g.fa.gz
    anvi-gen-contigs-database -f $g.fa -o $g.db -L -1
    anvi-run-hmms -c $g.db
done
{% endhighlight %}

At the end of this, my directory looks like this:

{% highlight bash %}
-rw-r--r--  1 meren  staff  6778880 Jun 29 14:52 E_coli_BL21.db
-rw-r--r--  1 meren  staff  4624090 Jun 29 14:45 E_coli_BL21.fa
-rw-r--r--  1 meren  staff  4567241 Jun 29 14:46 E_coli_BL21.h5
-rw-r--r--  1 meren  staff  6904832 Jun 29 14:52 E_coli_B_REL606.db
-rw-r--r--  1 meren  staff  4695964 Jun 29 14:46 E_coli_B_REL606.fa
-rw-r--r--  1 meren  staff  4638100 Jun 29 14:48 E_coli_B_REL606.h5
-rw-r--r--  1 meren  staff  7996416 Jun 29 14:53 E_coli_O111.db
-rw-r--r--  1 meren  staff  5447818 Jun 29 14:48 E_coli_O111.fa
-rw-r--r--  1 meren  staff  5379365 Jun 29 14:49 E_coli_O111.h5
{% endhighlight %}

Now all anvi'o contigs databases are ready, it is time to prepare the `--external-genomes` file:

|name|contigs_db_path|
|:--|:--|
|BL21|E-coli-BL21.db|
|B_REL606|E-coli-B_REL606.db|
|O111|E-coli-O111.db|

And run `anvi-pan-genome`:

{% highlight bash %}
$ anvi-pan-genome --external external_genomes.txt -o pan-output --num-threads 20
{% endhighlight %}

At the end of this, I can run the `anvi-interactive` to see my results:

{% highlight bash %}
$ cd pan-output/pan-output
$ anvi-interactive -s samples.db -p profile.db -t tree.txt -d view_data.txt -A additional_view_data.txt --manual --title "Pangenome of three E. coli's"
{% endhighlight %}

And this is what pops up in my browser:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2015-11-14-pan-genomics/e-coli-pan.png"><img src="{{ site.url }}/images/anvio/2015-11-14-pan-genomics/e-coli-pan.png" width="80%" /></a>
</div>

This doesn't look too fancy, but your genomes will look be much more informative and interesting than three random *E. coli*'s. Additionally, the output directory will contain two files, `additional_view_data.txt` and `anvio-samples-order.txt`. You can edit these fies to add more contextual information about your genomes, and generate an updated `samples.db` for your interactive interface ([more information on anvi'o samples databases]({% post_url anvio/2015-11-10-samples-db %})). For instance, here is a slightly more rich analysis one with 30 genomes:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2015-11-14-pan-genomics/pan-genome-2.png"><img src="{{ site.url }}/images/anvio/2015-11-14-pan-genomics/pan-genome-2.png" width="80%" /></a>
</div>

## Case II: *I want to work with my genome bins in anvi'o collections*

If you have been doing your metagenomic analyses using anvi'o, it is much more easier to run the pangenomic workflow. All you need to do is to create a `--internal-genomes` file like this one,

|name|bin_id|collection_id|profile_db_path|contigs_db_path|
|:--|:--:|:--:|:--|:--|
|MAG_01|DWH_SAND_Bin_01|SAND|DWH-SAND-MERGED/PROFILE.db|SAND-CONTIGS.db|
|MAG_02|DWH_SAND_Bin_04|SAND|DWH-SAND-MERGED/PROFILE.db|SAND-CONTIGS.db|
|MAG_03|SML_SEWAGE_Bin_08|SEWAGE|SML-SEWAGE-MERGED/PROFILE.db|SEWAGE-CONTIGS.db|

And run `anvi-pan-genome` using it instead:

{% highlight bash %}
$ anvi-pan-genome --internal internal_genomes.txt -o pan-output --num-threads 20
{% endhighlight %}

The rest is the same!

## Case III: Case I + Case II (*I want it all, because science*)

You are lucky! Just mix your internal and external genomes in the same command line:

{% highlight bash %}
$ anvi-pan-genome --internal internal_genomes.txt --external external_genomes.txt -o pan-output --num-threads 20
{% endhighlight %}


## Final words

We are realizing that there is a lot to explore in this front, and excited to work with the community. Please let us know if you have any suggestions using our [discussion group]({% post_url anvio/2015-10-14-anvio-discussion-group %}), the comments section below. You can also keep an eye on our [public code repository](http://github.com/meren/anvio) for new releases.
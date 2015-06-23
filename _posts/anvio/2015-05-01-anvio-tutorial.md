---
layout: post
title: "Anvi'o User Manual"
excerpt: "Sweets for people who managed to install the platform."
modified: 2015-05-01 
tags: []
categories: [anvio]
comments: true
---

{% include _toc.html %}

This article gives a brief overview of anvi'o metagenomic worklow. If you run into any issues, please post a comment down below, or open an <a href="https://github.com/meren/anvio/issues">issue</a>.

If you are here, you must have already installed the platform, and have [run the "mini test"]({% post_url anvio/2015-05-01-installation %}) succesfully.


# Preparation

To run anvio, you need these files (_at least_):

* __A FASTA file__ of your contigs. I will call it `contigs.fa` throughout this manual.
* And __BAM files__ for your samples. Let's say, for the sake of brevity, you have two samples in your analysis, `X` and `Y`, and the BAM files for these samples are named `X-raw.bam` and `Y-raw.bam`.

Contig names in `contigs.fa` must match names found in your BAM files. This is very important. To make sure that is the case, please run these two commands in your terminal:

{% highlight bash %}
grep '>' contigs.fa | head
anvi-profile -i X-raw.bam --list-contigs | head
{% endhighlight %}

Do they look identical? They better be. If they don't you need to fix that before going forward.

If you exported your FASTA file and BAM files using CLC, type these two commands, and you are going to be fine:

{% highlight bash %}
sed -i 's/ .*$//g' contigs.fa
sed -i 's/_contig_/contig/g' contigs.fa
{% endhighlight %}


# Programs to analyze contigs FASTA

Following sections will describe each program you will frequently use while walking you through an example analysis. To see a complete list of programs that are distributed with the platform, type `anvi-` in your terminal, and then press the `TAB` key twice.

## anvi-gen-annotation-db

Annotation database will be the essential ingredient everything you will do with anvio. It will process `contigs.fa`, and store it in a better formatted way. You can decorate your annotation database, but the one command you have to run is this:

    anvi-gen-annotation-database -f contigs.fa -o annotation.db

Each anvi-* program has a help menu that explains available parameters. Don't forget to check them for everything you run. If something is not clearly explained, let me know:

    anvi-gen-annotation-database --help

Once you have the annotation database, you can decorate it with stuff that would be very useful later.

## anvi-populate-genes-table

The program makes populates relevant tables in the annotation database with information on annotation. I will give an example using RAST annotation.

If you have MyRAST installed, you can run these two commands to store the annotation of your contigs in the annotation database (the first line will query RAST server, which may take a while depending on the number of contigs you have), the second line will incorporate the returning info into anvio's annotation database:

{% highlight bash %}
svr_assign_to_dna_using_figfams < contigs.fa > svr_assign_to_dna_using_figfams.txt 
anvi-populate-genes-table annotation.db -p myrast_cmdline_dont_use -i svr_assign_to_dna_using_figfams.txt
{% endhighlight %}

Once you have RAST annotations, I suggest you to run the following command to export the information from the table into a more native matrix form and store it separately (if you take a look at the exported matrix, you will see that it is a simpler, and more standard form of what we got from RAST):

    anvi-export-genes-table annotation.db -o rast_annotation.txt

If you have to re-create the annotation database for some reason, you can now use this newly generated file (`rast_annotation.txt`) to repopulate the annotation database instead of querying RAST again using the following command:

    anvi-populate-genes-table annotation.db -p default_matrix -i rast_annotation.txt

## anvi-populate-search-table

The program populates tables that holds HMM search results in the annotation database.

Anvi'o can do wonders with HMM models. To decorate your annotation database with hits from HMM models that ship with the platform (which, at this point, constitute published single-copy gene collections), run this command:

    anvi-populate-search-table annotation.db


## anvi-populate-collections-table

This program populates tables that holds information about different organizations of contigs (collection of bins).

For instance, if you have some CONCOCT results for your merged profile (which you DON'T have at this point if you are reading this tutorial the first time, but you will have them later if you keep reading), you can add them into the annotation database like this:

    anvi-populate-collections-table annotation.db --parser concoct -i concoct.txt

In this case `concoct.txt` is the clustering.csv file CONCOCT generates. To see the format you can take a look at [`anvio/tests/mini_test/concoct.txt`](https://github.com/meren/anvio/blob/master/tests/mini_test/concoct.txt)


# Programs to analyze BAM files

## anvi-init-bam

anvio likes BAM files when they are sorted and indexed. This is why I named your BAM file for sample X you got from CLC or Bowtie as `X-raw.bam` instead of `X.bam`. To initialize this BAM file you can run this command:

    anvi-init-bam X-raw.bam -o X

But of course it is not fun to do every BAM file you have one by one. So what to do?

A slightly better way to do would require you to do it in a `for` loop. First create a file called, say, `SAMPLE_IDs`. For your samples (`X` and `Y`) it will look like this:

{% highlight bash %}
$ cat SAMPLE_IDs
X
Y
{% endhighlight %}

Then, you can run `anvi-init-bam` on all of them by typing this:

{% highlight bash %}
for sample in `cat SAMPLE_IDs`; do anvi-init-bam $sample-raw.bam -o $sample; done
{% endhighlight %}

Good.

## anvi-profile

Profiling step makes sense of each BAM file separately by utilizing the information stored in the annotation database. The result of the profiling step is a special file that describes the run (`RUNONFO.cp`), and a profile database (`PROFILE.db`).

The minimal command to profile a BAM file looks like this:

anvi-profile -i X.bam -a annotation.db

But I encourage you to take a look at the default paramers. One of the most critical parameter is `-M` (`--min-contig-length`) parameter. The default is 10,000. Which means the profiling step will take into consideration only the contigs that are longer than 10Kb. This may be too large for your analysis. But clustering and visualization steps in anvio have some limitations, so you can't really say `-M 0` in most cases. The rule of thumb is to keep the number of contigs anvio will deal to a maximum of 20,000. How can you know how many contigs are there for a given `-M` value? Well, one thing to find that out is this:

    sqlite3 annotation.db 'select count(*) from contigs_basic_info where length > 10000;'

This command will print out the number of contigs longer than 10Kb in your dataset. You can try different values until the output is about 20,000, and use that value for `-M`. But I will not recommend you to go below 1Kb. The main reason to that is the fact that anvio relies on k-mer frequencies to better cluster contigs, and tetra-nucleotides (the default way for anvio to make sense of the sequence composotion) become very unstable very quickly.

Once you know what you `-M` is, you can, again, profile multiple samples using the similar approach we used for initializing BAM files:

    for sample in `cat SAMPLE_IDs`; do anvi-profile -i $sample.bam -M YOUR_M_VALUE -a annotation.db; done

__Note__: If you are planning to work with and visualize single profiles (without merging), use `--cluster-contigs` flag (more explanation on this will come later). You don't need to use this flag if you are planning to merge multiple profiles (i.e., if you have more than one BAM files to work with).


# Programs to work with anvi'o profiles

## anvi-merge

Once you have your BAM files profiled, the next logical step is to merge all the profiles that should be analyzed together.

It goes without saying that every profiling step must have used the same parameters for analysis. If profiles have been generated with different annotation databases or with different parameters will not get merged, and you will get angry error messages from anvio.

In an ideal case, this should be enough to merge your stuff:

    anvi-merge X/RUNINFO.cp Y/RUNINFO.cp -o XY-MERGED -a annotation.db

Or alternatively you can run this:

    anvi-merge */RUNINFO.cp -o XY-MERGED -a annotation.db

Before you run these commands in your real-world data, you must understand the details of clustering.

### Clustering during merging

A default step in the merging process is to generate a hierarchical clustering of all splits using anvi'o's default _clustering configurations_. One of these clustering configurations clusters contigs using only k-mer frequencies (if you generated your annotation database with default parameters, the k is 4, and your k-mer frequencies will be 'tetranucleotide frequencies'), another one of them mixes k-mer frequencies with distribution patterns across samples for clustering, etc. The hierarchical clustering result is necessary for __visualization__, and __supervised binning__. Therefore, by default, anvi'o will attempt to cluster your contigs using these configurations. However, if you have, say, more than 25,000 splits, clustering step will be very time consuming (multiple hours to even days), and visualization of this data will be very challenging. There are other solutions to this problem that will be discussed later, but if you would like to skip hierarchical clustering, you can use `--skip-hierarchical-clustering` flag.

During merging, anvi'o will also use [CONCOCT](http://www.nature.com/nmeth/journal/v11/n11/full/nmeth.3103.html) for unsupervised binning. CONCOCT can deal with hundreds of thousands of splits. Which means, regardless of the number of splits you have, and even if you skip the hierarchical clustering step, there will be a collection in the merged profile database (collection id of which will be 'CONCOCT') with genome bins identified by CONCOCT in an unsupervised manner, from which you can generate a summary. But if you would like to skip CONCOCT clustering, you can use `--skip-concoct-binning` flag. 

## anvi-interactive

Interactive interface is good for a couple of things: it allows you to browse your data in an intuitive way, shows you multiple aspects of your data, visualize the results of unsupervised binning, and most importantly, it allows you to perform supervised binning. In order to run the interactive interface on a run, hierarchical clustering must have been performed, and stored in the profile database (otherwise you will hear anvi'o complaining about that).

 Once the merging is done you can finally run the interactive interface:

    anvi-interactive -p XY-MERGED/PROFILE.db -a annotation.db

If you had too many contigs and had to skip the hierarchical clustering step during merging, and if you only have unsupervised binning done by CONCOCT, don't be worried and keep reading on.

But just to give you an idea, here is a screenshot from the interactive interface (taken while we were analyzing a random interface with anvi'o version 0.8.1):

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/misc/infant-gut-screenshot.png"><img src="{{ site.url }}/images/anvio/misc/infant-gut-screenshot.png" width="50%" /></a>
</div>

## anvi-summarize

Once you have a collection, it is time for you to summarize your results.

Collection could be genome bins identified by CONCOCT in an unsupervised manner during anvi-merge, or it could be genome bins you identified through anvi-interactive. anvi-summary will take the collection id as a parameter, and generate a static HTML output for you, during which, among other things,

* All your splits will be merged back to contigs and stored as FASTA files,
* Completion and contamination scores for your genome bins will be computed and stored into files,
* Matrix files will be generated for your genome bins across your samples with respect to their mean coverage, variability, etc.

You can run this to summarize the bins saved under the collection id 'CONCOCT' into MY_SUMMARY directory:

    anvi-summarize -p XY-MERGED/PROFILE.db -a annotation.db -o XY-MERGED-SUMMARY -c CONCOCT

If you are not sure which collections are available to you, you can always see a list of them by running this command:

    anvi-summarize -p XY-MERGED/PROFILE.db -a annotation.db -o XY-MERGED-SUMMARY --list-collections


## anvi-refine

After running anvi-summarize, you may realize that you are  not happy with one or more of your bins. This often is the case when you are working with very large datasets and when you are forced to skip the supervised clustering step. anvi-summarize gives you the ability to make finer adjustements in a bin that is not complete.

Please read [this article]({% post_url anvio/2015-05-11-anvi-refine %}) for a comprehensive introduction to refinement capacity of anvi'o.

---

# FAQ

### I run into an issue, who's fault is it?

It is Tom's. But you can always enter a [bug report](https://github.com/meren/anvio/issues), and it would be very useful for us to remember to fix that. The codebase is young, and there will be issues.

### I have suggestions!

This is great. That's exactly why you have accesss to this codebase. Either file an [issue](https://github.com/meren/anvio/issues), or send an e-mail (a.murat.eren@gmail.com).


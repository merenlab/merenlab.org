---
layout: post
title: "Combining metagenomics with metatranscriptomics"
excerpt: "Tricks for people who like to go deeper."
modified: 2015-06-10
tags: []
categories: [anvio]
comments: true
author: meren
---

{% include _toc.html %}

[The anvi'o metagenomics workflow]({% post_url anvio/2015-05-02-anvio-tutorial %}) is quite straighforward. But what if you have multiple sources of data that are mapped to the same contigs? If you have metagenomic, metatranscriptomic data from the same environment, and if you want to do a combined analysis, continue reading.

_.. and please do not hesitate to get in touch if this recipe does not work for you .._

## Context

Just to make things a bit easier, lets assume this mock example: You have two environments, and for each environment you have metagenomic and metatranscriptomic data.

You co-assembled your metagenomic data, and get your contigs. Then, you mapped your short reads back to these contigs from both metagenomic and metatranscriptomic samples. As a result of which, you have four BAM files:

* Env-A-Metagenome.bam
* Env-A-Metatranscriptome.bam
* Env-B-Metagenome.bam
* Env-B-Metatranscriptome.bam

You created your contigs database, and you profiled each of your BAM file  the way it is described in the standard metagenomic workflow. So far so good.

And it is time to merge.

## Merging without clustering

Normally you would have merged all your profiles by typing,

{% highlight bash %}
anvi-merge */RUNINFO.cp -c contigs.db -o MERGED/
{% endhighlight %}

Merging process creates multiple clusterings. Two of the default clusterings takes into account the distribution of each contig across samples. This means, if you merge your metagenome and metatranscriptomes, the clustering of each split will not only depend on their coverage in metagenomic samples, but also in metatranscriptomes. Which is not the most appropriate way to cluster contigs if we want to see an organization that reflect genome bins.

Therefore, we run our merging without any clustering, by adding these two parameters:

{% highlight bash %}
anvi-merge */RUNINFO.cp -c contigs.db -o MERGED \
           --skip-hierarchical-clustering \
           --skip-concoct-binning
{% endhighlight %}

Please see the help menu to see the details of any flag you are not familiar with.

This will result in a merged profile (which will be in the `MERGED/` directory), that will not contain any clustering. Therefore, `anvi-interactive` will be useless on it. But don't despair!

## Creating a precise clustering configuration

So, our goal is to do a clustering analysis only with respect to the metagenomic samples, as we don't want the coverage of our contigs in metatranscriptomic samples to affect their organization (I assume you can see why we don't want that).

Fortunately, it is easy to perform clustering with greater precision than just using all samples. To completely understand this you need to read about _clustering configurations_ concept we use in anvi'o. But let's move on.

This file is the default clustering configuration anvi'o uses for clustering contigs in merged samples based on sequence composition and coverage:

{% highlight ini %}
[general]

[TNF !CONTIGS.db::kmer_contigs]

[Coverage PROFILE.db::mean_coverage_contigs]

[Coverage_N PROFILE.db::mean_coverage_contigs]
normalize = False
log = True

[GC_content !CONTIGS.db::splits_basic_info]
columns_to_use = gc_content_parent,gc_content_parent,gc_content_parent
normalize = False
{% endhighlight %}

I will create a copy of this file in my work directory. You can copy-paste the text [from here](https://github.com/meren/anvio/blob/master/anvio/data/clusterconfigs/merged/tnf-cov) and put it in a file (lets call it `my_cluster_config.ini`), or you can use this command line to create a copy form the one you have on your disk right now:

{% highlight bash %}
cp `python -c 'import anvio; import os; \
    print os.path.dirname(anvio.__file__)'`/data/clusterconfigs/merged/tnf-cov my_cluster_config.ini
{% endhighlight %}

This config file creates a complex matrix, where there are multiple sources of information is mixed for each split. I will slightly change this clustering configuration so it employs the coverage information _only_ from the metagenomic mapping results.

To achieve this, I need to define a new `columns_to_use` directive with the proper sample names under `Coverage` and `Coverage_N` sections.

Well, I know my sample names, right? They must be "Env-A-Metagenome", and "Env-B-Metagenome". They may be. But it is **always** the best practice to make sure whether our sample names have been changed by anvi'o (yes, anvi'o _can_ change your sample names to make things more compatible with proper practices). Here is an easy way to see how does my sample names look like in the merged profile database:

{% highlight bash %}
$ sqlite3 OIL-PLUME-MERGED/PROFILE.db \
      'select value from self where key = "samples";' | sed 's/,/\n/g'
Env_A_Metagenome
Env_A_Metatranscriptome
Env_B_Metagenome
Env_B_Metatranscriptome
{% endhighlight %}

As you can see, anvi'o did change every `-` character with a `_` character. Because I should have never used `-` to name my files in the first place! Well done, anvi'o, thank you very much.

Now, I really know which sample names are relevant, so I can edit my fresh cluster configuration copy, `my_cluster_config.ini`, to reflect my desire to use only metagenomic samples for clustering:

{% highlight ini %}
[general]

[TNF !CONTIGS.db::kmer_contigs]

[Coverage PROFILE.db::mean_coverage_contigs]
columns_to_use = Env_A_Metagenome,Env_B_Metagenome

[Coverage_N PROFILE.db::mean_coverage_contigs]
columns_to_use = Env_A_Metagenome,Env_B_Metagenome
normalize = False
log = True

[GC_content !CONTIGS.db::splits_basic_info]
columns_to_use = gc_content_parent,gc_content_parent,gc_content_parent
normalize = False
{% endhighlight %}

## Clustering

Now I have a proper clustering configuration. It is time to let anvi'o do the clustering. For this, we will use a commandline program: `anvi-experimental-organization`

And here how I run it:

{% highlight bash %}
anvi-experimental-organization my_cluster_config.ini \
                               -i MERGED/            \
                               -c contigs.db      \
                               -p MERGED/PROFILE.db  \
                               -o new_newick_tree.txt
{% endhighlight %}


I know that this process can take a while depending on the number of splits I have (that's why it is always a good idea to run this in a screen session on a server somewhere instead of doing it on my laptop):

{% highlight bash %}
$ sqlite3 MERGED/PROFILE.db 'select value from self where key="num_splits";'
    22155
{% endhighlight %}

22,155... Sigh. Right at the edge. As you know, we never really want to go beyond 20,000 splits with the hirarchical clustering in anvi'o. CONCOCT is the only way to go if you have more splits than 20,000 - 25,000. There shall be another post to describe that in the future, but meanwhile, you set your `-M` during profiling appropriately to avoid too many splits.

## Calling interactive interface with the new clustering

The result of `anvi-experimental-organization` is this new file: `new_newick_tree.txt`, which contains the proper clustering of my splits that can help me visualize my merged run.

Here is how I start the interactive interface using this file:

{% highlight bash %}
anvi-interactive -c contigs.db -p MERGED/PROFILE.db -t new_newick_tree.txt
{% endhighlight %}

Voilà!-c
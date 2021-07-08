---
layout: post
title: "Combining metagenomics with metatranscriptomics"
excerpt: "Tricks for people who like to go deeper."
modified: 2015-06-10
tags: []
categories: [anvio]
comments: true
authors: [meren]
---


[The anvi'o metagenomics workflow]({% post_url anvio/2016-06-22-anvio-tutorial-v2 %}) is quite straighforward. But what if you have multiple sources of data that are mapped to the same contigs? If you have metagenomic, metatranscriptomic data from the same environment, and if you want to do a combined analysis, continue reading.

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
anvi-merge */PROFILE.db -c CONTIGS.db -o MERGED/
{% endhighlight %}

Merging process creates multiple clusterings. Two of the default clusterings takes into account the distribution of each contig across samples. This means, if you merge your metagenome and metatranscriptomes, the clustering of each split will not only depend on their coverage in metagenomic samples, but also in metatranscriptomes. Which is not the most appropriate way to cluster contigs if we want to see an organization that reflect genome bins.

Therefore, we run our merging without any clustering, by adding these two parameters:

{% highlight bash %}
anvi-merge */PROFILE.db -c CONTIGS.db -o MERGED \
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

Well, I know my sample names, right? They must be "Env-A-Metagenome", and "Env-B-Metagenome". They may be. But it is **always** the best practice to make sure whether our sample names have been changed by anvi'o (yes, anvi'o _can_ change your sample names (i.e., by replaceing `-` characters with `_` characters, or by replacing ` ` with `_`) to make things more compatible with proper practices). Here is an easy way to see how does my sample names look like in the merged profile database:

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
                               -c CONTIGS.db      \
                               -p MERGED/PROFILE.db  \
                               --skip-store-in-db \
                               -o new_newick_tree.txt
{% endhighlight %}

{:.notice}
You can also instruct `anvi-experimental-organization` to store the resulting tree directly into the profile database. If you would like that you can remove the `--skip-store-in-db`, and use the flag `--name` to specify under which name the tree should be stored.

This process will take some time depending on the number of splits. Check the number of splits you have:

{% highlight bash %}
$ sqlite3 MERGED/PROFILE.db 'select value from self where key="num_splits";'
    22155
{% endhighlight %}

And make sure you have less than 25,000. If you have more, either you can profile your files with a higher `-M` value, or you can push us by entering an issue to the GitHub page so we can try to sort this out in a better way.

## Using the new clustering in the interactive interface with the new clustering

If you did not use the flag `--skip-store-in-db`, you will find your tree in the combo box in the interface where you always find your trees. But if you did exactly how this tutorial did it, you will find your tree in this `new_newick_tree.txt`, which will contain the proper clustering of splits that can be used to visualize the merged run.

In that case, this is how you can start the interactive interface with this tree file, so it would also appear in the relevant combo box:

{% highlight bash %}
anvi-interactive -c CONTIGS.db -p MERGED/PROFILE.db -t new_newick_tree.txt
{% endhighlight %}

Voilà!

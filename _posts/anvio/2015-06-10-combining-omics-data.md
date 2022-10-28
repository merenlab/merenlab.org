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

If the recipe doesn't work for you, or you need more from it, please feel free to reach out to us:

{% include _join-anvio-discord.html %}

## Context

Just to make things a bit easier, lets assume this mock example: You have two environments, and for each environment you have metagenomic and metatranscriptomic data.

You assembled one of your metagenomes or co-assembled multiple of your metagenomic data, turned them into a single {% include ARTIFACT name="contigs-db" text="contigs database" %}, and mapped your short reads back to these contigs from both metagenomic and metatranscriptomic samples. As a result of which, you have four {% include ARTIFACT name="bam-file" text="bam files" %}:

* Env-A-Metagenome.bam
* Env-A-Metatranscriptome.bam
* Env-B-Metagenome.bam
* Env-B-Metatranscriptome.bam

Then, you used the program {% include PROGRAM name="anvi-profile" %} to profile each one of your BAM files and to get a {% include ARTIFACT name="single-profile-db" %} for each one of them. So far so good?

## Merging without clustering

Normally you would have merged all your single profile databases using {% include PROGRAM name="anvi-merge" %},

``` bash
anvi-merge */PROFILE.db -c CONTIGS.db -o MERGED/
```

Merging process automatically generates multiple clusterings of your contigs to organize them in the anvi'o {% include ARTIFACT name="interactive" %} interface. Two of the default clusterings takes into account the distribution of each contig across samples. This means, if you merge your metagenome and metatranscriptomes, the clustering of your contigs will be determined by both their coverage in metagenomes and metatranscriptomes. Which is not the most appropriate way to cluster contigs if we want to see an organization that reflect genome bins.

Therefore, we run our merging without any clustering, by adding these two parameters:

{% highlight bash %}
anvi-merge */PROFILE.db \
           -c CONTIGS.db \
           -o MERGED \
           --skip-hierarchical-clustering
{% endhighlight %}

Please see the help menu to see the details of any flag you are not familiar with.

This will result in a merged profile (which will be in the `MERGED/` directory), that will not contain any clustering and if you were to try to interactively visualize your merged profile, {% include PROGRAM name="anvi-interactive" %} will complain about it. But don't despair!

## Creating a precise clustering configuration

So, our goal is to do a clustering analysis only with respect to the metagenomic samples, as we don't want the coverage of our contigs in metatranscriptomic samples to affect their organization (I assume you can see why we don't want that).

Fortunately, it is easy to perform clustering with greater precision than just using all samples. Anvi'o uses clustering configurations to organize samples in a context dependent manner, which in return allows users to be extremely specific regarding which data source they would like to include in their clustering of contigs. For instance, you can find all default clustering configurations anvi'o uses [in these directories](https://github.com/merenlab/anvio/tree/master/anvio/data/clusterconfigs) under anvi'o codebase.

The following file, for example, is the default clustering configuration anvi'o uses for clustering contigs in merged samples based on sequence composition and coverage:

```ini
[general]

[TNF !CONTIGS.db::kmer_contigs]

[Coverage PROFILE.db::mean_coverage_contigs]
table_form = view

[Coverage_N PROFILE.db::mean_coverage_Q2Q3_contigs]
table_form = view
normalize = False
log = True

[GC_content !CONTIGS.db::splits_basic_info]
columns_to_use = gc_content_parent,gc_content_parent,gc_content_parent
normalize = False
```

We will use this file as a template to 'speficy' which samples we wish anvi'o to use for clustering.

For that, first I will create a copy of this file in my work directory. You can copy-paste the text [from here](https://github.com/meren/anvio/blob/master/anvio/data/clusterconfigs/merged/tnf-cov) and put it in a file (lets call it `my_cluster_config.ini`), or you can use this command line to create a copy form the one you have on your disk right now:

```
cp `python -c 'import anvio; import os; \
    print(os.path.dirname(anvio.__file__))'`/data/clusterconfigs/merged/tnf-cov my_cluster_config.ini
```

This config file creates a complex matrix, where there are multiple sources of information is mixed for each split. I will slightly change this clustering configuration so it employs the coverage information _only_ from the metagenomic mapping results.

To achieve this, I need to define a new `columns_to_use` directive with the proper sample names under `Coverage` and `Coverage_N` sections.

As far as this example goes, I can assume that I know my sample names, right? They must be `Env-A-Metagenome`, and `Env-B-Metagenome`. They certainly may be. But it is **always** the best practice to make sure whether our sample names have been changed by anvi'o (yes, anvi'o _can_ change your sample names (i.e., by replaceing `-` characters with `_` characters, or by replacing ` ` with `_`) to make things more compatible with proper practices). The best way to learn what anvi'o thinks your sample names are is to use the program {% include PROGRAM name="anvi-db-info" %}. When you run this program on any merged profile database, you will learn the precise sample names it contains from the output. Here is an example output for this mock example:

```
anvi-db-info MERGED/PROFILE.db

DB Info (no touch)
===============================================
Database Path ................................: MERGED/PROFILE.db
description ..................................: [Found, and it is 2020 characters long]
db_type ......................................: profile (variant: None)
version ......................................: 38


DB Info (no touch also)
===============================================
(...)
samples ......................................: Env_A_Metagenome, Env_A_Metatranscriptome, Env_B_Metagenome, Env_B_Metatranscriptome
(...)
```

As you can see, anvi'o did change every `-` character with a `_` character. Because I should have never used `-` to name my files in the first place! Well done, anvi'o, thank you very much.

Now, I really know which sample names are relevant, so I can edit my fresh cluster configuration copy, `my_cluster_config.ini`, to reflect my desire to use only metagenomic samples for clustering:

```ini
[general]

[TNF !CONTIGS.db::kmer_contigs]

[Coverage PROFILE.db::mean_coverage_contigs]
table_form = view
columns_to_use = Env_A_Metagenome,Env_B_Metagenome

[Coverage_N PROFILE.db::mean_coverage_Q2Q3_contigs]
table_form = view
columns_to_use = Env_A_Metagenome,Env_B_Metagenome
normalize = False
log = True

[GC_content !CONTIGS.db::splits_basic_info]
columns_to_use = gc_content_parent,gc_content_parent,gc_content_parent
normalize = False
```

## Clustering

Now I have a proper clustering configuration. It is time to let anvi'o do the clustering based on this recipe using the program {% include PROGRAM name="anvi-experimental-organization" %}.

This is how one could run to add a new clustering to the database using the current recipe:

```
anvi-experimental-organization my_cluster_config.ini \
                               --input-directory MERGED/ \
                               --contigs-db CONTIGS.db \
                               --profile-db MERGED/PROFILE.db \
                               --name CLUSTERING_BASED_ON_MGX
```

Now, when you run the program {% include PROGRAM name="anvi-interactive" %}, you will find your new clustering under the combo box called "items order", you can select it and re-draw your data to see an updated organization of your contigs.

Alternatively you can instruct {% include PROGRAM name="anvi-experimental-organization" %} to **not** store the resulting tree directly into the profile database, and instead report a newick-formatted clustering {% include ARTIFACT name="dendrogram" %}:

```
anvi-experimental-organization my_cluster_config.ini \
                               --input-directory MERGED/ \
                               --contigs-db CONTIGS.db \
                               --profile-db MERGED/PROFILE.db \
                               --skip-store-in-db \
                               --output-file CLUSTERING_BASED_ON_MGX.newick
```

Please note that this process may take a long time depending on the number of splits you have in your profile database. You can learn that number from the output of {% include PROGRAM name="anvi-db-info" %} or by running this on your terminal:

``` bash
$ sqlite3 MERGED/PROFILE.db 'select value from self where key="num_splits";'
    22155
```

If you have more than 25,000, you are in trouble. Ask us on Discord how to solve this problem:

{% include _join-anvio-discord.html %}

## Using the new clustering in the interactive interface with the new clustering

If you did not use the flag `--skip-store-in-db`, you will find your tree in the combo box in the interface where you always find your trees. But if you did exactly how this tutorial did it, you will find your tree in this `new_newick_tree.txt`, which will contain the proper clustering of splits that can be used to visualize the merged run.

In that case, this is how you can start the interactive interface with this tree file, so it would also appear in the relevant combo box:

```
anvi-interactive -c CONTIGS.db \
                 -p MERGED/PROFILE.db \
                 -t CLUSTERING_BASED_ON_MGX.txt
```

Voilà!

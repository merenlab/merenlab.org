---
layout: page
title: anvi-estimate-scg-taxonomy [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-estimate-scg-taxonomy
image:
  featurerelative: ../../../images/header.png
  display: true
---

Estimates taxonomy at genome and metagenome level. This program is the entry point to estimate taxonomy for a given set of contigs (i.e., all contigs in a contigs database, or contigs described in collections as bins). For this, it uses single-copy core gene sequences and the GTDB database.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/qclayssen.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Quentin Clayssen</span><div class="page-author-social-box"><a href="mailto:quentin.clayssen@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/ClayssenQ" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/qclayssen" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[scgs-taxonomy](../../artifacts/scgs-taxonomy) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[collection](../../artifacts/collection) <img src="../../images/icons/COLLECTION.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[bin](../../artifacts/bin) <img src="../../images/icons/BIN.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[metagenomes](../../artifacts/metagenomes) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[genome-taxonomy](../../artifacts/genome-taxonomy) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[genome-taxonomy-txt](../../artifacts/genome-taxonomy-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Usage


This program **quickly makes taxonomy estimates for genomes, metagenomes, or bins stored in your <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>.**

This is the final step in the scg-taxonomy workflow (described in its entirety [here](http://merenlab.org/2019/10/08/anvio-scg-taxonomy/)). Before running this program, you'll need to have run both <span class="artifact-n">[anvi-setup-scg-taxonomy](/software/anvio/help/7.1/programs/anvi-setup-scg-taxonomy)</span> and <span class="artifact-n">[anvi-run-scg-taxonomy](/software/anvio/help/7.1/programs/anvi-run-scg-taxonomy)</span> on the <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> you're working with for this project. 

[This tutorial](http://merenlab.org/2019/10/08/anvio-scg-taxonomy/#estimating-taxonomy-in-the-terminal) also includes a comprehensive overview of what this program can do. See that page for more information on all of the features described below. 

Keep in mind that this uses single-copy core genes and their hits in  [GTDB](https://gtdb.ecogenomic.org/), so it will not work well in bins with low completion or for Eukaryotic organisms. 

This program is implicitly run in the interactive interface, when you turn on "Realtime taxonomy estimation for bins (whenever possible)." So, if you've ever wondered where those estimates were coming from, now you know. 

So, what can this program do?

### 1. Estimate the taxonomy of a single genome

By default, this program wll assume your <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> contains only one genome, and will try to use the single-copy core genes (that were associated with taxonomy when you ran <span class="artifact-n">[anvi-run-scg-taxonomy](/software/anvio/help/7.1/programs/anvi-run-scg-taxonomy)</span>) to try to identify the taxonomy of your genome. 

When you run 

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;scg&#45;taxonomy &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span>
</div>

It will give you the best taxonomy hit for your genome. If you would like to see how it got there (by looking at the hits for each of the single-copy core genes), just use the `--debug` flag to see more information, as so: 

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;scg&#45;taxonomy &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                           &#45;&#45;debug 
</div>

### 2. Estimate the taxa within a metagenome 

By running this program in metagenome mode, it will assume that your <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> contains multiple genomes and will try to give you an overview of the taxa within it. To do this, it will determine which single-copy core gene has the most hits in your contigs (for example `Ribosomal_S6`), and then will look at the taxnomy hits for that gene across your contigs. The output will be this list of taxonomy results. 

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;scg&#45;taxonomy &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                           &#45;&#45;metagenome&#45;mode 
</div>

If you want to look at a specific gene (instead of the one with the most hits), you can also tell it to do that. For example, to tell it to look at Ribosomal_S9, run

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;scg&#45;taxonomy &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                           &#45;&#45;metagenome&#45;mode  \
                           &#45;&#45;scg&#45;name Ribosomal_S9
</div>

### 3. Look at relative abundance across multiple samples 

If you provide a <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span> or <span class="artifact-n">[single-profile-db](/software/anvio/help/7.1/artifacts/single-profile-db)</span>, then you'll be able to look at the relative abundance of your taxonomy hits (through a single-copy core gene) across your samples. Essentially, this adds additional columns to your output (one per sample) that descrbe the relative abundance of each hit in each sample. 

Running this will look something like this, 
<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;scg&#45;taxonomy &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                           &#45;&#45;metagenome&#45;mode  \
                           &#45;&#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7.1/artifacts/profile&#45;db)</span> \
                           &#45;&#45;compute&#45;scg&#45;coverages
</div>

For an example output, take a look at [this page](http://merenlab.org/2019/10/08/anvio-scg-taxonomy/#contigs-db--profile-db).

### 4. Estimate the taxonomy of your bins 

This program basically looks at each of the <span class="artifact-n">[bin](/software/anvio/help/7.1/artifacts/bin)</span>s in your <span class="artifact-n">[collection](/software/anvio/help/7.1/artifacts/collection)</span> as a single genome and tries to assign it taxonomy information. To do this, simply provide a collection, like this:

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;scg&#45;taxonomy &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                           &#45;&#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/7.1/artifacts/collection)</span>
</div>

You can also look at the relative abundances across your samples at the same time, by running something like this: 

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;scg&#45;taxonomy &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                           &#45;&#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/7.1/artifacts/collection)</span>  \
                           &#45;&#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7.1/artifacts/profile&#45;db)</span> \
                           &#45;&#45;compute&#45;scg&#45;coverages
</div>

### 5. Look at multiple metagenomes at the same time

You can even use this program to look at multiple metagenomes by providing a <span class="artifact-n">[metagenomes](/software/anvio/help/7.1/artifacts/metagenomes)</span> artifact. This is useful to get an overview of what kinds of taxa might be in your metagenomes, and what kinds of taxa they share. 

Running this

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;scg&#45;taxonomy &#45;&#45;metagenomes <span class="artifact&#45;n">[metagenomes](/software/anvio/help/7.1/artifacts/metagenomes)</span> \
                           &#45;&#45;output&#45;file&#45;prefix EXAMPLE
</div>

will give you an output file containing all taxonomic levels found and their coverages in each of your metagenomes. 

For a concrete example, check out [this page](http://merenlab.org/2019/10/08/anvio-scg-taxonomy/#many-contigs-dbs-for-many-metagenomes). 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-estimate-scg-taxonomy.md) to update this information.


## Additional Resources


* [Usage examples and warnings](http://merenlab.org/scg-taxonomy)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-estimate-scg-taxonomy) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

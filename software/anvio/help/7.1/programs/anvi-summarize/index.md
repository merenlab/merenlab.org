---
layout: page
title: anvi-summarize [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-summarize
image:
  featurerelative: ../../../images/header.png
  display: true
---

Summarizer for anvi&#x27;o pan or profile db&#x27;s. Essentially, this program takes a collection id along with either a profile database and a contigs database or a pan database and a genomes storage and generates a static HTML output for what is described in a given collection. The output directory will contain almost everything any downstream analysis may need, and can be displayed using a browser without the need for an anvi&#x27;o installation. For this reason alone, reporting summary outputs as supplementary data with publications is a great idea for transparency and reproducibility.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[collection](../../artifacts/collection) <img src="../../images/icons/COLLECTION.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[pan-db](../../artifacts/pan-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[genomes-storage-db](../../artifacts/genomes-storage-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[summary](../../artifacts/summary) <img src="../../images/icons/SUMMARY.png" class="artifact-icon-mini" /></span></p>


## Usage


Anvi-summarize lets you look at a **comprehensive overview of your <span class="artifact-n">[collection](/software/anvio/help/7.1/artifacts/collection)</span>** and its many statistics that anvi'o has calculated. 

It will create a folder called `SUMMARY` that contains many different summary files, including an HTML output that conviently displays them all for you. This folder will contain anything a future user might use to import your collection, so it's useful to send to others or transfer an entire anvi'o collection and all of its data. 

In a little more detail, this program will   
* generate <span class="artifact-n">[fasta](/software/anvio/help/7.1/artifacts/fasta)</span> files containing your original contigs.   
* estimate various stats about each of your bins, including competition, redundacy, and information about all of your <span class="artifact-n">[hmm-hits](/software/anvio/help/7.1/artifacts/hmm-hits)</span>    
* generate various tab-delimited matrix files with information about your bins across your samples, including various statistics.   

## Running anvi-summarize 

### Running on a profile database

A standard run of anvi-summarize on a <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span> will look something like this:

<div class="codeblock" markdown="1">
anvi&#45;summarize &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
               &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7.1/artifacts/profile&#45;db)</span> \
               &#45;o MY_SUMMARY \
               &#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/7.1/artifacts/collection)</span>
</div>

This will name the output directory `MY_SUMMARY` instead of the standard `SUMMARY`. 

When running on a profile database, you also have options to 
* output very accurate (but intensely processed) coverage and detection data for each gene (using `--init-gene-coverages`)
* edit your contig names so that they contain the name of the bin that the contig is in (using `--reformat-contig-names`)
* also display the amino acid sequeunces for your gene calls.  (using `--report-aa-seqs-for-gene-calls`)

### Running on a pan database

When running on a <span class="artifact-n">[pan-db](/software/anvio/help/7.1/artifacts/pan-db)</span>, you'll want to instead provide the associated genomes storage database. 

<div class="codeblock" markdown="1">
anvi&#45;summarize &#45;g <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/7.1/artifacts/genomes&#45;storage&#45;db)</span> \
               &#45;p <span class="artifact&#45;n">[pan&#45;db](/software/anvio/help/7.1/artifacts/pan&#45;db)</span> \
               &#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/7.1/artifacts/collection)</span> 
</div>

You can also choose to display DNA sequences for your gene clusters instead of amino acid sequences with the flag `--report-DNA-sequences`

### Other notes

If you're unsure what collections are in your database, you can run this program with the flag `--list-collections` or by running <span class="artifact-n">[anvi-show-collections-and-bins](/software/anvio/help/7.1/programs/anvi-show-collections-and-bins)</span>.

You can also use the flag `--quick-summary` to get a less comprehensive summary with a much shorter processing time. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-summarize.md) to update this information.


## Additional Resources


* [anvi-summarize in the metagenomic workflow tutorial](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-summarize)

* [anvi-summarize in the pangenomic workflow tutorial](http://merenlab.org/2016/11/08/pangenomics-v2/#summarizing-an-anvio-pan-genome)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-summarize) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

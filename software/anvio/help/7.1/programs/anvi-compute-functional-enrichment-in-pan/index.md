---
layout: page
title: anvi-compute-functional-enrichment-in-pan [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-compute-functional-enrichment-in-pan
image:
  featurerelative: ../../../images/header.png
  display: true
---

A program that computes functional enrichment within a pangenome..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/ivagljiva.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Iva Veseli</span><div class="page-author-social-box"><a href="mailto:iveseli@uchicago.edu" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://github.com/ivagljiva" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/adw96.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Amy D. Willis</span><div class="page-author-social-box"><a href="http://statisticaldiversitylab.com/team/amy-willis" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:adwillis@uw.edu" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/AmyDWillis" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/adw96" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[misc-data-layers](../../artifacts/misc-data-layers) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[pan-db](../../artifacts/pan-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[genomes-storage-db](../../artifacts/genomes-storage-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[functions](../../artifacts/functions) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[functional-enrichment-txt](../../artifacts/functional-enrichment-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Usage


This program computes functional enrichment within a pangenome and returns a <span class="artifact-n">[functional-enrichment-txt](/software/anvio/help/7.1/artifacts/functional-enrichment-txt)</span> file.

{:.warning}
For its sister programs, see <span class="artifact-n">[anvi-compute-metabolic-enrichment](/software/anvio/help/7.1/programs/anvi-compute-metabolic-enrichment)</span> and <span class="artifact-n">[anvi-compute-functional-enrichment-across-genomes](/software/anvio/help/7.1/programs/anvi-compute-functional-enrichment-across-genomes)</span>.

{:.notice}
Please also see <span class="artifact-n">[anvi-display-functions](/software/anvio/help/7.1/programs/anvi-display-functions)</span> which can both calculate functional enrichment, AND give you an interactive interface to display the distribution of functions.

## Enriched functions in a pangenome

For this to run, you must provide a <span class="artifact-n">[pan-db](/software/anvio/help/7.1/artifacts/pan-db)</span> and <span class="artifact-n">[genomes-storage-db](/software/anvio/help/7.1/artifacts/genomes-storage-db)</span> pair, as well as a <span class="artifact-n">[misc-data-layers](/software/anvio/help/7.1/artifacts/misc-data-layers)</span> that associates genomes in your pan database with categorical data. The program will then find functions that are enriched in each group (i.e., functions that are associated with gene clusters that are characteristic of the genomes in that group). 

{:.notice}
Note that your <span class="artifact-n">[genomes-storage-db](/software/anvio/help/7.1/artifacts/genomes-storage-db)</span> must have at least one functional annotation source for this to work.

This analysis will help you identify functions that are associated with a specific group of genomes in a pangenome and determine the functional core of your pangenome. For example, in the *Prochlorococcus* pangenome (the one used in [the pangenomics tutorial, where you can find more info about this program](http://merenlab.org/2016/11/08/pangenomics-v2/#making-sense-of-functions-in-your-pangenome)), this program finds that `Exonuclease VII` is enriched in the `low-light` genomes and not in `high-light` genomes. The output file provides various statistics about how confident the program is in making this association.

### How does it work?

What this program does can be broken down into three steps:

1. **Determine groups of genomes**. The program uses a <span class="artifact-n">[misc-data-layers](/software/anvio/help/7.1/artifacts/misc-data-layers)</span> variable (containing categorical, not numerical, data) to split genomes in a pangenome into two or more groups. For example, in the pangenome tutorial, the categorical variable name was `light` that partitioned genomes into `low-light` and `high-light `groups.

2.  **Determine the "functional associations" of gene clusters**. In short, this is collecting the functional annotations for all of the genes in each cluster and assigning the one that appears most frequently to represent the entire cluster.

3. **Quantify the distribution of functions in each group of genomes**. For this, the program determines to what extent a particular function is enriched in specific groups of genomes and reports it as a <span class="artifact-n">[functional-enrichment-txt](/software/anvio/help/7.1/artifacts/functional-enrichment-txt)</span> file. It does so by running the script `anvi-script-enrichment-stats`. 

{:.notice}
The script `anvi-script-enrichment-stats` was implemented by [Amy Willis](https://github.com/adw96), and described first in [this paper](https://doi.org/10.1186/s13059-020-02195-w).

{:.notice}
Check out [Alon's behind the scenes post](http://merenlab.org/2016/11/08/pangenomics-v2/#making-sense-of-functions-in-your-pangenome), which goes into a lot more detail.

### Basic usage

Here is the simplest way to run this program:

<div class="codeblock" markdown="1">
anvi&#45;compute&#45;functional&#45;enrichment&#45;in&#45;pan &#45;p <span class="artifact&#45;n">[pan&#45;db](/software/anvio/help/7.1/artifacts/pan&#45;db)</span>\
                                          &#45;g <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/7.1/artifacts/genomes&#45;storage&#45;db)</span> \
                                          &#45;o <span class="artifact&#45;n">[functional&#45;enrichment&#45;txt](/software/anvio/help/7.1/artifacts/functional&#45;enrichment&#45;txt)</span> \
                                          &#45;&#45;category&#45;variable CATEGORY \
                                          &#45;&#45;annotation&#45;source FUNCTION_SOURCE
</div>

The <span class="artifact-n">[pan-db](/software/anvio/help/7.1/artifacts/pan-db)</span> must contain at least one categorical data layer in <span class="artifact-n">[misc-data-layers](/software/anvio/help/7.1/artifacts/misc-data-layers)</span>, and you must choose one of these categories to define your pan-groups with the `--category-variable` parameter. You can see available variables with <span class="artifact-n">[anvi-show-misc-data](/software/anvio/help/7.1/programs/anvi-show-misc-data)</span> program with the parameters `-t layers --debug`.

Note that by default any genomes not in a category will be ignored; you can instead include these in the analysis by using the flag `--include-ungrouped`.

The <span class="artifact-n">[genomes-storage-db](/software/anvio/help/7.1/artifacts/genomes-storage-db)</span> must have at least one functional annotation source, and you must choose one of these sources with the `--annotation-source`. If you do not know which functional annotation sources are available in your <span class="artifact-n">[genomes-storage-db](/software/anvio/help/7.1/artifacts/genomes-storage-db)</span>, you can use the `--list-annotation-sources` parameter to find out.

### Additional options

By default, gene clusters with the same functional annotation will be merged. But if you provide the `--include-gc-identity-as-function` parameter and set the annotation source to be 'IDENTITY', anvi'o will treat gene cluster names as functions and enable you to investigate enrichment of each gene cluster independently. This is how you do it:

<div class="codeblock" markdown="1">
anvi&#45;compute&#45;functional&#45;enrichment&#45;in&#45;pan &#45;p <span class="artifact&#45;n">[pan&#45;db](/software/anvio/help/7.1/artifacts/pan&#45;db)</span>\
                                          &#45;g <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/7.1/artifacts/genomes&#45;storage&#45;db)</span> \
                                          &#45;o <span class="artifact&#45;n">[functional&#45;enrichment&#45;txt](/software/anvio/help/7.1/artifacts/functional&#45;enrichment&#45;txt)</span> \
                                          &#45;&#45;category&#45;variable CATEGORY \
                                          &#45;&#45;annotation&#45;source IDENTITY \
                                          &#45;&#45;include&#45;gc&#45;identity&#45;as&#45;function
</div>

To output a functional occurrence table, which describes the number of times each of your functional associations occurs in each genome you're looking at, use the `--functional-occurrence-table-output` parameter, like so:

<div class="codeblock" markdown="1">
anvi&#45;compute&#45;functional&#45;enrichment&#45;in&#45;pan &#45;p <span class="artifact&#45;n">[pan&#45;db](/software/anvio/help/7.1/artifacts/pan&#45;db)</span>\
                                          &#45;g <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/7.1/artifacts/genomes&#45;storage&#45;db)</span> \
                                          &#45;o <span class="artifact&#45;n">[functional&#45;enrichment&#45;txt](/software/anvio/help/7.1/artifacts/functional&#45;enrichment&#45;txt)</span> \
                                          &#45;&#45;category&#45;variable CATEGORY \
                                          &#45;&#45;annotation&#45;source FUNCTION_SOURCE \
                                          &#45;&#45;functional&#45;occurrence&#45;table&#45;output FUNC_OCCURRENCE.TXT
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-compute-functional-enrichment-in-pan.md) to update this information.


## Additional Resources


* [A description of the enrichment script run by this program can be found in Shaiber et al 2020](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02195-w)

* [An example of pangenome functional enrichment in the context of the Prochlorococcus metapangenome from Delmont and Eren 2018 is included in the pangenomics tutorial](http://merenlab.org/2016/11/08/pangenomics-v2/)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-compute-functional-enrichment-in-pan) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

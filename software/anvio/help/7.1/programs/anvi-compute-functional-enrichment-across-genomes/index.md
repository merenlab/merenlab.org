---
layout: page
title: anvi-compute-functional-enrichment-across-genomes [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-compute-functional-enrichment-across-genomes
image:
  featurerelative: ../../../images/header.png
  display: true
---

A program that computes functional enrichment across groups of genomes..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/ivagljiva.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Iva Veseli</span><div class="page-author-social-box"><a href="mailto:iveseli@uchicago.edu" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://github.com/ivagljiva" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/adw96.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Amy D. Willis</span><div class="page-author-social-box"><a href="http://statisticaldiversitylab.com/team/amy-willis" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:adwillis@uw.edu" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/AmyDWillis" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/adw96" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[groups-txt](../../artifacts/groups-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[genomes-storage-db](../../artifacts/genomes-storage-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[external-genomes](../../artifacts/external-genomes) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[internal-genomes](../../artifacts/internal-genomes) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[functions](../../artifacts/functions) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[functional-enrichment-txt](../../artifacts/functional-enrichment-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Usage


This program computes functional enrichment across groups of genomes and returns a <span class="artifact-n">[functional-enrichment-txt](/software/anvio/help/7.1/artifacts/functional-enrichment-txt)</span> file.

{:.warning}
For its sister programs, see <span class="artifact-n">[anvi-compute-functional-enrichment-in-pan](/software/anvio/help/7.1/programs/anvi-compute-functional-enrichment-in-pan)</span> and <span class="artifact-n">[anvi-compute-metabolic-enrichment](/software/anvio/help/7.1/programs/anvi-compute-metabolic-enrichment)</span>.

{:.notice}
Please also see <span class="artifact-n">[anvi-display-functions](/software/anvio/help/7.1/programs/anvi-display-functions)</span> which can both calculate functional enrichment, AND give you an interactive interface to display the distribution of functions.

## Functional enrichment

You can use this program by combining genomes described through <span class="artifact-n">[external-genomes](/software/anvio/help/7.1/artifacts/external-genomes)</span>, <span class="artifact-n">[internal-genomes](/software/anvio/help/7.1/artifacts/internal-genomes)</span>, and/or stored in a <span class="artifact-n">[genomes-storage-db](/software/anvio/help/7.1/artifacts/genomes-storage-db)</span>. In addition to sources for your genomes, you will need to provide a <span class="artifact-n">[groups-txt](/software/anvio/help/7.1/artifacts/groups-txt)</span> file to declare which genome belongs to which group for enrichment analysis to consider.

### How does it work?

1. **Aggregate functions from all sources**. Gene calls in each genome are tallied according to their functional annotations from the given annotation source.

2. **Quantify the distribution of functions in each group of genomes**. This information is then used by `anvi-script-enrichment-stats` to fit a GLM to determine (1) the level that a particular functional annotation is unique to a single group and (2) the percent of genomes it appears in in each group. This produces a <span class="artifact-n">[functional-enrichment-txt](/software/anvio/help/7.1/artifacts/functional-enrichment-txt)</span> file.

{:.notice}
The script `anvi-script-enrichment-stats` was implemented by [Amy Willis](https://github.com/adw96), and described first in [this paper](https://doi.org/10.1186/s13059-020-02195-w).


### Basic usage

You can use it with a single source of genomes:

<div class="codeblock" markdown="1">
anvi&#45;compute&#45;functional&#45;enrichment&#45;across&#45;genomes &#45;i <span class="artifact&#45;n">[internal&#45;genomes](/software/anvio/help/7.1/artifacts/internal&#45;genomes)</span> \
                                                  &#45;o <span class="artifact&#45;n">[functional&#45;enrichment&#45;txt](/software/anvio/help/7.1/artifacts/functional&#45;enrichment&#45;txt)</span> \
                                                  &#45;&#45;annotation&#45;source FUNCTION_SOURCE
</div>

or many:

<div class="codeblock" markdown="1">
anvi&#45;compute&#45;functional&#45;enrichment&#45;across&#45;genomes &#45;i <span class="artifact&#45;n">[internal&#45;genomes](/software/anvio/help/7.1/artifacts/internal&#45;genomes)</span>\
                                                  &#45;e <span class="artifact&#45;n">[external&#45;genomes](/software/anvio/help/7.1/artifacts/external&#45;genomes)</span> \
                                                  &#45;g <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/7.1/artifacts/genomes&#45;storage&#45;db)</span> \
                                                  &#45;o <span class="artifact&#45;n">[functional&#45;enrichment&#45;txt](/software/anvio/help/7.1/artifacts/functional&#45;enrichment&#45;txt)</span> \
                                                  &#45;&#45;annotation&#45;source FUNCTION_SOURCE
</div>

### Additional Parameters

You can get a tab-delimited matrix describing the occurrence (counts) of each function within each genome using the `--functional-occurrence-table-output` parameter:

<div class="codeblock" markdown="1">
anvi&#45;compute&#45;functional&#45;enrichment&#45;across&#45;genomes &#45;i <span class="artifact&#45;n">[internal&#45;genomes](/software/anvio/help/7.1/artifacts/internal&#45;genomes)</span> \
                                                  &#45;o <span class="artifact&#45;n">[functional&#45;enrichment&#45;txt](/software/anvio/help/7.1/artifacts/functional&#45;enrichment&#45;txt)</span> \
                                                  &#45;&#45;annotation&#45;source FUNCTION_SOURCE
                                                  &#45;&#45;functional&#45;occurrence&#45;table&#45;output FUNC_OCCURRENCE.TXT
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-compute-functional-enrichment-across-genomes.md) to update this information.


## Additional Resources


* [A description of the enrichment script run by this program can be found in Shaiber et al 2020](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02195-w)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-compute-functional-enrichment-across-genomes) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

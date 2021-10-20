---
layout: page
title: anvi-run-scg-taxonomy [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-run-scg-taxonomy
image:
  featurerelative: ../../../images/header.png
  display: true
---

The purpose of this program is to affiliate single-copy core genes in an anvi&#x27;o contigs database with taxonomic names. A properly setup local SCG taxonomy database is required for this program to perform properly. After its successful run, `anvi-estimate-scg-taxonomy` will be useful to estimate taxonomy at genome-, collection-, or metagenome-level).

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


<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[scgs-taxonomy-db](../../artifacts/scgs-taxonomy-db) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[hmm-hits](../../artifacts/hmm-hits) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[scgs-taxonomy](../../artifacts/scgs-taxonomy) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>


## Usage


This program **associates the single-copy core genes in your <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> with taxnomy information.**  

Once this information is stored in your <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> (in the form of a <span class="artifact-n">[scgs-taxonomy](/software/anvio/help/7.1/artifacts/scgs-taxonomy)</span> artifact), you can run <span class="artifact-n">[anvi-estimate-scg-taxonomy](/software/anvio/help/7.1/programs/anvi-estimate-scg-taxonomy)</span> or use the <span class="artifact-n">[anvi-interactive](/software/anvio/help/7.1/programs/anvi-interactive)</span> and enable "Realtime taxonomy estimate for bins." Check out [this tutorial](http://merenlab.org/2019/10/08/anvio-scg-taxonomy/) for more information. 

In order to run this program, you'll need a <span class="artifact-n">[scgs-taxonomy-db](/software/anvio/help/7.1/artifacts/scgs-taxonomy-db)</span>, which you can set up by running <span class="artifact-n">[anvi-setup-scg-taxonomy](/software/anvio/help/7.1/programs/anvi-setup-scg-taxonomy)</span>. 

### What does this program do? 

In short, this program searches all of the single-copy core genes that it uses for this workflow (which are the 22 listed on [this page](https://github.com/merenlab/anvio/tree/master/anvio/data/misc/SCG_TAXONOMY/GTDB/SCG_SEARCH_DATABASES)) against the [GTDB](https://gtdb.ecogenomic.org/) databases that you downloaded, and stores hits in your <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>. In other words, it finds your single-copy core genes and assigns them taxonomy. This way, it can use these single-copy core genes later to estimate the taxnomy of larger groups of contigs that include these single-copy core genes when you run <span class="artifact-n">[anvi-estimate-scg-taxonomy](/software/anvio/help/7.1/programs/anvi-estimate-scg-taxonomy)</span>. 

### Sweet. How do I run it? 

<div class="codeblock" markdown="1">
anvi&#45;run&#45;scg&#45;taxonomy &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span>
</div>

In case you're running this on a genome and not getting any hits, you have the option to try lowering the percent identity required for a hit (as long as you're careful with it). The default value is 90 percent. 

<div class="codeblock" markdown="1">
anvi&#45;run&#45;scg&#45;taxonomy &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                      &#45;&#45;min&#45;percent&#45;identity 70
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-run-scg-taxonomy.md) to update this information.


## Additional Resources


* [Usage examples and warnings](http://merenlab.org/scg-taxonomy)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-run-scg-taxonomy) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

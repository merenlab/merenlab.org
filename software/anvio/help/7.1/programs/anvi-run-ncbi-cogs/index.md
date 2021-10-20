---
layout: page
title: anvi-run-ncbi-cogs [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-run-ncbi-cogs
image:
  featurerelative: ../../../images/header.png
  display: true
---

This program runs NCBI&#x27;s COGs to associate genes in an anvi&#x27;o contigs database with functions. COGs database was been designed as an attempt to classify proteins from completely sequenced genomes on the basis of the orthology concept..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[cogs-data](../../artifacts/cogs-data) <img src="../../images/icons/DATA.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[functions](../../artifacts/functions) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>


## Usage


This program **associates genes in your <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> with functions using NCBI's [Clusters of Orthologus Groups (COGs) database](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102395/).**

Before you run this program, you'll have to set up the COGs database on your computer with the program <span class="artifact-n">[anvi-setup-ncbi-cogs](/software/anvio/help/7.1/programs/anvi-setup-ncbi-cogs)</span>.

To run, you'll need to provide a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>. If you stored the <span class="artifact-n">[cogs-data](/software/anvio/help/7.1/artifacts/cogs-data)</span> that you got from running <span class="artifact-n">[anvi-setup-ncbi-cogs](/software/anvio/help/7.1/programs/anvi-setup-ncbi-cogs)</span> in a custom location, you'll need to provide that path as well. The output is a <span class="artifact-n">[functions](/software/anvio/help/7.1/artifacts/functions)</span> artifact.

<div class="codeblock" markdown="1">
anvi&#45;run&#45;ncbi&#45;cogs &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
            &#45;&#45;cog&#45;data&#45;dir path/to/<span class="artifact&#45;n">[cogs&#45;data](/software/anvio/help/7.1/artifacts/cogs&#45;data)</span>
</div>

Without the flag `--cog-data-dir`, anvi'o will just search in the default location.

By default, this program uses DIAMOND in the "fast" setting for database searching. To instead run in "sensitive" mode, just call:

<div class="codeblock" markdown="1">
anvi&#45;run&#45;ncbi&#45;cogs &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
            &#45;&#45;sensitive
</div>

You can also use blastp to search, by running:

<div class="codeblock" markdown="1">
anvi&#45;run&#45;ncbi&#45;cogs &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
            &#45;&#45;search&#45;with blastp
</div>




{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-run-ncbi-cogs.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-run-ncbi-cogs) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

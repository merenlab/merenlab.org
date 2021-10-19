---
layout: page
title: anvi-script-gen_stats_for_single_copy_genes.py [program]
categories: [anvio]
comments: false
redirect_from: /m/anvi-script-gen_stats_for_single_copy_genes.py
image:
  featurerelative: ../../../images/header.png
  display: true
---

A simple script to generate info from search tables, given a contigs-db.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[genes-stats](../../artifacts/genes-stats) <img src="../../images/icons/STATS.png" class="artifact-icon-mini" /></span></p>


## Usage


This program **provides information about each of the single copy core genes in your <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span>**. 

Simply provide a <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span>, and it will create a <span class="artifact-n">[genes-stats](/software/anvio/help/main/artifacts/genes-stats)</span> file containing a variety of information about the single copy core genes in your database. 

<div class="codeblock" markdown="1">
anvi&#45;script&#45;gen_stats_for_single_copy_core_genes.py &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/main/artifacts/contigs&#45;db)</span> 
</div>

The console output will tell you the total number of contigs, splits, and nucleotides in your <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span>, while the text output will tell you the source, name, and e-value of each single-copy core gene. 

You can get information from only single-copy core genes from a specific source. To see what sources are availible in your <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span>, run 

<div class="codeblock" markdown="1">
anvi&#45;script&#45;gen_stats_for_single_copy_core_genes.py &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/main/artifacts/contigs&#45;db)</span> \
                                                    &#45;&#45;list&#45;sources
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-script-gen_stats_for_single_copy_genes.py.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-gen_stats_for_single_copy_genes.py) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

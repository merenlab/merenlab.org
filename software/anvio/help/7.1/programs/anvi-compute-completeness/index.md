---
layout: page
title: anvi-compute-completeness [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-compute-completeness
image:
  featurerelative: ../../../images/header.png
  display: true
---

A script to generate completeness info for a given list of _splits_.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[splits-txt](../../artifacts/splits-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[hmm-source](../../artifacts/hmm-source) <img src="../../images/icons/HMM.png" class="artifact-icon-mini" /></span></p>


## Can provide


This program does not seem to provide any artifacts. Such programs usually print out some information for you to see or alter some anvi'o artifacts without producing any immediate outputs.


## Usage


This program tells you the completeness and redundency of single-copy gene sources available for your <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>. 

For example, some of the defaults are collections of single-copy core genes named  `Protista_83`, `Archaea_76`, and `Bacteria_71`. This program will give you a rough estimate of how many Protist, Archaeal, and Bacterial genomes are included in your dataset using these single-copy core genes. 

You can use the following run to list available completeness sources in your <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>:

<div class="codeblock" markdown="1">
anvi&#45;compute&#45;completeness &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                          &#45;&#45;list&#45;completeness&#45;sources
</div>
                              
Then you can run this program on a specifc source as folows:

<div class="codeblock" markdown="1">
anvi&#45;compute&#45;completeness &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                          &#45;&#45;completeness&#45;source Bacteria_71
</div>
                              
You can also provide a <span class="artifact-n">[splits-txt](/software/anvio/help/7.1/artifacts/splits-txt)</span> to focus on a specific set of splits, or declare a minimum e-value for a gene to count as a hit. The default is `1e-15`.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-compute-completeness.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-compute-completeness) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

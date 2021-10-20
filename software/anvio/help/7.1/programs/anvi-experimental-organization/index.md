---
layout: page
title: anvi-experimental-organization [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-experimental-organization
image:
  featurerelative: ../../../images/header.png
  display: true
---

Create an experimental clustering dendrogram..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[clustering-configuration](../../artifacts/clustering-configuration) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[dendrogram](../../artifacts/dendrogram) <img src="../../images/icons/NEWICK.png" class="artifact-icon-mini" /></span></p>


## Usage


This program can use an anvi'o <span class="artifact-n">[clustering-configuration](/software/anvio/help/7.1/artifacts/clustering-configuration)</span> file to access various data sources in anvi'o databases to produce a hierarchical clustering dendrogram for items.

It is especially powerful when the user wishes to create a hierarchical clustering of contigs or gene clusters using only a specific set of samples. If you would like to see an example usage of this program see the article on [combining metagenomics with metatranscriptomics](https://merenlab.org/2015/06/10/combining-omics-data/).


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-experimental-organization.md) to update this information.


## Additional Resources


* [An example use of this program](https://merenlab.org/2015/06/10/combining-omics-data/)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-experimental-organization) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

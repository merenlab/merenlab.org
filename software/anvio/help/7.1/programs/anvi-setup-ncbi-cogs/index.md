---
layout: page
title: anvi-setup-ncbi-cogs [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-setup-ncbi-cogs
image:
  featurerelative: ../../../images/header.png
  display: true
---

Download and setup NCBI&#x27;s Clusters of Orthologous Groups database.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


This program seems to know what its doing. It needs no input material from its user. Good program.


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[cogs-data](../../artifacts/cogs-data) <img src="../../images/icons/DATA.png" class="artifact-icon-mini" /></span></p>


## Usage


This program **downloads and organizes a local copy of the data from NCBI's [COGs database](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102395/) for use in function annotation.** This program generates a <span class="artifact-n">[cogs-data](/software/anvio/help/7.1/artifacts/cogs-data)</span> artifact, which is required to run the program <span class="artifact-n">[anvi-run-ncbi-cogs](/software/anvio/help/7.1/programs/anvi-run-ncbi-cogs)</span>. 

### Set up COGs data
<div class="codeblock" markdown="1">
anvi&#45;setup&#45;ncbi&#45;cogs &#45;&#45;just&#45;do&#45;it
</div>

If you already have a <span class="artifact-n">[cogs-data](/software/anvio/help/7.1/artifacts/cogs-data)</span> artifact and are trying to redownload this data, run 

<div class="codeblock" markdown="1">
anvi&#45;setup&#45;ncbi&#45;cogs &#45;&#45;reset
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-setup-ncbi-cogs.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-setup-ncbi-cogs) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

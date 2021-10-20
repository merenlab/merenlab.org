---
layout: page
title: anvi-setup-scg-taxonomy [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-setup-scg-taxonomy
image:
  featurerelative: ../../../images/header.png
  display: true
---

The purpose of this program is to download necessary information from GTDB (https://gtdb.ecogenomic.org/), and set it up in such a way that your anvi&#x27;o installation is able to assign taxonomy to single-copy core genes using `anvi-run-scg-taxonomy` and estimate taxonomy for genomes or metagenomes using `anvi-estimate-scg-taxonomy`).

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


This program seems to know what its doing. It needs no input material from its user. Good program.


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[scgs-taxonomy-db](../../artifacts/scgs-taxonomy-db) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>


## Usage


This program **downloads and sets up the search databases used for the scg-taxonomy workflow** (from [GTDB](https://gtdb.ecogenomic.org/)) so that you can run <span class="artifact-n">[anvi-run-scg-taxonomy](/software/anvio/help/7.1/programs/anvi-run-scg-taxonomy)</span> and <span class="artifact-n">[anvi-estimate-scg-taxonomy](/software/anvio/help/7.1/programs/anvi-estimate-scg-taxonomy)</span>. This program generates a <span class="artifact-n">[scgs-taxonomy-db](/software/anvio/help/7.1/artifacts/scgs-taxonomy-db)</span> artifact, which is required to run both of those programs. 

For more information on that workflow, check out [this page](http://merenlab.org/2019/10/08/anvio-scg-taxonomy/)

You will only have to run this program once per anvi'o installation. 

Why is this not done by default? It just makes things easier downstream to build these databases with the DIAMOND installed on your computer to avoid incompatibility issues. Besides, it should take under a minute and is as simple as running

<div class="codeblock" markdown="1">
anvi&#45;setup&#45;scg&#45;taxonomy
</div>

If you have already already run this program and are trying to redownload this data, run 

<div class="codeblock" markdown="1">
anvi&#45;setup&#45;scg&#45;taxonomy &#45;&#45;reset
</div>

You can also download a specific release of this database by providing its URL with the flag `--scg-taxonomy-remote-database-url`. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-setup-scg-taxonomy.md) to update this information.


## Additional Resources


* [Usage examples and warnings](http://merenlab.org/scg-taxonomy)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-setup-scg-taxonomy) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

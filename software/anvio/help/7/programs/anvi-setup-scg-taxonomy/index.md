---
layout: page
title: anvi-setup-scg-taxonomy [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

The purpose of this program is to download necessary information from GTDB (https://gtdb.ecogenomic.org/), and set it up in such a way that your anvi&#x27;o installation is able to assign taxonomy to single-copy core genes using `anvi-run-scg-taxonomy` and estimate taxonomy for genomes or metagenomes using `anvi-estimate-scg-taxonomy`).

See **[program help menu](../../../../vignette#anvi-setup-scg-taxonomy)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[scgs-taxonomy-db](../../artifacts/scgs-taxonomy-db) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"></p>

## Usage


This program **downloads and sets up the search databases used for the scg-taxonomy workflow** (from [GTDB](https://gtdb.ecogenomic.org/)) so that you can run <span class="artifact-n">[anvi-run-scg-taxonomy](/software/anvio/help/7/programs/anvi-run-scg-taxonomy)</span> and <span class="artifact-n">[anvi-estimate-scg-taxonomy](/software/anvio/help/7/programs/anvi-estimate-scg-taxonomy)</span>. This program generates a <span class="artifact-n">[scgs-taxonomy-db](/software/anvio/help/7/artifacts/scgs-taxonomy-db)</span> artifact, which is required to run both of those programs. 

For more information on that workflow, check out [this page](http://merenlab.org/2019/10/08/anvio-scg-taxonomy/)

You will only have to run this program once per anvi'o installation. 

Why is this not done by default? It just makes things easier downstream to build these databases with the DIAMOND installed on your computer to avoid incompatibility issues. Besides, it should take under a minute and is as simple as running

<div class="codeblock" markdown="1">
anvi&#45;setup&#45;scg&#45;databases
</div>

If you have already already run this program and are trying to redownload this data, run 

<div class="codeblock" markdown="1">
anvi&#45;setup&#45;scg&#45;databases &#45;&#45;reset
</div>

You can also download a specific release of this database by providing its URL with the flag `--scg-taxonomy-remote-database-url`. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-setup-scg-taxonomy.md) to update this information.


## Additional Resources


* [Usage examples and warnings](http://merenlab.org/scg-taxonomy)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-setup-scg-taxonomy) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

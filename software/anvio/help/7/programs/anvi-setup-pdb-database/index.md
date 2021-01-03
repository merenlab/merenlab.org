---
layout: page
title: anvi-setup-pdb-database [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Setup or update an offline database of representative PDB structures clustered at 95%.

See **[program help menu](../../../../vignette#anvi-setup-pdb-database)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[pdb-db](../../artifacts/pdb-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"></p>

## Usage



## Basic usage 

This program creates a <span class="artifact-n">[pdb-db](/software/anvio/help/7/artifacts/pdb-db)</span> local database that holds PDB structures from [this sequence database](https://salilab.org/modeller/supplemental.html), which is hosted by the [Sali lab](https://salilab.org/).  Their database comprises all PDB RCSB sequences that have been clustered at 95% sequence similarity. They seem to update their database every couple of months (thank you guys!).


The purpose of <span class="artifact-n">[anvi-setup-pdb-database](/software/anvio/help/7/programs/anvi-setup-pdb-database)</span> to have a local copy of reference structures that can be used to, for example, get template structures for homology modelling when <span class="artifact-n">[anvi-gen-structure-database](/software/anvio/help/7/programs/anvi-gen-structure-database)</span> is ran.


Running this program is easy:

<div class="codeblock" markdown="1">
anvi&#45;setup&#45;pdb&#45;database &#45;&#45;just&#45;do&#45;it
</div>

If you already have a <span class="artifact-n">[pdb-db](/software/anvio/help/7/artifacts/pdb-db)</span> artifact and are trying to redownload this data, run 

<div class="codeblock" markdown="1">
anvi&#45;setup&#45;pdb&#45;database &#45;&#45;reset
</div>

Or if you just want to update your database, run 

<div class="codeblock" markdown="1">
anvi&#45;setup&#45;pdb&#45;database &#45;&#45;update
</div>

## Notes

The output <span class="artifact-n">[pdb-db](/software/anvio/help/7/artifacts/pdb-db)</span> database is ~20GB and its contents may take several hours to download.



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-setup-pdb-database.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-setup-pdb-database) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

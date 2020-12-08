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

See **[program help menu](../../../vignette#anvi-setup-pdb-database)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[pdb-db](../../artifacts/pdb-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>

## Requires or uses

<p style="text-align: left" markdown="1"></p>

## Usage


This program downloads a section of the [Protein Data Bank](https://www.rcsb.org/) which is required for strucutral analyses in anvi'o, such as running <span class="artifact-n">[anvi-gen-structure-database](/software/anvio/help/main/programs/anvi-gen-structure-database)</span> or <span class="artifact-n">[anvi-3dev](/software/anvio/help/main/programs/anvi-3dev)</span>. 

<div class="codeblock" markdown="1">
anvi&#45;setup&#45;pdb&#45;database &#45;&#45;just&#45;do&#45;it
</div>

If you already have a <span class="artifact-n">[pdb-db](/software/anvio/help/main/artifacts/pdb-db)</span> artifact and are trying to redownload this data, run 

<div class="codeblock" markdown="1">
anvi&#45;setup&#45;pdb&#45;database &#45;&#45;reset
</div>

Or if you just want to update your database, run 

<div class="codeblock" markdown="1">
anvi&#45;setup&#45;pdb&#45;database &#45;&#45;update
</div>



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-setup-pdb-database.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-setup-pdb-database) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

---
layout: page
title: anvi-setup-trna-taxonomy [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

The purpose of this program is to setup necessary databases for tRNA genes collected from GTDB (https://gtdb.ecogenomic.org/), genomes in your local anvi&#x27;o installation so taxonomy information for a given set of tRNA sequences can be identified using `anvi-run-trna-taxonomy` and made sense of via `anvi-estimate-trna-taxonomy`).

See **[program help menu](../../../../vignette#anvi-setup-trna-taxonomy)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[trna-taxonomy-db](../../artifacts/trna-taxonomy-db) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"></p>

## Usage


This program downloads a local copy of a subset of the databases from [GTDB](https://gtdb.ecogenomic.org/) (stored in a <span class="artifact-n">[trna-taxonomy-db](/software/anvio/help/7/artifacts/trna-taxonomy-db)</span>), so that tRNA sequences in your dataset can be associated with taxonomy information. It is required to run this program before you can run <span class="artifact-n">[anvi-run-trna-taxonomy](/software/anvio/help/7/programs/anvi-run-trna-taxonomy)</span> or <span class="artifact-n">[anvi-estimate-trna-taxonomy](/software/anvio/help/7/programs/anvi-estimate-trna-taxonomy)</span>.

Like other `anvi-setup-` programs, this only needs to be run once per anvi'o version. The default path is `anvio/data/misc/TRNA-TAXONOMY`. You can store the resulting <span class="artifact-n">[trna-taxonomy-db](/software/anvio/help/7/artifacts/trna-taxonomy-db)</span> in a custom location if desired), but then you'll need to provide the path to it whenever you run <span class="artifact-n">[anvi-run-trna-taxonomy](/software/anvio/help/7/programs/anvi-run-trna-taxonomy)</span>. 

To run this program, you can simply run

<div class="codeblock" markdown="1">
anvi&#45;setup&#45;trna&#45;taxonomy 
</div>

If you are trying to redownload these databases, run: 

<div class="codeblock" markdown="1">
anvi&#45;setup&#45;trna&#45;taxonomy &#45;&#45;reset
</div>

Alternatively, you can use `--redo-databases` if you just want to update the database version without redownloading the data. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-setup-trna-taxonomy.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-setup-trna-taxonomy) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

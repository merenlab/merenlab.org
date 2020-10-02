---
layout: page
title: anvi-setup-trna-taxonomy [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

The purpose of this program is to setup necessary databases for tRNA genes collected from GTDB (https://gtdb.ecogenomic.org/), genomes in your local anvi&#39;o installation so taxonomy information for a given set of tRNA sequences can be identified using `anvi-run-trna-taxonomy` and made sense of via `anvi-estimate-trna-taxonomy`)..

See **[program help menu](../../../vignette#anvi-setup-trna-taxonomy)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[trna-taxonomy-db](../../artifacts/trna-taxonomy-db)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"></p>

## Usage


{:.notice}
**No one has described the usage of this program** :/ If you would like to contribute, please see previous examples [here](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs), and feel free to add a Markdown formatted file in that directory named "anvi-setup-trna-taxonomy.md". For a template, you can use the markdown file for `anvi-gen-contigs-database`. THANK YOU!


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-setup-trna-taxonomy) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

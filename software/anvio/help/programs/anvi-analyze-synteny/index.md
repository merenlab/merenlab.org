---
layout: page
title: anvi-analyze-synteny [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

This program helps you extract synteny patterns from a group of similar loci or genomes. This is accomplished by deconstructing contigs into Ngrams (group of neighboring genes of which N is the number of genes) using a sliding window method. Once completed, anvi&#39;o will export a table with Ngrams for you to work with. If more than one annotation source is provided, please designate which annotation_source anvio should use to make the Ngrams with using the --ngram-source parameter (e.g. gene_clusters, COG_FUNCTION). By default anvi&#39;o will ignore Ngrams that contain genes without annotations. If you would like to override this, you can use the --analyze-unknown-functions flag..

See **[program help menu](../../../vignette#anvi-analyze-synteny)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[ngrams](../../artifacts/ngrams)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span> <span class="artifact-r">[genomes-storage-db](../../artifacts/genomes-storage-db)</span></p>

## Usage


{:.notice}
**No one has described the usage of this program** :/ If you would like to contribute, please see previous examples [here](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs), and feel free to add a Markdown formatted file in that directory named "anvi-analyze-synteny.md". For a template, you can use the markdown file for `anvi-gen-contigs-database`. THANK YOU!


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-analyze-synteny) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

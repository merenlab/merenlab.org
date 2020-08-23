---
layout: page
title: anvi-script-calculate-pn-ps-ratio [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

This program calculates for each gene the ratio of pN/pS (the metagenomic analogy of dN/dS; see doi:10.1038/nature11711 and doi:10.7717/peerj.2959) based on metagenomic read recruitment, however, unlike standard pN/pS calculations, it relies on codons rather than nucleotides for accurate estimations of synonimity.

See **[program help menu](../../../vignette#anvi-script-calculate-pn-ps-ratio)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span> <span class="artifact-r">[variability-profile-txt](../../artifacts/variability-profile-txt)</span></p>

## Usage


{:.notice}
**No one has described the usage of this program** :/ If you would like to contribute, please see previous examples [here](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs), and feel free to add a Markdown formatted file in that directory named "anvi-script-calculate-pn-ps-ratio.md". For a template, you can use the markdown file for `anvi-gen-contigs-database`. THANK YOU!


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-calculate-pn-ps-ratio) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

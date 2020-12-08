---
layout: page
title: anvi-script-gen-distribution-of-genes-in-a-bin [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Quantify the detection of genes in genomes in metagenomes to identify the environmental core. This is a helper script for anvi&#x27;o metapangenomic workflow.

See **[program help menu](../../../vignette#anvi-script-gen-distribution-of-genes-in-a-bin)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span> <span class="artifact-r">[profile-db](../../artifacts/profile-db)</span> <span class="artifact-r">[collection](../../artifacts/collection)</span> <span class="artifact-r">[bin](../../artifacts/bin)</span></p>

## Usage


{:.notice}
**No one has described the usage of this program** :/ If you would like to contribute, please see previous examples [here](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs), and feel free to add a Markdown formatted file in that directory named "anvi-script-gen-distribution-of-genes-in-a-bin.md". For a template, you can use the markdown file for `anvi-gen-contigs-database`. THANK YOU!


## Additional Resources


* [This program in action as part of the metapangenomic workflow](http://merenlab.org/data/prochlorococcus-metapangenome/#classification-of-genes-as-ecgs-and-eags-by-the-distribution-of-genes-in-a-genome-across-metagenomes)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-gen-distribution-of-genes-in-a-bin) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

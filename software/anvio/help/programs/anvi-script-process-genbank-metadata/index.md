---
layout: page
title: anvi-script-process-genbank-metadata [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

This script takes the &#39;metadata&#39; output of the program `ncbi-genome-download` (see https://github.com/kblin/ncbi-genome-download for details), and processes each GenBank file found in the metadata file to generate a FASTA file, as well as genes and functions files for each entry. Plus, it autmatically generates a FASTA TXT file descriptor for anvi&#39;o snakemake workfloes. So it is a multi-talented program like that..

See **[program help menu](../../../vignette#anvi-script-process-genbank-metadata)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[contigs-fasta](../../artifacts/contigs-fasta)</span> <span class="artifact-p">[functions-txt](../../artifacts/functions-txt)</span> <span class="artifact-p">[external-gene-calls](../../artifacts/external-gene-calls)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"></p>

## Usage


{:.notice}
**No one has described the usage of this program** :/ If you would like to contribute, please see previous examples [here](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs), and feel free to add a Markdown formatted file in that directory named "anvi-script-process-genbank-metadata.md". For a template, you can use the markdown file for `anvi-gen-contigs-database`. THANK YOU!


## Additional Resources


* [A tutorial on using this program to access NCBI genomes for &#39;omics analyses in Anvi&#39;o](http://merenlab.org/2019/03/14/ncbi-genome-download-magic/)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-process-genbank-metadata) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

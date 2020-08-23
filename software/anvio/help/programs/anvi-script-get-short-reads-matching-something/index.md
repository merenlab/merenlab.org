---
layout: page
title: anvi-script-get-short-reads-matching-something [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

You give this program one or more FASTQ files and a short sequence, and it returns all short reads from the FASTQ file that matches to it. The purpose of this is to get back short reads that may be extending into hypervariable regions of genomes, resulting a decreased mappability of short reads in the metagenome given a reference. You often see those areas of genomes as significant dips in coverage, and in most cases with a large number of SNVs. When you provide the downstream conserved sequence, this program allows you to take a better look at those regions at the short read level without any mapping.

See **[program help menu](../../../vignette#anvi-script-get-short-reads-matching-something)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[short-reads-fasta](../../artifacts/short-reads-fasta)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[fasta](../../artifacts/fasta)</span></p>

## Usage


{:.notice}
**No one has described the usage of this program** :/ If you would like to contribute, please see previous examples [here](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs), and feel free to add a Markdown formatted file in that directory named "anvi-script-get-short-reads-matching-something.md". For a template, you can use the markdown file for `anvi-gen-contigs-database`. THANK YOU!


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-get-short-reads-matching-something) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

---
layout: page
title: anvi-script-compute-ani-for-fasta [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Run ANI between contigs in a single FASTA file.

See **[program help menu](../../../../vignette#anvi-script-compute-ani-for-fasta)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[genome-similarity](../../artifacts/genome-similarity) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[fasta](../../artifacts/fasta) <img src="../../images/icons/FASTA.png" class="artifact-icon-mini" /></span></p>

## Usage


This program computes the average nucleotide identity between reads in a single fasta file (using PyANI). 

To compute the ANI (or other genome distance metrics) between two genomes in different fasta files, use <span class="artifact-n">[anvi-compute-genome-similarity](/software/anvio/help/7/programs/anvi-compute-genome-similarity)</span>. 

A default run of this program looks like this: 

<div class="codeblock" markdown="1">
anvi&#45;script&#45;compute&#45;ani&#45;for&#45;fasta &#45;f <span class="artifact&#45;n">[fasta](/software/anvio/help/7/artifacts/fasta)</span> \ 
                                  &#45;o path/to/output \
                                  &#45;&#45;method ANIb
</div>

By default, the PyANI method is ANIb (which aligns 1020 nt fragments of your sequences using BLASTN+). You can switch to ANIm, ANIblastall, or TETRA if desired. See the [PyANI documentation](https://github.com/widdowquinn/pyani) for more informaiton. 

You also have the option to change the distance metric (from the default "euclidean") or the linkage method (from the default "ward") or provide a path to a log file for debug messages. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-script-compute-ani-for-fasta.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-compute-ani-for-fasta) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

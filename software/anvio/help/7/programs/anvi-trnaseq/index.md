---
layout: page
title: anvi-trnaseq [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

A program to process raw tRNA-seq dataset, which is the sequencing of tRNA transcripts in a given sample, to generate an anvi&#x27;o tRNA-seq database.

See **[program help menu](../../../../vignette#anvi-trnaseq)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[trnaseq-db](../../artifacts/trnaseq-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[trnaseq-fasta](../../artifacts/trnaseq-fasta) <img src="../../images/icons/FASTA.png" class="artifact-icon-mini" /></span></p>

## Usage


The input for this program is a **properly formatted** <span class="artifact-n">[trnaseq-fasta](/software/anvio/help/7/artifacts/trnaseq-fasta)</span>, containing sequences from a tRNA-seq sample or split.

The program identifies tRNA among the input sequences, profiles the tRNA primary sequence and secondary structure, and filters single nucleotide variants from modified nucleotides.

The primary output of the program is a <span class="artifact-n">[trnaseq-db](/software/anvio/help/7/artifacts/trnaseq-db)</span>. Supplemental files are also produced: an analysis summary file, a tab-separated file of unique sequences not identified as tRNA, and a tab-separated file showing the range of 5' and 3' variants trimmed from tRNA sequences.

We encourage you to read the list of options in the `anvi-trnaseq --help` menu to understand how the user can manipulate the multifaceted analyses performed by the program.

The program can generate a .ini file for tRNA feature parameterization using an alternate command, `anvi-trnaseq --default-feature-param-file <param.ini>`. The default parameterizations in the file can be modified by the user, and the file can be used as the `--feature-param-file` argument in the main mode of the program. `anvi-trnaseq --print-default-feature-params` can also be used to quickly and neatly display the defaults in the terminal.


### Create a tRNA-seq database from a FASTA file, using 16 cores

<div class="codeblock" markdown="1">
anvi&#45;trnaseq &#45;f <span class="artifact&#45;n">[trnaseq&#45;fasta](/software/anvio/help/7/artifacts/trnaseq&#45;fasta)</span> \
             &#45;S example_sample_name \
             &#45;o example_empty_output_directory_path \
             &#45;T 16
</div>

### Create a tRNA-seq database from a sample identified as a demethylase split, overwriting the output directory if it already exists

<div class="codeblock" markdown="1">
anvi&#45;trnaseq &#45;f <span class="artifact&#45;n">[trnaseq&#45;fasta](/software/anvio/help/7/artifacts/trnaseq&#45;fasta)</span> \
             &#45;S example_sample_name \
             &#45;o example_empty_output_directory_path \
             &#45;T 16 \
             &#45;&#45;treatment demethylase \
             &#45;W
</div>

{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-trnaseq.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-trnaseq) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

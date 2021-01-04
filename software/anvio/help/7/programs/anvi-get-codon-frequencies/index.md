---
layout: page
title: anvi-get-codon-frequencies [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Get amino acid or codon frequencies of genes in a contigs database.

See **[program help menu](../../../../vignette#anvi-get-codon-frequencies)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[codon-frequencies-txt](../../artifacts/codon-frequencies-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[aa-frequencies-txt](../../artifacts/aa-frequencies-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>

## Usage


This program **calculates the frequency of each codon or amino acid of every gene in your <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span>**. 

To run with all standard parameters, simply provide a <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> and path for the output file as follows: 

<div class="codeblock" markdown="1">
anvi&#45;get&#45;codon&#45;frequencies &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> \ 
                &#45;o name/of/output_file.txt 
</div>

The output of this is a <span class="artifact-n">[codon-frequencies-txt](/software/anvio/help/7/artifacts/codon-frequencies-txt)</span> that counts the number of times each codon appears in all of your genes.

If instead you want to calculate the data for the amino acids, run 

<div class="codeblock" markdown="1">
anvi&#45;get&#45;codon&#45;frequencies &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> \ 
                &#45;o name/of/output_file.txt  \
                &#45;&#45;return&#45;AA&#45;frequencies&#45;instead \
                &#45;&#45;gene&#45;caller&#45;id MY_FAVORITE_GENE
</div>

In this example, the flag `gene-caller-id` means that it will only count the amino acid frequencies of a single gene, namely `MY_FAVORITE_GENE`.

You can also return the data as a percent of the total number of codons or amino acids in the gene (with the flag `--percent-normalize`) or calculate the percent that each codon encoding the same amino acid appears in the gene (for example, 0.4 GCT and 0.6 GCC for alanine) (with the flag `--merens-codon-normalization`). 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-get-codon-frequencies.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-get-codon-frequencies) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

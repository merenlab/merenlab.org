---
layout: page 
title: anvi-gen-contigs-database [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/PROGRAM.png" alt="PROGRAM" style="width:100px; border:none" />

Generate a new anvi&#39;o contigs database.

[Back to help main page](../../)

## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[contigs-db](../../artifacts/contigs-db)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-fasta](../../artifacts/contigs-fasta)</span> <span class="artifact-r">[external-gene-calls](../../artifacts/external-gene-calls)</span></p>

## Usage


### Create a contigs database from a FASTA file

<div class="codeblock" markdown="1">
anvi-gen-contigs-database -f [contigs-fasta](/software/anvio/help/artifacts/contigs-fasta) \
                          -o [contigs-db](/software/anvio/help/artifacts/contigs-db)
</div>

A FASTA file that contains one or more sequences. These sequences may belong
to a single genome, or could be contigs obtained from an assembly.

### Create a contigs database with external gene calls

<div class="codeblock" markdown="1">
anvi-gen-contigs-database -f [contigs-fasta](/software/anvio/help/artifacts/contigs-fasta) \
                          -o [contigs-db](/software/anvio/help/artifacts/contigs-db) \
                          -e [external-gene-calls](/software/anvio/help/artifacts/external-gene-calls)
</div>




## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-gen-contigs-database) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

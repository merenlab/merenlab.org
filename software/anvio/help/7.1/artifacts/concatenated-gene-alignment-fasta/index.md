---
layout: page
title: concatenated-gene-alignment-fasta [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/concatenated-gene-alignment-fasta
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/FASTA.png" alt="FASTA" style="width:100px; border:none" />

A FASTA-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-get-sequences-for-gene-clusters](../../programs/anvi-get-sequences-for-gene-clusters)</span> <span class="artifact-p">[anvi-get-sequences-for-hmm-hits](../../programs/anvi-get-sequences-for-hmm-hits)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-gen-phylogenomic-tree](../../programs/anvi-gen-phylogenomic-tree)</span></p>


## Description

This file **contains the alignment information for multiple genes across different organisms**.

Basically, a single gene alignment compares a single gene's sequence across multiple organisms. For example, you could align some specific rRNA sequence across all of the organisms in your sample. This alignment highlights both mutations and insertions and deletions (indicated with dashes). 

Clustal programs do a great job of visualizing this data, by color coding it. Here is an example from Anvi'o's pangenome display: 

![A lovely clustal-like alignment from the anvi'o pangenome display](../../images/example_alignment.png)

A concatenated gene alignment fasta contains multiple of these gene alignments, in order to generate a tree based off of multiple genes. 

This information can then be used to generate a phylogenomic tree using <span class="artifact-n">[anvi-gen-phylogenomic-tree](/software/anvio/help/7.1/programs/anvi-gen-phylogenomic-tree)</span> or through programs like [FastTree](http://www.microbesonline.org/fasttree/). 

In Anvi'o, this is an output of <span class="artifact-n">[anvi-get-sequences-for-gene-clusters](/software/anvio/help/7.1/programs/anvi-get-sequences-for-gene-clusters)</span> (for generating a tree based off of gene clusters in your workflow) as well as <span class="artifact-n">[anvi-get-sequences-for-hmm-hits](/software/anvio/help/7.1/programs/anvi-get-sequences-for-hmm-hits)</span> (for generating a tree based off of the genes that got HMM hits). 



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/concatenated-gene-alignment-fasta.md) to update this information.


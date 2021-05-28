---
layout: page
title: trnaseq-db [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/DB.png" alt="DB" style="width:100px; border:none" />

A DB-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-trnaseq](../../programs/anvi-trnaseq)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-convert-trnaseq-database](../../programs/anvi-convert-trnaseq-database)</span></p>


## Description

A tRNA-seq database is an anvi'o database that **contains tRNA sequences predicted from a <span class="artifact-n">[trnaseq-fasta](/software/anvio/help/7/artifacts/trnaseq-fasta)</span> for a sample, and associated information**.

This database is the key output of **<span class="artifact-n">[anvi-trnaseq](/software/anvio/help/7/programs/anvi-trnaseq)</span>**, which predicts which input sequences are tRNA, clusters them into discrete biological sequences, and predicts positions in the sequences that are the sites of nucleotide modifications. The database therefore contains tables of information produced by each part of the process.

* Unique sequences predicted to be tRNA, including read counts
* Primary sequence and secondary structural features (stems and loops) predicted in each unique tRNA
* Unconserved nucleotides in the primary sequence that differ from expectation
* Unpaired nucleotides in the stems
* "Trimmed" tRNA sequences, formed from unique sequences only differing by 5' nucleotides beyond the acceptor stem and 3' nucleotides of the CCA acceptor region
* "Normalized" tRNA sequences, formed by dereplicating trimmed tRNA sequences that are 3' fragments produced by incomplete reverse transcription and by mapping biological tRNA fragments
* Potentially modified tRNA sequences, formed by clustering normalized tRNA sequences and retaining those clusters that differ by 3-4 nucleotides at potentially modified positions

This database is the key input to **<span class="artifact-n">[anvi-convert-trnaseq-database](/software/anvio/help/7/programs/anvi-convert-trnaseq-database)</span>**, which takes one or more databases comprising the samples of an experiment and generates a <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> of tRNA seed sequences and <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span>. These can then be displayed and manipulated in anvi'o like other 'omics data.

{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/trnaseq-db.md) to update this information.


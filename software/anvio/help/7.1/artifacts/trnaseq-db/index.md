---
layout: page
title: trnaseq-db [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/trnaseq-db
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/DB.png" alt="DB" style="width:100px; border:none" />

A DB-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-trnaseq](../../programs/anvi-trnaseq)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-merge-trnaseq](../../programs/anvi-merge-trnaseq)</span></p>


## Description

A tRNA-seq database **contains information on tRNA sequences predicted from a single tRNA-seq sample**.

This database is the key output of **<span class="artifact-n">[anvi-trnaseq](/software/anvio/help/7.1/programs/anvi-trnaseq)</span>**. That program predicts which reads are tRNA through structural profiling, clusters tRNA reads into discrete biological sequences, and predicts the positions of nucleotide modifications.

The series of steps implemented in <span class="artifact-n">[anvi-trnaseq](/software/anvio/help/7.1/programs/anvi-trnaseq)</span> sequentially adds the following information to the database.

* Unique sequences predicted to be tRNA, including read counts
* Primary sequence and secondary structural features (stems and loops) predicted in each profiled tRNA
* Unconserved nucleotides in the primary sequence that differ from expectation
* Unpaired nucleotides in the stems
* "Trimmed" tRNA sequences, formed from unique sequences only differing by 3' nucleotides of the CCA acceptor region and 5' nucleotides beyond the acceptor stem
* "Normalized" tRNA sequences, formed by dereplicating trimmed tRNA sequences that are 3' fragments from incomplete reverse transcription and by mapping biological 5' and interior tRNA fragments
* Potentially modified tRNA sequences, formed by clustering normalized tRNA sequences and retaining those clusters that differ by 3-4 nucleotides at aligned positions

This database is the key input to **<span class="artifact-n">[anvi-merge-trnaseq](/software/anvio/help/7.1/programs/anvi-merge-trnaseq)</span>**, which takes one or more databases comprising the samples in an experiment and generates a <span class="artifact-n">[trnaseq-contigs-db](/software/anvio/help/7.1/artifacts/trnaseq-contigs-db)</span> of tRNA seed sequences and <span class="artifact-n">[trnaseq-profile-db](/software/anvio/help/7.1/artifacts/trnaseq-profile-db)</span>s. These tRNA-seq variant contigs and profile databases can then be manipulated and displayed in anvi'o like normal <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>s and <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span>s.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/trnaseq-db.md) to update this information.


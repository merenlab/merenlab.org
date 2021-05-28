---
layout: page
title: trnaseq-fasta [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/FASTA.png" alt="FASTA" style="width:100px; border:none" />

A FASTA-type anvi'o artifact. This artifact is typically provided **by the user** for anvi'o to import into its databases, process, and/or use.

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


There are no anvi'o tools that generate this artifact, which means it is most likely provided to the anvi'o ecosystem by the user.


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-trnaseq](../../programs/anvi-trnaseq)</span></p>


## Description

A <span class="artifact-n">[trnaseq-fasta](/software/anvio/help/7/artifacts/trnaseq-fasta)</span> is a <span class="artifact-n">[fasta](/software/anvio/help/7/artifacts/fasta)</span> file of sequences from a single tRNA-seq sample of split that is suitable to be used by <span class="artifact-n">[anvi-trnaseq](/software/anvio/help/7/programs/anvi-trnaseq)</span> to create a <span class="artifact-n">[trnaseq-db](/software/anvio/help/7/artifacts/trnaseq-db)</span>.

Like <span class="artifact-n">[contigs-fasta](/software/anvio/help/7/artifacts/contigs-fasta)</span> files, this file **is required to have simple deflines**. Take a look at your deflines prior to mapping, and remove anything that is not a digit, an ASCII letter, an underscore, or a dash character. The program <span class="artifact-n">[anvi-script-reformat-fasta](/software/anvio/help/7/programs/anvi-script-reformat-fasta)</span> can do this automatically for you with the flag `--simplify-names`.

We recommend using <span class="artifact-n">[anvi-run-workflow](/software/anvio/help/7/programs/anvi-run-workflow)</span> to create this file from paired-end tRNA-seq reads. The <span class="artifact-n">[trnaseq-workflow](/software/anvio/help/7/artifacts/trnaseq-workflow)</span> uses [illumina-utils](https://github.com/merenlab/illumina-utils) to merge FASTQ files that may contain a mixture of fully and partially overlapping reads, which both occur using 100 bp (or shorter) reads containing barcodes due to the length of tRNA. Even with 150 bp reads, there may be pre-tRNA covered by partially but not fully overlapping reads. <span class="artifact-n">[anvi-script-reformat-fasta](/software/anvio/help/7/programs/anvi-script-reformat-fasta)</span> comes after illumina-utils in the workflow.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/trnaseq-fasta.md) to update this information.


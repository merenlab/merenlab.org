---
layout: page
title: short-reads-fasta [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/FASTA.png" alt="FASTA" style="width:100px; border:none" />

A FASTA-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-get-short-reads-from-bam](../../programs/anvi-get-short-reads-from-bam)</span> <span class="artifact-p">[anvi-get-short-reads-mapping-to-a-gene](../../programs/anvi-get-short-reads-mapping-to-a-gene)</span> <span class="artifact-p">[anvi-script-gen-pseudo-paired-reads-from-fastq](../../programs/anvi-script-gen-pseudo-paired-reads-from-fastq)</span> <span class="artifact-p">[anvi-script-gen-short-reads](../../programs/anvi-script-gen-short-reads)</span> <span class="artifact-p">[anvi-script-get-primer-matches](../../programs/anvi-script-get-primer-matches)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-script-gen-pseudo-paired-reads-from-fastq](../../programs/anvi-script-gen-pseudo-paired-reads-from-fastq)</span></p>


## Description

Similarly to a <span class="artifact-n">[genes-fasta](/software/anvio/help/7/artifacts/genes-fasta)</span>, a short reads fasta is what it sounds like: a <span class="artifact-n">[fasta](/software/anvio/help/7/artifacts/fasta)</span> file containing short reads.

Short reads usually refer to the initial pieces of sequencing data that you had before you assembled them into longer contigs. In other words, these are the kinds of reads you could get out of a technique like Sanger sequencing. Knowing how those short reads align to your contigs is vital for analysis (In fact, that's a lot of the functionality of a <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span>!).

In anvi'o, you can currently get short reads from four sources:

* from a <span class="artifact-n">[bam-file](/software/anvio/help/7/artifacts/bam-file)</span> by running the program <span class="artifact-n">[anvi-get-short-reads-from-bam](/software/anvio/help/7/programs/anvi-get-short-reads-from-bam)</span>
* from a gene in a <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> by running the program <span class="artifact-n">[anvi-get-short-reads-mapping-to-a-gene](/software/anvio/help/7/programs/anvi-get-short-reads-mapping-to-a-gene)</span>
* from something (like a contig) within a <span class="artifact-n">[fasta](/software/anvio/help/7/artifacts/fasta)</span> file by running the program <span class="artifact-n">[anvi-script-get-primer-matches](/software/anvio/help/7/programs/anvi-script-get-primer-matches)</span>
* by generating mock short read data from assembled contigs (as done in the program <span class="artifact-n">[anvi-script-gen-short-reads](/software/anvio/help/7/programs/anvi-script-gen-short-reads)</span>).


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/short-reads-fasta.md) to update this information.


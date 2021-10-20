---
layout: page
title: trnaseq-profile-db [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/trnaseq-profile-db
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/DB.png" alt="DB" style="width:100px; border:none" />

A DB-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-merge-trnaseq](../../programs/anvi-merge-trnaseq)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-tabulate-trnaseq](../../programs/anvi-tabulate-trnaseq)</span></p>


## Description

A tRNA-seq profile database is a **<span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span> variant containing tRNA seed coverage information from one or more samples**.

This database is created by the program, <span class="artifact-n">[anvi-merge-trnaseq](/software/anvio/help/7.1/programs/anvi-merge-trnaseq)</span>, which is part of the <span class="artifact-n">[trnaseq-workflow](/software/anvio/help/7.1/artifacts/trnaseq-workflow)</span>. This program also creates a <span class="artifact-n">[trnaseq-contigs-db](/software/anvio/help/7.1/artifacts/trnaseq-contigs-db)</span>.

## Specific and nonspecific coverage

The coverage of tRNA seeds by tRNA-seq reads is determined differently than the coverage of contigs by metagenomic reads. Metagenomic contigs are constructed by an assembly tool and reads are assigned to contigs by a mapping tool. Reads that map to multiple contigs are randomly assigned to one contig. tRNA-seq seeds and coverages are found simultaneously by <span class="artifact-n">[anvi-trnaseq](/software/anvio/help/7.1/programs/anvi-trnaseq)</span> and <span class="artifact-n">[anvi-merge-trnaseq](/software/anvio/help/7.1/programs/anvi-merge-trnaseq)</span>. Two types of coverage are tracked: **specific** coverage of reads unique to seeds and **nonspecific** coverage of reads in multiple seeds. tRNA-seq reads are often short fragments found in numerous tRNAs; random assignment of these reads would distort tRNA abundances and coverage patterns.

Separate tRNA-seq profile databases are produced for specific and nonspecific coverages. A "combined" database containing both sets of data is produced by default for convenience, allowing specific and nonspecific coverages to be compared side-by-side in the <span class="artifact-n">[anvi-interactive](/software/anvio/help/7.1/programs/anvi-interactive)</span> interface. A "summed" database of specific + nonspecific coverage can optionally be produced. The `--nonspecific-output` option of <span class="artifact-n">[anvi-merge-trnaseq](/software/anvio/help/7.1/programs/anvi-merge-trnaseq)</span> controls the production of nonspecific, combined, and summed databases.

## Modifications versus SNVs

The other significant difference between a tRNA-seq profile database and a normal <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span> is that variable nucleotides are restricted to tRNA modification positions predicted from mutation signatures. Single nucleotide variants are purposefully excluded, though they can be mistaken for modifications, especially with permissive parameterization of <span class="artifact-n">[anvi-merge-trnaseq](/software/anvio/help/7.1/programs/anvi-merge-trnaseq)</span> (see that artifact for more information).

## Uses

Tabulation of tRNA-seq data by <span class="artifact-n">[anvi-tabulate-trnaseq](/software/anvio/help/7.1/programs/anvi-tabulate-trnaseq)</span> takes a specific and optionally nonspecific profile database in addition to a <span class="artifact-n">[trnaseq-contigs-db](/software/anvio/help/7.1/artifacts/trnaseq-contigs-db)</span>.

Interactive visualization of tRNA-seq data in <span class="artifact-n">[anvi-interactive](/software/anvio/help/7.1/programs/anvi-interactive)</span> requires a specific, nonspecific, combined, or summed profile database in addition to a <span class="artifact-n">[trnaseq-contigs-db](/software/anvio/help/7.1/artifacts/trnaseq-contigs-db)</span>.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/trnaseq-profile-db.md) to update this information.


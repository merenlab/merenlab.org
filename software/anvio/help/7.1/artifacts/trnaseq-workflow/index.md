---
layout: page
title: trnaseq-workflow [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/trnaseq-workflow
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/WORKFLOW.png" alt="WORKFLOW" style="width:100px; border:none" />

A WORKFLOW-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-run-workflow](../../programs/anvi-run-workflow)</span></p>


## Required or used by


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

The tRNA-seq workflow **is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow run by the anvi'o script, <span class="artifact-n">[anvi-run-workflow](/software/anvio/help/7.1/programs/anvi-run-workflow)</span>**.

The workflow can run the following programs in order.
- [Illumina-utils](https://github.com/merenlab/illumina-utils), for merging paired-end reads and quality control
- <span class="artifact-n">[anvi-script-reformat-fasta](/software/anvio/help/7.1/programs/anvi-script-reformat-fasta)</span>, for making FASTA deflines anvio-compliant
- <span class="artifact-n">[anvi-trnaseq](/software/anvio/help/7.1/programs/anvi-trnaseq)</span>, for predicting tRNA sequences, structures, and modification sites in each sample
- <span class="artifact-n">[anvi-merge-trnaseq](/software/anvio/help/7.1/programs/anvi-merge-trnaseq)</span>, for predicting tRNA seed sequences and their modification sites from the set of samples
- <span class="artifact-n">[anvi-run-trna-taxonomy](/software/anvio/help/7.1/programs/anvi-run-trna-taxonomy)</span>, for assigning taxonomy to tRNA seeds
- <span class="artifact-n">[anvi-tabulate-trnaseq](/software/anvio/help/7.1/programs/anvi-tabulate-trnaseq)</span>, for generating tables of seed and modification information that are easily manipulated

## Input

The tRNA-seq workflow requires two files to run: a json-formatted config file and <span class="artifact-n">[samples-txt](/software/anvio/help/7.1/artifacts/samples-txt)</span>. Generate the default config file, here called `config.json`, with the following command.

<div class="codeblock" markdown="1">
anvi&#45;run&#45;workflow &#45;w trnaseq &#45;&#45;get&#45;default&#45;config config.json
</div>

Different "rules," or steps, of the workflow can be turned on and off as needed in the config file. The workflow can be restarted at intermediate rules without rerunning prior rules that have already completed.

<span class="artifact-n">[samples-txt](/software/anvio/help/7.1/artifacts/samples-txt)</span> will contain a list of FASTQ or FASTA files and associated information on each library. FASTQ files contain unmerged paired-end tRNA-seq reads. Reads are merged in the workflow by [Illumina-utils](https://github.com/merenlab/illumina-utils). FASTA files contain merged reads, and the initial read-merging steps in the workflow are skipped.

Here is an example tRNA-seq samples file with FASTQ inputs.

| sample | treatment | r1 | r2 | r1_prefix | r2_prefix |
| --- | --- | --- | --- | --- | --- |
| ecoli_A1_noDM | untreated | FASTQ/ecoli_A1_noDM.r1.fq.gz | FASTQ/ecoli_A1_noDM.r2.fq.gz | NNNNNN | TTCCAGT |
| ecoli_A1_DM | demethylase | FASTQ/ecoli_A1_DM.r1.fq.gz | FASTQ/ecoli_A1_DM.r2.fq.gz | NNNNNN | TCTGAGT |
| ecoli_B1_noDM | untreated | FASTQ/ecoli_B1_noDM.r1.fq.gz | FASTQ/ecoli_B1_noDM.r2.fq.gz | NNNNNN | TGGTAGT |
| ecoli_B1_DM | demethylase | FASTQ/ecoli_B1_DM.r1.fq.gz | FASTQ/ecoli_B1_DM.r2.fq.gz | NNNNNN | CTGAAGT |

The treatment column is optional. The treatment indicates a chemical application, such as demethylase, and can be used to have a bearing on seed sequence determination in <span class="artifact-n">[anvi-merge-trnaseq](/software/anvio/help/7.1/programs/anvi-merge-trnaseq)</span>. In the absence of a treatment column, all samples are assigned the same treatment, which can be specified in the `anvi_trnaseq` section of the workflow config file and defaults to `untreated`.

Read 1 and 2 prefix columns are also optional. These represent sequences that Illumina-utils should identify and trim from the start of the read. In the example, the read 1 prefix is a unique molecular identifier (UMI) of 6 random nucleotides, and the read 2 prefix is a sample barcode. Illumina-utils will discard the paired-end read if the prefix is not found. In the example, the read 1 UMI will always be found, but the read 2 barcode must match exactly.

Here is an equivalent tRNA-seq samples file with FASTA inputs.

| sample | treatment | fasta |
| --- | --- | --- |
| ecoli_A1_noDM | untreated | FASTA/ecoli_A1_noDM.fa.gz |
| ecoli_A1_DM | demethylase | FASTA/ecoli_A1_DM.fa.gz |
| ecoli_B1_noDM | untreated | FASTA/ecoli_B1_noDM.fa.gz |
| ecoli_B1_DM | demethylase | FASTA/ecoli_B1_DM.fa.gz |

Note that barcodes and other sequence prefixes should already be trimmed from FASTA sequences.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/trnaseq-workflow.md) to update this information.


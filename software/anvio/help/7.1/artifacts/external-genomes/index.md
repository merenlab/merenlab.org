---
layout: page
title: external-genomes [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/external-genomes
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact is typically provided **by the user** for anvi'o to import into its databases, process, and/or use.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-script-gen-genomes-file](../../programs/anvi-script-gen-genomes-file)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-compute-functional-enrichment-across-genomes](../../programs/anvi-compute-functional-enrichment-across-genomes)</span> <span class="artifact-r">[anvi-compute-genome-similarity](../../programs/anvi-compute-genome-similarity)</span> <span class="artifact-r">[anvi-compute-metabolic-enrichment](../../programs/anvi-compute-metabolic-enrichment)</span> <span class="artifact-r">[anvi-dereplicate-genomes](../../programs/anvi-dereplicate-genomes)</span> <span class="artifact-r">[anvi-display-functions](../../programs/anvi-display-functions)</span> <span class="artifact-r">[anvi-estimate-genome-completeness](../../programs/anvi-estimate-genome-completeness)</span> <span class="artifact-r">[anvi-estimate-metabolism](../../programs/anvi-estimate-metabolism)</span> <span class="artifact-r">[anvi-gen-genomes-storage](../../programs/anvi-gen-genomes-storage)</span> <span class="artifact-r">[anvi-get-sequences-for-hmm-hits](../../programs/anvi-get-sequences-for-hmm-hits)</span> <span class="artifact-r">[anvi-script-gen-functions-per-group-stats-output](../../programs/anvi-script-gen-functions-per-group-stats-output)</span> <span class="artifact-r">[anvi-script-gen-hmm-hits-matrix-across-genomes](../../programs/anvi-script-gen-hmm-hits-matrix-across-genomes)</span></p>


## Description

An external genome is any genome assembly that was converted into a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> from its original FASTA file format using the program <span class="artifact-n">[anvi-gen-contigs-database](/software/anvio/help/7.1/programs/anvi-gen-contigs-database)</span>. You can obtain one of these in a variety of ways, the most common being 1) downloading a genome from a database such as NCBI and 2) assembling a genome yourself from sequencing reads. The key thing is that the sequences in the <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> represent a _single_ microbial population (or species, if you are not working with microbes) - ie, it is not a metagenome.

The external genomes file format enables anvi'o to work with one or more external genomes. A TAB-delimited external genomes file will be composed of at least the following two columns:

|name|contigs_db_path|
|:--|:--|
|Name_01|/path/to/contigs-01.db|
|Name_02|/path/to/contigs-02.db|
|Name_03|/path/to/contigs-03.db|
|(...)|(...)|

{:.warning}
Please make sure names in the `name` column does not include any special characters (underscore is fine). It is also a good idea to keep these names short and descriptive as they will appear in various figures in downstream analyses.

Also see **<span class="artifact-n">[internal-genomes](/software/anvio/help/7.1/artifacts/internal-genomes)</span>** and **<span class="artifact-n">[metagenomes](/software/anvio/help/7.1/artifacts/metagenomes)</span>**.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/external-genomes.md) to update this information.


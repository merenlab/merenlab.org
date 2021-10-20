---
layout: page
title: fasta-txt [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/fasta-txt
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact is typically provided **by the user** for anvi'o to import into its databases, process, and/or use.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


There are no anvi'o tools that generate this artifact, which means it is most likely provided to the anvi'o ecosystem by the user.


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-run-workflow](../../programs/anvi-run-workflow)</span></p>


## Description

This is a file used by <span class="artifact-n">[anvi-run-workflow](/software/anvio/help/7.1/programs/anvi-run-workflow)</span> that lists the name and path of all of the input <span class="artifact-n">[fasta](/software/anvio/help/7.1/artifacts/fasta)</span> files.

As of now, this file is used in the <span class="artifact-n">[contigs-workflow](/software/anvio/help/7.1/artifacts/contigs-workflow)</span>, <span class="artifact-n">[pangenomics-workflow](/software/anvio/help/7.1/artifacts/pangenomics-workflow)</span>, and [the reference mode](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#references-mode) of the <span class="artifact-n">[metagenomics-workflow](/software/anvio/help/7.1/artifacts/metagenomics-workflow)</span>.

In its simplest form, a <span class="artifact-n">[fasta-txt](/software/anvio/help/7.1/artifacts/fasta-txt)</span> is a TAB-delmited file with two columns for `name` and `path`. Here is an example:

|name|path|
|:--|:--|
|SAMPLE_01|path/to/sample_01.fa|
|SAMPLE_02|path/to/sample_02.fa|

Paths can be absolute or relative, and FASTA files can be compressed or not. That's all up to you.

One of the primary users of the <span class="artifact-n">[fasta-txt](/software/anvio/help/7.1/artifacts/fasta-txt)</span> is the [anvi'o snakemake workflows](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/), and to make it more compatible to complex workflow needs, <span class="artifact-n">[fasta-txt](/software/anvio/help/7.1/artifacts/fasta-txt)</span> supports the following additional columns to provide more information for each FASTA file when available, such as <span class="artifact-n">[external-gene-calls](/software/anvio/help/7.1/artifacts/external-gene-calls)</span> file and/or a <span class="artifact-n">[functions-txt](/software/anvio/help/7.1/artifacts/functions-txt)</span>.

Here is an example with those additional columns:

|name|path|external_gene_calls|gene_functional_annotation|
|:--|:--|:--|:--|
|SAMPLE_01|path/to/sample_01.fa|<span class="artifact-n">[external-gene-calls](/software/anvio/help/7.1/artifacts/external-gene-calls)</span>_01.txt|<span class="artifact-n">[functions-txt](/software/anvio/help/7.1/artifacts/functions-txt)</span>_01.txt|
|SAMPLE_02|path/to/sample_02.fa|<span class="artifact-n">[external-gene-calls](/software/anvio/help/7.1/artifacts/external-gene-calls)</span>_02.txt|<span class="artifact-n">[functions-txt](/software/anvio/help/7.1/artifacts/functions-txt)</span>_02.txt|

For more information, check out the [anvi'o workflow tutorial](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#fastatxt)


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/fasta-txt.md) to update this information.


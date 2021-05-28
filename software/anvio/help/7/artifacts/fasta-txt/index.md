---
layout: page
title: fasta-txt [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact is typically provided **by the user** for anvi'o to import into its databases, process, and/or use.

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


There are no anvi'o tools that generate this artifact, which means it is most likely provided to the anvi'o ecosystem by the user.


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-run-workflow](../../programs/anvi-run-workflow)</span></p>


## Description

This is a file used by <span class="artifact-n">[anvi-run-workflow](/software/anvio/help/7/programs/anvi-run-workflow)</span> that lists the name and path of all of the input <span class="artifact-n">[fasta](/software/anvio/help/7/artifacts/fasta)</span> files. 

As of now, this file is used in the <span class="artifact-n">[contigs-workflow](/software/anvio/help/7/artifacts/contigs-workflow)</span>, <span class="artifact-n">[pangenomics-workflow](/software/anvio/help/7/artifacts/pangenomics-workflow)</span>, and [the reference mode](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#references-mode) of the <span class="artifact-n">[metagenomics-workflow](/software/anvio/help/7/artifacts/metagenomics-workflow)</span>.

This file will look something like this: 

    name        path
    SAMPLE_01   path/to/sample_01.fa
    SAMPLE_02   path/to/sample_02.fa
    
Note that the paths can be either absolute or relative, and the fasta files can be either compressed or not. That's all up to you. 

To input more information, for each file you can also specify an <span class="artifact-n">[external-gene-calls](/software/anvio/help/7/artifacts/external-gene-calls)</span> file and/or a <span class="artifact-n">[functions-txt](/software/anvio/help/7/artifacts/functions-txt)</span>. Just provide those with additional columns, as so: 

    name        path                    external_gene_calls             gene_functional_annotation
    SAMPLE_01   path/to/sample_01.fa    <span class="artifact-n">[external-gene-calls](/software/anvio/help/7/artifacts/external-gene-calls)</span>_01.txt  <span class="artifact-n">[functions-txt](/software/anvio/help/7/artifacts/functions-txt)</span>_01.txt
    SAMPLE_02   path/to/sample_02.fa    <span class="artifact-n">[external-gene-calls](/software/anvio/help/7/artifacts/external-gene-calls)</span>_02.txt  <span class="artifact-n">[functions-txt](/software/anvio/help/7/artifacts/functions-txt)</span>_02.txt

For more information, check out the [anvi'o workflow tutorial](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#fastatxt)


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/fasta-txt.md) to update this information.


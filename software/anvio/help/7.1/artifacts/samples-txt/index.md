---
layout: page
title: samples-txt [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/samples-txt
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


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-run-workflow](../../programs/anvi-run-workflow)</span> <span class="artifact-r">[anvi-script-get-primer-matches](../../programs/anvi-script-get-primer-matches)</span></p>


## Description

A **TAB-delimited** file to describe samples and paired-end FASTQ files associated with them. By doing so, this file type links sample names to raw sequencing reads.

This file type includes required and optional columns.

The following three columns are **required** for this file type:

* `sample`: a single-word sample name,
* `r1`: path to the FASTQ file for pair one, and
* `r2`: path to the FASTQ file for pair two.

While you can use relative paths for `r1` and `r2`, it is always better to have absolute paths to improve reproducibility.

The following is an **optional** column:

* `group`: A single-word categorical variable that assigns two or more samples into two or more groups. This is useful to co-assemble multiple samples so that you can bin them later. 

For more information, see the [anvi'o workflow tutorial](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#samplestxt)

### Examples samples.txt file

Here is an example file:

|sample|group|r1|r2|
|:--|:--|:--|:--|
|Sample_01|WARM|/path/to/XXX-01-R1.fastq.gz|/path/to/XXX-01-R2.fastq.gz|
|Sample_02|COLD|/path/to/YYY-02-R1.fastq.gz|/path/to/YYY-02-R2.fastq.gz|
|Sample_03|COLD|/path/to/ZZZ-03-R1.fastq.gz|/path/to/ZZZ-03-R2.fastq.gz|


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/samples-txt.md) to update this information.


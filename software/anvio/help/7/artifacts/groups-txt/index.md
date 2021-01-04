---
layout: page
title: groups-txt [artifact]
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


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-compute-functional-enrichment](../../programs/anvi-compute-functional-enrichment)</span></p>


## Description

This is a two-column, TAB-delimited file that describes which samples belong to which groups. It is a required input for computing enrichment scores of metabolic modules in <span class="artifact-n">[anvi-compute-functional-enrichment](/software/anvio/help/7/programs/anvi-compute-functional-enrichment)</span>.

Here is an example:

|sample|group|
|:--|:--|
|sample_01|GROUP-A|
|sample_02|GROUP-B|
|sample_03|GROUP-A|
|(...)|(...)|

{:.warning}
The names in the `sample` column must match those in the "modules" mode output file that you provide to the <span class="artifact-n">[anvi-compute-functional-enrichment](/software/anvio/help/7/programs/anvi-compute-functional-enrichment)</span> program. If you know that the sample names match but you are still getting errors, you might need to specify which column in the "modules" mode output contains those sample names using the `--sample-header` parameter of <span class="artifact-n">[anvi-compute-functional-enrichment](/software/anvio/help/7/programs/anvi-compute-functional-enrichment)</span>. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/groups-txt.md) to update this information.


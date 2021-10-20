---
layout: page
title: quick-summary [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/quick-summary
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-summarize-blitz](../../programs/anvi-summarize-blitz)</span></p>


## Required or used by


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

The output of <span class="artifact-n">[anvi-summarize-blitz](/software/anvio/help/7.1/programs/anvi-summarize-blitz)</span>.

<span class="artifact-n">[anvi-summarize-blitz](/software/anvio/help/7.1/programs/anvi-summarize-blitz)</span> summarizes read-recruitment statistics for a collection of bins across multiple samples. It produces long format output in which each row contains the (weighted) average statistics of a bin in a sample. Each statistic is summarized in a different column of the file.

Here is an example output file from this program, summarizing `detection` and `mean_coverage_Q2Q3` data (the default statistics) for 3 bins across multiple samples:

unique_id | bin_name | sample | detection | mean_coverage_Q2Q3
|:---|:---|:---|:---|:---|
0 | bin_1 | sample_1 | 0.015553023620503776 | 1.0272713907674214
1 | bin_2 | sample_1 | 0.0004871607502275562 | 0.0
2 | bin_3 | sample_1 | 0.0023636043452898497 | 0.0
3 | bin_1 | sample_2 | 0.015767421346662747 | 1.1101759286484367
4 | bin_2 | sample_2 | 0.0004871607502275562 | 0.0
5 | bin_3 | sample_2 | 0.001595914458984989 | 0.0
[...] | [...] |[...] |[...] |[...]

The `unique_id` column is just a unique index for each row. Each column after the `sample` column contains a different statistic (to learn how to include different or additional statistics in this output, read the <span class="artifact-n">[anvi-summarize-blitz](/software/anvio/help/7.1/programs/anvi-summarize-blitz)</span> page.)


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/quick-summary.md) to update this information.


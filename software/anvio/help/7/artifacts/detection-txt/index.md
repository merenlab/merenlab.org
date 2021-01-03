---
layout: page
title: detection-txt [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-export-gene-coverage-and-detection](../../programs/anvi-export-gene-coverage-and-detection)</span></p>


## Required or used by


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

This is a text file containing **the detection value for each gene in each sample** that was in the <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span> and <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> that you used when you ran  <span class="artifact-n">[anvi-export-gene-coverage-and-detection](/software/anvio/help/7/programs/anvi-export-gene-coverage-and-detection)</span>. 

This is a tab-delimited file where each row describes a specific gene and each column describes one of your samples. Each cell contains the detection of that gene in that sample. 

### Example

    key       sample_1    sample_2    sample_3 ...
    13947     0.291093    0.984394    0.9289432         
    13948     0.895842    0.828481    0.3721947
    23026     0.949383    0.983923    1.0000000
    ...






{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/detection-txt.md) to update this information.


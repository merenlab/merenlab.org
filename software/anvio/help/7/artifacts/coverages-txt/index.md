---
layout: page
title: coverages-txt [artifact]
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


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-export-gene-coverage-and-detection](../../programs/anvi-export-gene-coverage-and-detection)</span> <span class="artifact-p">[anvi-export-splits-and-coverages](../../programs/anvi-export-splits-and-coverages)</span> <span class="artifact-p">[anvi-get-split-coverages](../../programs/anvi-get-split-coverages)</span> <span class="artifact-p">[anvi-script-get-coverage-from-bam](../../programs/anvi-script-get-coverage-from-bam)</span></p>


## Required or used by


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

This is a text file containing **the average coverage for each contig in each sample** that was in the <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span> and <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> that you used when you ran <span class="artifact-n">[anvi-export-splits-and-coverages](/software/anvio/help/7/programs/anvi-export-splits-and-coverages)</span> or <span class="artifact-n">[anvi-export-gene-coverage-and-detection](/software/anvio/help/7/programs/anvi-export-gene-coverage-and-detection)</span>. 

This is a tab-delimited file where each row describes a specific split/gene and each column describes one of your samples. Each cell contains the average coverage of that contig in that sample. 

This artifact is really only used when taking information out of anvi'o, so enjoy your coverage information :) 

### Example for splits

(the type of output you would get from <span class="artifact-n">[anvi-export-splits-and-coverages](/software/anvio/help/7/programs/anvi-export-splits-and-coverages)</span>)

    contig                  sample_1    sample_2    sample_3 ...
    Day1_contig1_split1     5.072727    4.523432    1.2343243         
    Day1_contig1_split2     6.895844    5.284812    9.3721947
    Day1_contig2_split1     2.357049    3.519150    8.2385691
    ...


### Example for genes

(the type of output you would get from <span class="artifact-n">[anvi-export-gene-coverage-and-detection](/software/anvio/help/7/programs/anvi-export-gene-coverage-and-detection)</span>)

    key       sample_1    sample_2    sample_3 ...
    13947     10.29109    1.984394    6.8289432         
    13948     34.89584    6.284812    3.3721947
    23026     23.94938    9.239235    13.238569
    ...






{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/coverages-txt.md) to update this information.


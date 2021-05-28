---
layout: page
title: vcf [artifact]
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


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-script-variability-to-vcf](../../programs/anvi-script-variability-to-vcf)</span></p>


## Required or used by


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

This represents a file in the Variant Call Format, which is a standard format for storing sequence variations like SNVs. 

You can convert the information in a <span class="artifact-n">[variability-profile-txt](/software/anvio/help/7/artifacts/variability-profile-txt)</span> to <span class="artifact-n">[vcf](/software/anvio/help/7/artifacts/vcf)</span> with the program <span class="artifact-n">[anvi-script-variability-to-vcf](/software/anvio/help/7/programs/anvi-script-variability-to-vcf)</span>. 

### What's in this file? 

For more details, you can check out the [VCF wikipedia page](https://en.wikipedia.org/wiki/Variant_Call_Format). 

#### Header

Briefly, this file's header (marked by `##` at the beginning of each line) contains various metadata. This includes the date, link to the reference file, contig information, etc. It also contains 
- what information will be reported (denoted by `INFO`). 
- what additional filters will be run on each SNV (denoted by `FILTER`). For example, marking which variants are below a certain quality threshold. 
- what format to display additional data in (denoted by `FORMAT`). 

#### Body

The body of the file contains identifying information for the variation (the chromosome, position and ID), the identity of the position in the reference and alternative alleles present in your data. Following this is a quality score for your data and the additional information specified by the header. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/vcf.md) to update this information.


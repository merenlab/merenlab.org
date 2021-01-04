---
layout: page
title: pn-ps-data [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/CONCEPT.png" alt="CONCEPT" style="width:100px; border:none" />

A CONCEPT-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-script-calculate-pn-ps-ratio](../../programs/anvi-script-calculate-pn-ps-ratio)</span></p>


## Required or used by


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

This describes the output of <span class="artifact-n">[anvi-script-calculate-pn-ps-ratio](/software/anvio/help/7/programs/anvi-script-calculate-pn-ps-ratio)</span>, which calculates the pN/pS ratio for each gene in a <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span>. 

{:.notice}
See the page for <span class="artifact-n">[anvi-script-calculate-pn-ps-ratio](/software/anvio/help/7/programs/anvi-script-calculate-pn-ps-ratio)</span> for an explanation of the pN/pS ratio 

This describes a directory that contains the following three files: 

`pN_pS_ratio.txt`: a matrix of the pN/pS values where each row is a gene and each column is a metagenome. 

    gene_caller_id  ANE_004_05M             ANE_004_05M             ANE_150_05M           ...
    1248            0.024404785965479868    0.024404785965479868    0.019650249960943812
    1249            0.014399205561072536    0.014399205561072536    0.0
    ...

`SAAV_counts.txt`: a matrix that describes the number of SAAVs for each gene in each metagenome.  
        
    gene_caller_id  ANE_004_05M     ANE_004_05M     ANE_150_05M      ...
    1248            12              11              10    
    1249            1               1               0
    ...
        
`SCV_counts.txt`: a matrix that describes the number of SCVs for each gene in each metagenome.  

    gene_caller_id  ANE_004_05M     ANE_004_05M     ANE_150_05M      ...
    1248            143             154             148    
    1249            19              19              14
    ...


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/pn-ps-data.md) to update this information.


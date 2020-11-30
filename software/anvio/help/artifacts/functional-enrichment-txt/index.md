---
layout: page
title: functional-enrichment-txt [artifact]
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


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-get-enriched-functions-per-pan-group](../../programs/anvi-get-enriched-functions-per-pan-group)</span></p>


## Required or used by


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

This matrix describes the gene cluster level functional associations that are enriched within specific groups of your pangenome. It is the output of <span class="artifact-n">[anvi-get-enriched-functions-per-pan-group](/software/anvio/help/programs/anvi-get-enriched-functions-per-pan-group)</span> and is described in more detail [in the pangenomics tutorial](http://merenlab.org/2016/11/08/pangenomics-v2/#making-sense-of-functions-in-your-pangenome). 

Each row in this matrix describes a functional assocation that is associated with one or more of your gene clusters. These are listed in order of their pan-group, with the highest enrichment scores displayed first. 

The following columns of information are listed 

- your functional annotation source (as the column header): describes the functional assocation that this row is about. 
- enrichment_score: a measure of much this particular function is enriched in the pan-group it is associated with (i.e., measures how unique this function [see column 1] is to this group [see column 5]) 
- unadjusted_p_value 
- adjusted_q_value 
- associated groups: the list of pan groups that this function is associated with (for example, the low light group in the pangenomic tutorial)
- function_accession number 
- gene_cluster_ids: lists the gene clusters that are included in this functional assocation. 
-  p values for each group: gives the portion of the group member genomes where this functional association was found. 
-N values for each group: gives the total number of genomes in each group. 

Here is a more concrete example (the same example as the [pangenomic tutorial](http://merenlab.org/2016/11/08/pangenomics-v2/#making-sense-of-functions-in-your-pangenome)). Note that that tutorial uses `COG_FUNCTION` as the functional annotation source, and has `LL` (low light) and `HL` (high light) as the two pan-groups. 

|COG_FUNCTION | enrichment_score | unadjusted_p_value | adjusted_q_value | associated_groups | function_accession | gene_clusters_ids | p_LL | p_HL | N_LL | N_HL|
|-- | -- | -- | -- | -- | -- | -- | -- | -- | --| --|
|Proteasome lid subunit RPN8/RPN11, contains Jab1/MPN domain metalloenzyme (JAMM) motif | 31.00002279 | 2.58E-08 | 1.43E-06 | LL | COG1310 | GC_00002219, GC_00003850, GC_00004483 | 1 | 0 | 11 | 20|
|Adenine-specific DNA glycosylase, acts on AG and A-oxoG pairs | 31.00002279 | 2.58E-08 | 1.43E-06 | LL | COG1194 | GC_00001711 | 1 | 0 | 11 | 20|
|Periplasmic beta-glucosidase and related glycosidases | 31.00002279 | 2.58E-08 | 1.43E-06 | LL | COG1472 | GC_00002086, GC_00003909 | 1 | 0 | 11 | 20|
|Single-stranded DNA-specific exonuclease, DHH superfamily, may be involved in archaeal DNA replication intiation | 31.00002279 | 2.58E-08 | 1.43E-06 | LL | COG0608 | GC_00002752, GC_00003786, GC_00004838, GC_00007241 | 1 | 0 | 11 | 20|
|Ser/Thr protein kinase RdoA involved in Cpx stress response, MazF antagonist | 31.00002279 | 2.58E-08 | 1.43E-06 | LL | COG2334 | GC_00002783, GC_00003936, GC_00004631, GC_00005468 | 1 | 0 | 11 | 20|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|
|Signal transduction histidine kinase | -7.34E-41 | 1 | 1 | NA | COG5002 | GC_00000773, GC_00004293 | 1 | 1 | 11 | 20|
|tRNA A37 methylthiotransferase MiaB | -7.34E-41 | 1 | 1 | NA | COG0621 | GC_00000180, GC_00000851 | 1 | 1 | 11 | 20|


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/functional-enrichment-txt.md) to update this information.


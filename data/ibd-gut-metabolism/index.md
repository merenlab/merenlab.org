---
layout: page
title: A reproducible workflow for Veseli et al, 2023
modified: 2023-04-15
excerpt: "A bioninformatics workflow for our study on microbial metabolism in the IBD gut environment."
comments: true
authors: [iva]
---

FIXME: Introduction with link to pre-print and summary of figshare/zenodo links


## Study description

## Obtaining our dataset of public fecal metagenomes

### Criteria for sample selection and sample groups

### Downloading the public metagenomes used in this study

### Metagenome processing: single assemblies and annotations


## Intermediate metagenome processing

### Analyzing number of populations per sample

### Generating Figure 1: 
scatterplot of sequencing depth vs number of populations

### Removal of samples with low-sequencing depth

### Final set of samples and their contigs DBs


## Metabolism analyses (metagenomes)

### Metabolism estimation

### Normalization of pathway copy numbers to PPCN

### Enrichment analysis for IBD-enriched pathways

### Generating Figure 2


## Obtaining a dataset of gut microbial genomes from the GTDB

### Genome processing: the (snakemake) contigs workflow

### Using the EcoPhylo workflow for quick identification of relevant gut microbes
starting from 3 phyla of gut microbes

### Subsetting gut genomes by detection in our sample groups

### Read recruitment from metagenome samples to gut microbes

### Percent abundance calculations

### Metabolism estimation in gut genomes

### The 'HMI score': Classifying genomes by level of metabolic independence

### Generating the genome phylogeny

### Generating Figure 3


## Machine learning analyses

### Classifying gut metagenomes as IBD vs Healthy

### Classifying antibiotic time-series metagenomes from Palleja et al

### Generating Figure 4


## Supplementary Analyses

### Exploring annotation efficiency (SF 2)

### Additional comparisons of metabolic pathways (SF 3)

### Examining cohort effect (SF 4)

### Testing classifier generalizability
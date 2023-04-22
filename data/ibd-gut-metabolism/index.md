---
layout: page
title: A reproducible workflow for Veseli et al, 2023
modified: 2023-04-15
excerpt: "A bioinformatics workflow for our study on microbial metabolism in the IBD gut environment."
comments: true
authors: [iva]
---

<div class="extra-info" markdown="1">

<span class="extra-info-header">Summary</span>

**The purpose of this page** is to provide access to reproducible data products and analyses for the study "**High metabolic independence is a determinant of microbial resilience in the face of gut stress**" by Veseli et al.

Here is a list of links for quick access tohe data described in our manuscript and on this page:

* [doi:10.6084/m9.figshare.22679080](https://doi.org/10.6084/m9.figshare.22679080): Supplementary Tables.


</div>

{:.notice}
If you have any questions and/or if you are unable to find an important piece of information here, please feel free to leave a comment down below, send an e-mail to [us]({{ site.url }}/people/), or get in touch with us through Discord:

{% include _join-anvio-discord.html %}


## Study description

This study is a follow-up to our [previous study](https://doi.org/10.1186/s13059-023-02924-x) on microbial colonization of the gut following fecal microbiota transplant, in which we introduced the concept of high metabolic independence as a determinant of microbial resilience for populations colonizing new individuals or living in individuals with inflammatory bowel disease (IBD). In the current work, we sought to (1) confirm our prior observations in a high-throughput comparative analysis of the gut microbiomes of healthy individuals and individuals with IBD, and (2) demonstrate that high metabolic independence is a robust marker of general stress experienced by the gut microbiome. To do this, we

* created single assemblies of a large dataset of **publicly*available fecal metagenomes**
* computed the community*level copy numbers of metabolic pathways in each sample, and **normalized** these copy numbers with the estimated number of populations in each sample to obtain **per-population copy numbers (PPCNs)**
* determined which **metabolic pathways are enriched** in the IBD sample group
* identified **bacterial reference genomes** associated with the human gut environment
* scored the reference genomes as having **high or low metabolic independence** based upon completeness scores of the IBD-enriched metabolic pathways in each genome
* analyzed the **distribution of each group of genomes** within the healthy fecal metagenomes and those from individuals with IBD
* **trained a machine learning classifier** with the IBD-enriched pathway PPCN data to differentiate between metagenome samples from individuals with IBD and those from healthy individuals
* tested the classifier on an independent time-series dataset **tracking decline and recovery** of the gut microbiome following **antibiotic treatment** (it worked quite well!)

This webpage offers access to the details of our computational methods (occasionally with helpful context not given in the methods section of the manuscript) for the steps outlined above as well as to the datasets needed to reproduce our work. The workflow is organized into several large sections, each of which covers a set of related steps.

## Obtaining our initial dataset of public fecal metagenomes

This section covers the steps for acquiring and processing our initial set of publicly-available gut metagenomes. We downloaded, assembled and annotated 2,893 samples from 13 different studies, in order to evaluate our metabolic competency hypothesis across a wide diversity of cohorts from different geographical locations, age groups, hospital systems, and degress of healthiness. Note that this extensive dataset was later filtered to remove samples with low-sequencing depth (as described in the next section), and as a result, not all of these samples were utilized for the main analyses in our study. However, you can access the full list of the 2,893 samples that we considered in sheet (c) of Supplementary Table 1 (FIXME: link). 

### Criteria for sample selection and sample groups
We sought to obtain a large number of fecal metagenomes from healthy individuals and from individuals with IBD. We used the following criteria to search for studies offering such samples:

1) The study provided publicly-available shotgun metagenomes of fecal matter. These were usually stool metagenomes, but we also accepted luminal aspirate samples from the ileal pouch of patients that have undergone a colectomy with an ileal pouch anal anastomosis (from [Vineis et al. 2016](https://doi.org/10.1128/mBio.01713-16)), since these samples also represent fecal matter.
2) The study sampled from people living in industrialized countries. These countries have a higher incidence of IBD and share dietary and lifestyle tendencies that result in a similar gut microbiome composition when compared to developing countries. We found one study that sampled from both industrialized and developing countries ([Rampelli et al. 2015](https://doi.org/10.1016/j.cub.2015.04.055)), and from this study we used only those samples from industrialized areas.
3) The study included samples from people with IBD and/or they included samples from people without gastrointestinal (GI) disease or inflammation. The latter were often healthy controls from studies of diseases besides IBD, or from studies of treatments such as dietary interventions and antibiotics, and in these cases we included only the control samples in our dataset.
4) The study provided clear metadata for differentiating between case and control samples, so that we could accurately assign samples to the appropriate group.

These criteria led us to 13 studies of the human gut microbiome, which are summarized in sheet (a) of Supplementary Table 1 (FIXME: link). We almost certainly did not find all possible studies that fit our requirements, but we were sufficiently satisfied with the large number of samples and the breadth of human diversity encompassed by these studies, so we stopped there.

Each of the samples had to be assigned to a group based upon the general health status of the sample donor. In order to do this, we had to make some decisions about what to include (or not) in our characterization of "healthy" people. Not being clinical GI experts, we did our best with the metadata and cohort descriptions that were provided by each study (though there is room for disagreement here). We decided upon three groups: a healthy group of samples from people without gastrointestinal disease or inflammation, which contained the control samples from most of the studies; an IBD group of samples from people with a confirmed diagnosis of ulcerative colitis (UC), Crohn's disease (CD), or unclassified IBD; and an intermediate 'non-IBD' group of samples from people without a definite IBD diagnosis but who nevertheless may be presenting symptoms of GI distress or inflammation. [Lloyd-Price et al. 2019](https://doi.org/10.1038/s41586-019-1237-9) describes the criteria and justification for this last group quite eloquently:

<blockquote>
Subjects not diagnosed with IBD based on endoscopic and histopathologic findings were classified as ‘non-IBD’ controls, including the aforementioned healthy individuals presenting for routine screening, and those with more benign or non-specific symptoms. This creates a control group that, while not completely ‘healthy’, differs from the IBD cohorts specifically by clinical IBD status.

<div class="blockquote-author">Lloyd-Price et al. 2019</div>
</blockquote>

For the studies characterizing their controls as 'non-IBD' samples (there were 3), we applied the same grouping of their samples within our dataset in order to be consistent. We also decided to include samples from people with colorectal cancer (CRC) (from [Feng et al. 2015](https://doi.org/10.1038/ncomms7528)) in the 'non-IBD' group, since CRC represents a source of potential inflammation in the GI tract.

We note that we did not exclude samples from individuals with high BMI from the healthy group, because the study that included such samples ([Le Chatelier et al. 2013](https://doi.org/10.1038/nature12506)) had already excluded individuals with GI disease, diabetes, and other such conditions.

### Downloading the public metagenomes used in this study

### Metagenome processing: single assemblies and annotations


## Selecting our final dataset of metagenomes

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
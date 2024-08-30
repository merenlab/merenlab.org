---
layout: page
title: A reproducible workflow for Meyer & Eren manuscript 2024
modified: 2024-08-30
excerpt: "A bioinformatics workflow for our metagenomics data integration across observatories"
comments: true
authors: [raissa]
---

<div class="extra-info" markdown="1">

<span class="extra-info-header">Summary</span>

**The purpose of this page** is to provide access to reproducible data products and analyses for the manuscript "**Marine Microbial Observatories for the Future: From Samples to Data to Legacy using integrated omics strategies**" by Meyer & Eren


In this study, we 
-  integrate metagenomics metadata and data from observatories (Hawaii Ocean Time Series and Bermuda Atlantic Time-series Study), sampling expeditions (Bio-GO-SHIP, bioGEOTRACES, Malaspina, and Tara Oceans), and citizen science initiatives (Ocean Sampling Day)
-  generate an anvi'o contigs database describing 51 SAR11 isolate genomes (reference genomes)
-  recruite reads from the above listed projects' metagenomes to the SAR11 reference genomes, and profil the recruitment results
-  investigate the patterns in genes across metagenomes recruited to the individual SAR11 reference genomes

Sections in this document will detail all the steps of downloading and processing SAR11 genomes and metagenomes, mapping metagenomic reads onto the SAR11 genomes, as well as analysing and visualising the outcomes.

For the curation of metadata, please consult the [public-marine-omics-metadata](https://github.com/merenlab/public-marine-omics-metadata/tree/main) GitHub repository, where all steps of metadata gathering and standardisation are described in detail.

{:.notice}
If you have any questions, notice an issue, and/or are unable to find an important piece of information here, please feel free to leave a comment down below, send an e-mail to [us](/people/), or get in touch with us through Discord:

{% include _join-anvio-discord.html %}

Information for the projects included in this analysis (information only for depth up to 100 m)
Project | Acronym | Accession | Years covered | # Samples 
| - | - | - | - | - |
[Bermuda Atlantic Time-series Study](https://bats.bios.asu.edu/about/) | BATS | PRJNA385855 | 2003 - 2004, 2009 | 40 
[Bio-GO-SHIP](https://biogoship.org) | BGS | PRJNA656268 | 2011, 2013, 2014, 2016 - 2018 | 444 
[bioGEOTRACES](https://www.nature.com/articles/sdata2018176) | BGT | PRJNA385854 | 2010, 2011 | 323 
[Hawaii Ocean Time-series](http://hahana.soest.hawaii.edu/hot/hot_jgofs.html) | HOT1 \| HOT3 |  PRJNA385855 \| PRJNA352737 | 2003, 2004 \| 2014 - 2017 | 33 \| 230
[Malaspina](https://www.nature.com/articles/s41597-024-02974-1) | MAL | PRJEB52452 | 2011 | 16
OSD
[Tara Oceans](https://fondationtaraocean.org/en/home/) | TARA | PRJEB1787 | 2009 - 2012 | 96


</div>




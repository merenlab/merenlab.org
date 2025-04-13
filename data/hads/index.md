---
layout: page
title: Hawaiʻi Diel Sampling (HADS)
modified: 2024-12-18
excerpt: ""
comments: true
authors: []
---

The purpose of this page is to provide access to the raw data and reproducible data products generated from the Hawaiʻi Diel Sampling (HADS) Project.

We are still in the process of preparing and updating the contents of this page. Please keep an eye on this space for more soon.

{% include IMAGE path="images/hawaii-diel-2021-74.jpg" width=80 caption="Transiting to our offshore sampling site through the Sampan Channel of Kāneʻohe Bay, Oʻahu, Hawaiʻi on the first morning of sample collection." %}

## Motivation

Microbial communities experience environmental fluctuations across timescales ranging from seconds to seasons, and their responses are evident at multiple levels -- from changes in community composition to the physiological reactions of individual cells or from diel cycles to seasonal variations ([Fuhrman et al. 2015](https://www.nature.com/articles/nrmicro3417)).

Time series studies in marine systems have largely focused on resolving changes in microbial community composition from seasonal ([Bunse and Pinhassi 2017](https://doi.org/10.1016/j.tim.2016.12.013); [Giovannoni and Vergin 2012](https://doi.org/10.1126/science.1198078)) to daily timescales ([Needham and Fuhrman 2016](https://doi.org/10.1038/nmicrobiol.2016.5)). While the physiological responses of individual microbial populations offer important insights into their ecology and evolution, such population level responses, especially at short timescaless, are less well understood in complex environments. Responses to short-term fluctuations that occur on timescales that span from seconds to hours are mostly reflected in changes at the level of transcription and translational regulation without any immediate impact on community composition. We generated HADS dataset to contribute an interlinked 'omics resource that lends itself to studies of subtle and population-resolved responses of microbes to environmental change.

HaDS is a collection of metagenomes, metatranscriptomes and metaepitransriptomes generated over a 48 hours period at 3 hour intervals at two sampling sites in Kāneʻohe Bay, Hawaiʻi. The spatiotemporal dynamics of the two surface water sampling stations (so-called HP1 and STO1) are well characterized through the Kāneʻohe Bay Time-series ([Tucker et al. 2021](https://doi.org/10.7717/peerj.12274)), an ongoing monthly time-series sampling program of surface ocean biogeochemistry and microbial communities. Our high-resolution multi-omics approach, combined with concurrent measurements of biogeochemical parameters (chlorophyll, temperature, and nutrient concentration) combined with long-term microbial community and biogeochemistry data at both sampling sites, enables the exploration of microbial population responses to environmental fluctuations and long-term change.

{% include IMAGE path="images/hawaii-diel-2021-80.jpg" width=80 caption="Sample processing next to the docks of Hawaiʻi Institute of Marine Biology (HIMB)." %}

## Data

At both the coastal Kāneʻohe Bay station (HP1) and the adjacent offshore station (STO1), we sampled at 33 time-points across 48 hours, and subsequently produced 59 metatranscriptomes, 65 short-read metagenomes, 8 long-read metagenomes, and 66 metaepitranscriptomes. We also generated four deeply-sequenced short-read metagenomes from samples collected in the late fall and spring prior to HaDS through routine Kāneʻohe Bay Time-series sampling. The following data items give access to RAW sequencing results as well as processed data items that are reproducible with anvi'o v8.0 or later.

* NCBI Project id PRJNA1201851 offers access to all raw data for short-read and long-read metagenomes, as well as metatranscriptomes and metaepitranscriptomes.
* doi:(pending URL from BCO-DMO). Biogeochemical data that covers the sampling period.
* doi:[10.6084/m9.figshare.28784717](https://doi.org/10.6084/m9.figshare.28784717) serves anvi'o {% include ARTIFACT name="contigs-db" %} files for the indivdiual co-assemblies of short-read (SR) as well as long-read (LR) sequencing of metagenomes. Please note that an anvi'o {% include ARTIFACT name="contigs-db" %} includes gene calls, functional annotations, HMM hits, and other information about each contig, and you can always use the program {% include PROGRAM name="anvi-export-contigs" %} to get a FASTA file for sequences. The following screenshot of the {% include PROGRAM name="anvi-display-contigs-stats" %} output gives an idea about the contents of each {% include ARTIFACT name="contigs-db" %} file:

{% include IMAGE path="images/anvi-display-contigs-stats.png" width=80 caption="Sample processing next to the docks of Hawaiʻi Institute of Marine Biology (HIMB)." %}

* doi:[10.6084/m9.figshare.28784762](https://doi.org/10.6084/m9.figshare.28784762) serves FASTA files for metagenome-assembled genomes (MAGs) we have reconstructed from short-read and long-read sequencing of the metagenomes. They are the outputs of quite a preliminary effort, thus secondary attempts to recover genomes from the co-assemblies are most welcome (and very much encouraged).
* doi:[10.6084/m9.figshare.28784765](https://doi.org/10.6084/m9.figshare.28784765). The [EcoPhylo](https://anvio.org/help/main/workflows/ecophylo/) output that describes the phylogeography of ribosomal protein XXX

Please feel free to reach out to us if you have any questions regarding access and/or processing of these datasets.

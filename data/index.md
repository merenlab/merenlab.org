---
layout: page
title: Public Data
modified: 2015-07-30
excerpt: "Publicly available data and other information on publications from our lab"
comments: true
---

This page serves the publicly available data mentioned in our publications, or in publications our group is involved. Please do not hesitate to get in touch if something is missing.

{% include _toc.html %}

## Vineis et al. (2016)

{:.notice}
**The paper itself is here:** [Patient-Specific *Bacteroides* Genome Variants in Pouchitis](http://mbio.asm.org/content/7/6/e01713-16.full).

{:.notice}
**Here is a blog post on it:** [*Bacteroides* Genome Variants, and a reproducible science exercise with anvi'o]({% post_url miscellaneous/2016-11-21-bacteroides-genome-variants %}).

Data for the paper:

* Anvi'o profiles: [https://dx.doi.org/10.6084/m9.figshare.3851364.v1](https://dx.doi.org/10.6084/m9.figshare.3851364.v1){:target="_blank"}

* Primary and supplementary figures: [https://dx.doi.org/10.6084/m9.figshare.3851481.v2](https://dx.doi.org/10.6084/m9.figshare.3851481.v2){:target="_blank"}

* Supplementary tables: [https://dx.doi.org/10.6084/m9.figshare.3851478.v2](https://dx.doi.org/10.6084/m9.figshare.3851478.v2){:target="_blank"}

---

The anvi'o profiles article in this data collection contains 22 items:

[![Abstract]({{site.url}}/images/miscellaneous/2016-11-21-bacteroides-genome-variants/profiles.png)]({{site.url}}/images/miscellaneous/2016-11-21-bacteroides-genome-variants/profiles.png){:.center-img .width-60}

For an example on how to re-analyze these anvi'o profiles, please [click here]({% post_url miscellaneous/2016-11-21-bacteroides-genome-variants %}/#displaying-a-patient-metagenome-from-vineis-et-al){:target="_blank"}.


## Delmont & Eren (2016)

{:.notice}
**The paper itself is here:** [Identifying contamination with advanced visualization and analysis practices: metagenomic approaches for eukaryotic genome assemblies](https://peerj.com/articles/1839/).

* [This link](https://ndownloader.figshare.com/files/3677700) will download the archive file for anvi'o profile for merged datasets (221 Mb). The run script in this archive will automatically start the anvi'o interactive interface, and draw **Figure 1** (compatible with anvi'o v1.2.2, for which a [docker container]({% post_url anvio/2015-08-22-docker-image-for-anvio %}) is available).

* Visit [this adress](http://anvio.org/data/Delmont_et_al-2016-Tardigrade-Summary/) to view the anvi'o summary for genomic selections in Figure 1. Alternatively, you can use [this address](https://ndownloader.figshare.com/files/3677703) to download this static HTML output to view on your own computer.

* Following links give access to the individual genome files: 

  - [Draft genome for *H. dujardini*](http://anvio.org/data/Delmont_et_al-2016-Tardigrade-Summary/bin_by_bin/Tardigrade_draft_genome_01/Tardigrade_draft_genome_01-contigs.fa)
  - [Bacterial draft genome 01](http://anvio.org/data/Delmont_et_al-2016-Tardigrade-Summary/bin_by_bin/Bacterial_draft_genome_01/Bacterial_draft_genome_01-contigs.fa)
  - [Bacterial draft genome 02](http://anvio.org/data/Delmont_et_al-2016-Tardigrade-Summary/bin_by_bin/Bacterial_draft_genome_02/Bacterial_draft_genome_02-contigs.fa)
  - [Bacterial draft genome 03](http://anvio.org/data/Delmont_et_al-2016-Tardigrade-Summary/bin_by_bin/Bacterial_draft_genome_03/Bacterial_draft_genome_03-contigs.fa)

<div style="display: block; height: 0px;">&nbsp;</div>

* [This link](https://ndownloader.figshare.com/files/3677706) will download everything necessary to recreate -an unpolished version of- **Figure 2** appears in the manuscript, including the run script that will run the process automatically (compatible both with v1 and v2 branches of anvi'o).

* Following links give access to media files and supplementary tables:
  - [**Figure 1**](https://ndownloader.figshare.com/files/3677709). *Holistic assessment of the tardigrade genome release from Boothby et al. (2015). Dendrogram in the center organizes scaffolds based on sequence composition and coverage values in data from 11 DNA libraries. Scaffolds larger than 40 kbp were split into sections of 20 kbp for visualization purposes. Splits are displayed in the first inner circle and GC-content (0-71%) in the second circle. In the following 11 layers, each bar represents the portion of scaffolds covered by short reads in a given sample. The next layer shows the same information for RNA-Seq data. Scaffolds harboring genes used by Boothby et al. to support the expended HGT hypothesis is shown in the next layer. Finally, the most outer layer shows our selections of scaffolds as draft genome bins: the curated tardigrade genome (selection number 1), as well as three near-complete bacterial genomes originating from various contamination sources (selection number 2, 3, and 4).*

  - [**Figure 2**](https://ndownloader.figshare.com/files/3677715). *Occurrence of the 139 bacterial single-copy genes reported by Campbell et al. (2013) across scaffold collections. The top two plots display the frequency and distribution of single-copy genes in the raw tardigrade genomic assembly generated by Boothby et al. (2015), and Koutsovoulos et al. (2015), respectively. The bottom two plots display the same information for each of the curated tardigrade genomes. Each bar represents the squared-root normalized number of significant hits per single-copy gene. The same information is visualized as box-plots on the left side of each plot.*

  - [**Supplementary Figire 1**](https://ndownloader.figshare.com/files/3677712). *Visualization and curation of the raw tardigrade genome assembly from Koutsovoulos et al. (2015). In the left panel (curation step I), 24,841 scaffolds that were longer than 1 kbp from the raw assembly were clustered based on sequence composition and coverage values in data from the two Illumina sequencing libraries (the inner dendrogram). Scaffolds longer than 40 kbp were split into sections of 20 kbp for visualization purposes. The second layer shows the GC-content for each scaffold. Next two view layers represent the log-normalized mean coverage values for scaffolds in the two sequencing datasets. Finally, our scaffold selections (tardigrade draft 01 and six bacterial draft genomes) are displayed in the outer layer. In the right panel (curation step II), the 15,839 scaffolds from the tardigrade selection from step I were clustered based on sequence composition only for a more precise curation. Additional scaffold selections (tardigrade draft 02 and two bacterial draft genomes) are displayed in the outer layer.*

  - [**Supplementary Table 1**](https://ndownloader.figshare.com/files/3677718). *Summary of H. dujardini and bacterial genomes identified from the raw assembly results of Boothby et al. (2015) and Koutsovoulos et al. (2015). * Inferred from Boothby et al. (2015) and Koutsovoulos et al. (2015) publications. ** Scores were calculated using bacterial single copy genes from Campbell et al. (2013) and are only used to assess bacterial contamination levels in the eukaryotic assembly results.

  - [**Supplementary Table 2**](https://ndownloader.figshare.com/files/3677721). *Summary of functions identified by RAST in the bacterial draft genome #2 (selection #3 in Fig. 1).*

  - [**Supplementary Table 3**](https://ndownloader.figshare.com/files/3677724). *Summary of HMM hits for each bacterial single-copy gene (collection of 139 from Campbell et al. (2013)) identified in 1) the raw assembly by Boothby et al. (2015), 2) the raw assembly by Koutsovoulos et al. (2015), 3) the curated draft genome of Hypsibius dujardini  from Boothby et al. assembly in this study, and 4) the curated draft genome of H. dujardini from Koutsovoulos et al. (2015).*


*Everything mentioned on this page can be cited using doi [10.6084/m9.figshare.2067057](https://dx.doi.org/10.6084/m9.figshare.2067057).*


## Eren et al. (2015)

{:.notice}
**The paper itself is here:** [Anvi’o: an advanced analysis and visualization platform for ‘omics data](https://peerj.com/articles/1319/).

{:.notice}
The anvi'o profiles here will run with a much earlier version of anvi'o. If you would like to work with them, please checkout your anvi'o codebase to [this commit](https://github.com/meren/anvio/commit/265318856301ab4f7a71911aeaefb524339efbb3). Please don't hesitate to write us if you need assistance.

---

**Daily Infant Gut Samples by Sharon _et al._**. Raw data and anvi'o results for the section on supervised binning and the analysis of the variability in genome bins.

* Visit [this address](http://umkk2268fc06.merenbey.koding.io:8080) to **try the anvi'o interactive interface** on the infant gut data, which provides the basis for **Figure 2**.
* You can [view the summary of the 13 bins](http://anvio.org/data/INFANT-CLC-SUMMARY-SUPERVISED), or you can [download the browsable output](http://anvio.org/data/INFANT-CLC-SUMMARY-SUPERVISED.tar.gz).
* [This Github repository](https://github.com/meren/anvio-methods-paper-analyses) gives access to the code that generates **Figure 3** (see the [relevant](https://github.com/meren/anvio-methods-paper-analyses/tree/master/SHARON_et_al/VARIABILITY_REPORTS) directory).
* You can download the output of `anvi-merge` (the merged profile db, and the annotation db) for the infant gut metagenomes [from here](http://dx.doi.org/10.6084/m9.figshare.1499236).

---

**Pensacola Beach Samples by Overholt _et al._ and Rodriguez-R _et al._**. Raw data and anvi'o results for the section on linking cultivar genomes with metagenomes.

* While [this address](http://anvio.org/data/OVERHOLT-CULTIVARS-SUMMARY) gives access to the anvi'o summary of the ten cultivar genomes ([download](http://anvio.org/data/OVERHOLT-CULTIVARS-SUMMARY.tar.gz)), [this one](http://anvio.org/data/RODRIGUEZ-R-MG-SUMMARY) serves the 56 metagenomic bins ([download](http://anvio.org/data/RODRIGUEZ-R-MG-SUMMARY.tar.gz)) shown in **Figure 4**.
* You can download the output of `anvi-merge` for the mapping of metagenomes to Overholt cultivars [from here](http://dx.doi.org/10.6084/m9.figshare.1499234), and the output for `anvi-merge` for metagenomic bins is available [here](http://dx.doi.org/10.6084/m9.figshare.1499248).

---

**Gulf of Mexico Samples by Mason _et al._, and Yergeau _et al._**. Results for the section on linking metagenomes, metatranscriptomes, and single-cell genomes.

* You can [view the summary of the two bins](http://anvio.org/data/MASON-SAGs-SUMMARY-SUPERVISED) identified in the assembly of the single-cell genomes ([download](http://anvio.org/data/MASON-SAGs-SUMMARY-SUPERVISED.tar.gz)) (**Figure 5 panel A**).
* You can [view the summary of the three bins](http://anvio.org/data/MASON-YERGEAU-MG-SUMMARY-SUPERVISED) identified in the metagenomic assembly ([download](http://anvio.org/data/MASON-YERGEAU-MG-SUMMARY-SUPERVISED.tar.gz)) (**Figure 5 panel B**).
* [This Github repository](https://github.com/meren/anvio-methods-paper-analyses) also gives access to the code that generates **Figure 5 panel C** (see the [relevant](https://github.com/meren/anvio-methods-paper-analyses/tree/master/MASON_et_al/SCATTER_PLOTS) directory).
* You can download the output of `anvi-merge` for the mapping of all samples against the assembly of SAGs [from here](http://dx.doi.org/10.6084/m9.figshare.1499235), and the output for `anvi-merge` for mapping to metagenomic contigs is available [here](http://dx.doi.org/10.6084/m9.figshare.1499246).

---

**Media and Supplementary files**.

* Tables: [Additional file 1](http://dx.doi.org/10.6084/m9.figshare.1499237), [Additional file 2](http://dx.doi.org/10.6084/m9.figshare.1499238), [Additional file 3](http://dx.doi.org/10.6084/m9.figshare.1499239).
* Figures: [Figure 1](http://dx.doi.org/10.6084/m9.figshare.1499240), [Figure 2](http://dx.doi.org/10.6084/m9.figshare.1499241), [Figure 3](http://dx.doi.org/10.6084/m9.figshare.1499242), [Figure 4](http://dx.doi.org/10.6084/m9.figshare.1499243), [Figure 5](http://dx.doi.org/10.6084/m9.figshare.1499244), [Figure S1](http://dx.doi.org/10.6084/m9.figshare.1499245), [Figure S2](http://dx.doi.org/10.6084/m9.figshare.1499247).


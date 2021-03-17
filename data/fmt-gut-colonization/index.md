---
layout: page
title: The Fecal Microbiota Transplantation study
modified: 2021-03-16
excerpt: "by Andrea R. Watson et al 2021"
comments: true
image:
  featurerelative: images/header.png
  display: true
---

{% include _toc.html %}


<div class="extra-info" markdown="1">

<span class="extra-info-header">Summary</span>

**The purpose of this page** is to discuss key findings of the study "**[Adaptive ecological processes and metabolic independence drive microbial colonization and resilience in the human gut](https://doi.org/10.1101/2021.03.02.433653)**" by Watson et al.

On this page you will find data and ad hoc scripts for important analyses, results of which discussed in the paper.

</div>

{:.notice}
If you have any questions and/or if you are unable to find an important piece of information here, please feel free to leave a comment down below, send an e-mail to [us]({{ site.url }}/people/), or get in touch with us through Slack:

{% include _join-anvio-slack.html %}


## Study description

Our study involves the investigations of human gut metagenomes to understand the ecology of microbes before and after fecal microbiota transplantation (FMT) experiments. To do that, we 

* Reconstruct metagenome-assembled genomes (MAGs) from FMT donor gut metagenomes,
* Using these genomes investigate the distribution of donor populations in donor and recipient gut metagenomes before and after FMT,
* Investigate the prevalence of donor populations in global gut metagenomes,
* Define a set of quantitative and stringent criteria to determine 'colonization' events,
* Identify 'high-fitness' and 'low-fitness' donor populations based on their colonization abilities and prevalence in global gut metagenomes,
* Determine metabolic pathways that are enriched in either of these high-fitness or low-fitness groups,
* Study microbial genomes reconstructed from healthy individuals and individuals with IBD to show that high-fitness genomes are differentially enriched in IBD.

Most of these steps are detailed in [our manuscript](https://doi.org/10.1101/2021.03.02.433653) and on this page you will find details of intermediate steps and anvi'o reusable data currencies to reproduce our main findings.

## Reconstructing donor genomes

*To be described.*

{:.notice}
[doi:](){:target="_blank"} serves anvi'o {% include ARTIFACT name="contigs-db" text="contigs databases" %} for donor assemblies, and merged anvi'o {% include ARTIFACT name="profile-db" text="profile databases" %} that show the distribution of genomes across donor and recipient metagenomes. Use XXX for {% include ARTIFACT name="collection" %} name.


## Investigating metabolic competence among microbial genomes reconstructed from healthy individuals and individuals with IBD

{:.warning}
[doi:10.6084/m9.figshare.14225840](https://doi.org/10.6084/m9.figshare.14225840) serves each genome used in the analysis below.

This section will describe how to go from FASTA files for a given set of genomes to the Figure 4 in our study:

[![Figure 04](images/Figure_04.png)](images/Figure_04.png){:.center-img .width-50}

Our analyses in the previous chapters of our study showed that while the healthy donor environment could support both high-fitness and low-fitness populations, challenging microbes to colonize a new environment or to withstand massive ecosystem perturbation during FMT selects for high-fitness populations, suggesting that metabolic competence is a more critical determinant of fitness during stress than during homeostasis. Based on these observations, we hypothesized that,

* A gut environment in homeostasis will support a range of microbial populations with a wide spectrum of metabolic competency, and
* A gut environment under stress will select for high metabolic competency in microbial populations.

### Publicly available datasets

To test these hypotheses, we proceeded to compared genomes reconstructed from a cohort of healthy individuals by [Pasolli et al. (2019)](https://doi.org/10.1016/j.cell.2019.01.001) to genomes reconstructed from individuals who were diagnosed with inflammatory bowel disease (IBD). Our IBD dataset was composed of two cohorts: a set of patients with pouchitis from [Vineis et al (2016)](https://doi.org/10.1128/mBio.01713-16), and a set of pediatric Crohn's disease patients described in [Quince et al. (2015)](https://doi.org/10.1038/ajg.2015.357).

### FASTA files to anvi'o contigs databases

At this stage we have multiple FASTA files for each individual in three cohorts (healthy, pouchitis, and Crohn's) for downstream analyses. We first converted each FASTA file into an anvi'o {% include ARTIFACT name="contigs-db" %} using the program {% include PROGRAM name="anvi-gen-contigs-database" %} and then processed each resulting contigs-db using,

* The program {% include PROGRAM name="anvi-run-hmms" %} to run hidden Markov models for the identification of single-copy core genes,
* The program {% include PROGRAM name="anvi-run-scg-taxonomy" %} to identify and store taxonomic information,
* and the program {% include PROGRAM name="anvi-run-kegg-kofams" %} to identify KEGG orthologs in each of our genomes for metabolic reconstruction.

The resulting contigs databases are publicly available. If you would like to reproduce the steps below, you should first download the following file [943 Mb] to your work directory,

``` bash
curl -L https://ndownloader.figshare.com/files/26842652 \
     -o WATSON_ET_AL_CONTIGS_DBS_FOR_METABOLIC_COMPARISON.tar.gz
```

And unpack it [2.7 Gb]:

``` bash
tar -zxvf WATSON_ET_AL_CONTIGS_DBS_FOR_METABOLIC_COMPARISON.tar.gz
```

After having generated a contigs database for each genome, we used all the information about these genomes that were available to us to prepared an {% include ARTIFACT name="external-genomes" %} file. This is simply a TAB-delimited text file that describes the location of each contigs database for each genome on the disk, along with optional metadata columns to describe their properties. Here are a few lines from this file:

|name|individual|cohort|from|accession ID|genome_completion|genome_redundancy|num_splits|total length|t_phylum|t_class|t_order|t_family|t_genus|t_species|contigs_db_path|
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
|POU_P204_MAG_00001|POU_P204|POUCHITIS|Vineis et al. (doi:10.1128/mBio.01713-16)|phs000262.v3.p2 (dbGaP)|94.37|0.00|104|2081795|Firmicutes|Clostridia|Oscillospirales|Acutalibacteraceae|UMGS1071|None|WATSON_ET_AL_CONTIGS_DBS_FOR_METABOLIC_COMPARISON/POU_P204_MAG_00001.db|
|POU_P204_MAG_00002|POU_P204|POUCHITIS|Vineis et al. (doi:10.1128/mBio.01713-16)|phs000262.v3.p2 (dbGaP)|100.00|0.00|129|2564734|Firmicutes|Clostridia|Lachnospirales|Anaerotignaceae|Anaerotignum|Anaerotignum sp001304995|WATSON_ET_AL_CONTIGS_DBS_FOR_METABOLIC_COMPARISON/POU_P204_MAG_00002.db|
|POU_P204_MAG_00003|POU_P204|POUCHITIS|Vineis et al. (doi:10.1128/mBio.01713-16)|phs000262.v3.p2 (dbGaP)|98.59|1.41|202|4105023|Bacteroidota|Bacteroidia|Bacteroidales|Bacteroidaceae|Bacteroides|Bacteroides fragilis|WATSON_ET_AL_CONTIGS_DBS_FOR_METABOLIC_COMPARISON/POU_P204_MAG_00003.db|
|POU_P204_MAG_00004|POU_P204|POUCHITIS|Vineis et al. (doi:10.1128/mBio.01713-16)|phs000262.v3.p2 (dbGaP)|98.59|5.63|263|5146152|Bacteroidota|Bacteroidia|Bacteroidales|Bacteroidaceae|Bacteroides|None|WATSON_ET_AL_CONTIGS_DBS_FOR_METABOLIC_COMPARISON/POU_P204_MAG_00004.db|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|
|HMP_1061_SRS011061_bin.13|HMP_1061|HEALTHY|Passoli et al. (doi:10.1016/j.cell.2019.01.001)|SRS011061|76.06|40.85|233|1846034|Bacteroidota|Bacteroidia|Bacteroidales|Rikenellaceae|Alistipes|Alistipes putredinis|WATSON_ET_AL_CONTIGS_DBS_FOR_METABOLIC_COMPARISON/HMP_2012__SRS011061__bin.13.db|
|HMP_1061_SRS011061_bin.14|HMP_1061|HEALTHY|Passoli et al. (doi:10.1016/j.cell.2019.01.001)|SRS011061|98.59|0.00|219|2831497|Firmicutes|Clostridia|Oscillospirales|Ruminococcaceae|CAG-353|CAG-353 sp900066885|WATSON_ET_AL_CONTIGS_DBS_FOR_METABOLIC_COMPARISON/HMP_2012__SRS011061__bin.14.db|
|HMP_1061_SRS011061_bin.20|HMP_1061|HEALTHY|Passoli et al. (doi:10.1016/j.cell.2019.01.001)|SRS011061|85.92|0.00|300|1255151|Firmicutes|Clostridia|TANB77|CAG-508|CAG-273|CAG-273 sp000435755|WATSON_ET_AL_CONTIGS_DBS_FOR_METABOLIC_COMPARISON/HMP_2012__SRS011061__bin.20.db|
|HMP_1061_SRS011061_bin.23|HMP_1061|HEALTHY|Passoli et al. (doi:10.1016/j.cell.2019.01.001)|SRS011061|77.46|4.23|391|1510418|Firmicutes|Clostridia|Oscillospirales|Acutalibacteraceae|UMGS172|UMGS172 sp900539855|WATSON_ET_AL_CONTIGS_DBS_FOR_METABOLIC_COMPARISON/HMP_2012__SRS011061__bin.23.db|
|HMP_1061_SRS011061_bin.26|HMP_1061|HEALTHY|Passoli et al. (doi:10.1016/j.cell.2019.01.001)|SRS011061|100.00|2.82|232|1903353|Firmicutes|Negativicutes|Acidaminococcales|Acidaminococcaceae|Succiniclasticum|Succiniclasticum sp900544275|WATSON_ET_AL_CONTIGS_DBS_FOR_METABOLIC_COMPARISON/HMP_2012__SRS011061__bin.26.db|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|

You can also download a copy of it on your disk and explore it the following way:

``` bash
curl -L https://merenlab.org/data/fmt-gut-colonization/files/external-genomes.txt \
     -o external-genomes.txt
```

The next step is to estimate the completion of each metabolic module described in KEGG database for each genome using the program {% include PROGRAM name="anvi-estimate-metabolism" %}:

### Computing metabolic completion

``` bash
anvi-estimate-metabolism -e external-genomes.txt \
                         --matrix-format \
                         -O metabolism
```

On a laptop computer this command will take up-to 25 minutes to characterize all metabolic modules in 604 genomes, and will generate the following files:

* `metabolism-completeness-MATRIX.txt`
* `metabolism-ko_hits-MATRIX.txt`
* `metabolism-presence-MATRIX.txt`

The output file we will use for downstream analyses is `metabolism-completeness-MATRIX.txt` (and if you don't want to wait for the {% include PROGRAM name="anvi-estimate-metabolism" %}, you can download this file into your work directory from [here](files/metabolism-completeness-MATRIX.txt)). In this particular file each column represents a genome, each row represents a KEGG module, and each data point is the estimated completeness of a given KEGG module in a given genome. Here is a few lines and columns from this file:

|module|POU_P204_MAG_00001|POU_P204_MAG_00002|POU_P204_MAG_00003|POU_P204_MAG_00004|POU_P204_MAG_00005|(...)|
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
|M00001|0.8|0.9|1.0|1.0|0.9|(...)|
|M00002|0.66|1.0|1.0|1.0|1.0|(...)|
|M00003|0.87|0.87|1.0|1.0|0.75|(...)|
|M00004|0.57|0.71|1.0|1.0|1.0|(...)|
|M00005|1.0|1.0|1.0|1.0|1.0|(...)|
|M00006|0.0|0.0|1.0|1.0|1.0|(...)|
|M00007|0.75|1.0|1.0|1.0|1.0|(...)|
|M00008|0.0|0.0|0.75|0.75|0.5|(...)|
|M00009|0.25|0.5|0.62|0.75|0.56|(...)|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|

### Visualizing data in R

The rest of our analysis will take place in R. You can download the entire R script we used to generate relevant figure, and run it in your work directory get all the figures in a single step the following way:

``` bash
# download the script
curl -L https://merenlab.org/data/fmt-gut-colonization/files/plot-metabolic-competency.R \
     -o plot-metabolic-competency.R

# turn the executable bit on
chmod +x plot-metabolic-competency.R

# run it, and investigate the resulting PDF files:
./plot-metabolic-competency.R
```

The following steps highlight and describe some of the relevant steps in the script. As a reminder, we will need the following libraries to be initialized for everything downstream to run:

``` R
#!/usr/bin/env Rscript
library(ggplot2)
library(MASS)
library(ggridges)
library(gridExtra)
library(reshape2)
library(reshape)
```

We first read in the metabolic completion and external genomes files:

```
module_completion <- read.table(file='metabolism-completeness-MATRIX.txt', header = TRUE, sep = "\t")
external_genomes <- read.table(file='external-genomes.txt', header = TRUE, sep = "\t")
```

A quick check to make sure the number of genomes per group makes sense based on what we know about our groups:

```
   CROHNS  FMT_HIGH_FITNESS  FMT_LOW_FITNESS   HEALTHY    POUCHITIS       256                20               20       264           44
```

We then turn the matrix formatted file into a data frame, and add some relevant information per genome using the external genomes file:

``` R
# turn the boring matrix format into a data frame
df <- melt(module_completion)

# set some meaningful column names
colnames(df) <- c('module', 'genome', 'completion')

# use the external genomes file to associate each genome with a 'group',
# and an individual:
df$group <- external_genomes$cohort[match(df$genome, external_genomes$name)]
df$individual <- external_genomes$individual[match(df$genome, external_genomes$name)]
```

Our matrix contains 303 modules. But our statistical enrichment analysis of the high-fitness versus low-fitness genomes had revealed a set of 33 modules. Thus, we want to investigate only those modules across these new genomes to investigate whether there is a difference in genomes across cohorts:

``` R
df <- df[df$module %in% modules_of_interest, ]
```

Finally, some boring steps of defining explicit orders for x-axes in the visualizations that will follow:

``` R
# some boring steps of defining explicit orders for x-axes in boxplots
HEALTHY_subset <- aggregate(df[df$group == "HEALTHY", 3], list(df[df$group == "HEALTHY", ]$individual), median)
HEALTHY_order <- HEALTHY_subset[order(-HEALTHY_subset$x),]$Group.1
POUCHITIS_subset <- aggregate(df[df$group == "POUCHITIS", 3], list( df[df$group == "POUCHITIS", ]$individual), median)
POUCHITIS_order <- POUCHITIS_subset[order(-POUCHITIS_subset$x),]$Group.1
CROHNS_subset <- aggregate(df[df$group == "CROHNS", 3], list( df[df$group == "CROHNS", ]$individual), median)
CROHNS_order <- CROHNS_subset[order(-CROHNS_subset$x), ]$Group.1
individuals_order <- c(c('FMT_HIGH_FITNESS', 'FMT_LOW_FITNESS'), HEALTHY_order, POUCHITIS_order, CROHNS_order)

# set explicit group orders, and assign some group colors
groups_order <- c("FMT_HIGH_FITNESS", "FMT_LOW_FITNESS", "HEALTHY", "POUCHITIS", "CROHNS")
group_colors <- c("#ec5f1c", "#034f84", "#feb236", "#86af49", "#ff0202")
```

After these preparations, we are ready to run the following section,

``` R
# plot the boxplots
pdf(file = "boxplots.pdf",  width = 13, height = 5)
plot_individuals <- ggplot(data=df, aes(x=individual, y=completion, group=individual)) +
  geom_boxplot(aes(fill=group), alpha=0.35, outlier.shape = NA, color='#808080') +
  geom_jitter(colour='#222222', width = 0.20, height = 0.02, size=0.1, alpha=0.5) +
  theme_bw() +
  theme(legend.position="bottom", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  xlab("Individuals") +
  scale_x_discrete(limits = individuals_order) +
  scale_fill_manual(values = group_colors)

plot_groups <- ggplot(data=df, aes(x=group, y=completion, group=group)) +
  geom_boxplot(aes(fill=group), alpha=0.35, outlier.shape = NA, color=NA) +
  geom_violin(fill="#505050", alpha=0.35, width=1.3, colour = '#505050') +
  geom_jitter(colour='#222222', width = 0.3, height = 0.02, size=0.1, alpha=0.05) +
  theme_bw() +
  theme(legend.position="bottom") +
  ylab("Metabolic module completion") +
  xlab("Groups") +
  scale_x_discrete(limits = groups_order) +
  scale_fill_manual(values = group_colors)

grid.arrange(plot_groups, plot_individuals, ncol=2, widths=c(1, 3))
dev.off()
```

which generates a PDF file in our work directory that looks like this:

[![Figure 04 top panel](images/boxplots.png)](images/boxplots.png){:.center-img .width-100}

To visualize the ridge-line plots, we first selected a set of modules that differed significantly between the genomes from the healthy cohort and from those who were diagnosed with IBD to show some examples:


``` R
# modules_to_display <- c("M00924","M00122", "M00023", "M00028", "M00570", "M00082", "M00844", "M00015",  "M00526", "M00022") 
dfx <- df[df$module %in% modules_to_display, ]
dfx$module = factor(dfx$module, levels=modules_to_display)
```

Then, running the following section,

``` R
pdf(file = "ridges-for-metabolisms.pdf",  width = 13, height = 5)
ggplot(data=dfx, aes(x = completion, y = individual, fill = group)) +
  geom_density_ridges(scale = 10, size = 0.25, rel_min_height = 0.01, alpha=0.35, colour="#222222") +
  theme_ridges() +
  scale_y_discrete(limits = individuals_orderx) +
  scale_fill_manual(values = group_colors) +
  scale_x_continuous(breaks=c(0, 1)) +
  coord_cartesian(xlim = c(-0.5, 1.5)) +
  theme(legend.position="bottom", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  facet_wrap(. ~ module, nrow = 2)
dev.off()
```

which generates a PDF file in our work directory that looks like this:

[![Figure 04 bottom panel](images/ridges-for-metabolism.png)](images/ridges-for-metabolism.png){:.center-img .width-100}

Finally we run the following test to estimate the statistical significance of differences to report in our study,

``` R
printf <- function(...) invisible(print(sprintf(...)))
printf("DIFFERENCES BETWEEN HEALTHY vs CROHNS + POUCHITIS for RIDGELINE PLOTS")
for(module in modules_to_display){
    w <- wilcox.test(dfx[dfx$group == "HEALTHY" & dfx$module == module, ]$completion, dfx[dfx$group %in% c("CROHNS", "POUCHITIS") & dfx$module == module, ]$completion)
    printf("%s: %f", module, w$p.value)
}
```

And combined both panels in [Inkscape](https://inkscape.org/) to finalize the figure for publication.

<div style="display: block; height: 200px;">&nbsp;</div>


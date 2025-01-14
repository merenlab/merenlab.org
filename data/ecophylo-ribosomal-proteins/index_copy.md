---
layout: page
title: Phylogeography of ribosomal proteins with EcoPhylo
modified: 2021-10-21
excerpt: "by Schechter et al, 2025"
comments: true
authors: [matt]
---

- [ ] FIXME: change title and link to preprint
- [ ] FIXME: Short notice that most of this was run with clusterize on an HPC
- [ ] FIXME: remove clusterize from all the commands

**The purpose of this page** is to provide access to reproducible data products that underlie our key findings in the study "**[Phylogeography of ribosomal proteins reveals missing genomes from genome collections]()**" by [Matt Schechter]() et al.

{:.notice}
If you have any questions and/or if you are unable to find an important piece of information here, please feel free to leave a comment down below, send an e-mail to [us]({{ site.url }}/people/), or get in touch with us through Discord:

{% include _join-anvio-discord.html %}

## Set up environment

Download and open up this

``` bash
cd 2024_Schechter_EcoPhylo

open 2024_Schechter_EcoPhylo.Rproj
```

Open the file `REPRODUCIBLE_WORKFLOW.md` to run everything from the
notebook itself. This only things that will be missing will be merenlab
website markdown modifications.

**Load packages**

These are the `r` packages you will need to load and install before
running the analysis

``` r
packages <- c("tidyverse", "ggpubr", "fs", "ape", "treeio", "glue", "plotly", "readxl", "ggridges")
suppressWarnings(suppressMessages(lapply(packages, library, character.only = TRUE)))

# set path for rendering notebook
knitr::opts_knit$set(here::here()) 

# load custom functions
source(here::here("SCRIPTS", "plotting.R"))
source(here::here("SCRIPTS", "utils.R"))

# set path for running notebook in Rstudio
setwd(here::here())
```

## Reproducible / Reusable Data Products

Here we will put a link to the repo and any of the figshare data

- [ ] Human oral cavity `contigs-db` and `profile-db`
- [ ] Global surface ocean `contigs-db` and `profile-db`
- [ ] Hadza tribe gut `contigs-db` and `profile-db`

## Distribution HMM alignment coverage and SCG detection across GTDB

> The majority of commands using
> [anvi-run-workflow](https://anvio.org/help/main/programs/anvi-run-workflow/)
> were run of high performance compute clusters leveraging the tool
> powerful tool [clusterize](https://github.com/ekiefl/clusterize).
> Without access to compute nodes, these commands would take a VERY long
> time.

To identify a threshold to remove spurious hits, we examined the
distribution of `hmm-hit` alignment coverage across genomes across
[GTDB](https://gtdb.ecogenomic.org/). Here are the steps to reproduce
this:

**Step 1.** Run ecophylo over GTDB Bacteria and Archaea collections with
the `Bacteria_71` and `Archaea_76` HMM collections respectively with NO
HMM alignment coverage cutoff.

This allowed us to record all possible alignment coverage values and
examine the distribution.

``` bash
# Bacteria
# get default config
anvi-run-workflow -w ecophylo --get-default-config ecophylo_config_GTDB_reps.json

# Make log dir
mkdir -p ECOPHYLO_WORKFLOW_GTDB_REPS/00_LOGS

# Run workflow
anvi-run-workflow -w ecophylo -c ecophylo_config_GTDB_reps.json --skip-dry-run --additional-params --cluster \"clusterize -j={rule} -o={log} -n={threads} --exclude \'\' -x \" --jobs=120 --resources nodes=120 --latency-wait 100 --keep-going --rerun-incomplete

# Archaea
# cp config from bacteria and modify
cp ecophylo_config_GTDB_reps.json ecophylo_config_GTDB_reps_archaea.json

# Make log dir
mkdir -p ECOPHYLO_WORKFLOW_GTDB_REPS_ARCHAEA/00_LOGS

# Run workflow
anvi-run-workflow -w ecophylo -c ecophylo_config_GTDB_reps_archaea.json --skip-dry-run --additional-params --cluster \"clusterize -j={rule} -o={log} -n={threads} --exclude \'\' -x \" --jobs=120 --resources nodes=120 --latency-wait 100 --keep-going --rerun-incomplete
```

**Step 2.** Extract `hmmsearch`
[domtblout](http://eddylab.org/software/hmmer/Userguide.pdf) files from
EcoPhylo dir structure.

Here we will explore the distribution of HHM alignment coverages
(`Bacteria_71`, `Archaea_76`) across a subset of bacterial GTDB genomes.

``` bash
cd /project2/meren/PROJECTS/ECOPHYLO/Surface_Ocean

# Bacteria
for file in `ls ECOPHYLO_WORKFLOW_GTDB_REPS/01_REFERENCE_PROTEIN_DATA/*/Ribosomal_L16-dom-hmmsearch/hmm.domtable`;
do
    fname=$(dirname "${file}" | sed 's|/Ribosomal_L16-dom-hmmsearch||' | sed 's|ECOPHYLO_WORKFLOW_GTDB_REPS/01_REFERENCE_PROTEIN_DATA/||')
    echo -e "${fname}"
    sed "s/^/"${fname}"\t&/g" "${file}" >> hmm.domtable.GTDB.txt;
done

# Archaea
for file in `ls ECOPHYLO_WORKFLOW_GTDB_REPS_ARCHAEA/01_REFERENCE_PROTEIN_DATA/*/Ribosomal_L6-dom-hmmsearch/hmm.domtable`;
do
    fname=$(dirname "${file}" | sed 's|/Ribosomal_L6-dom-hmmsearch||' | sed 's|ECOPHYLO_WORKFLOW_GTDB_REPS_ARCHAEA/01_REFERENCE_PROTEIN_DATA/||')
    echo -e "${fname}"
    sed "s/^/"${fname}"\t&/g" "${file}" >> hmm.domtable.GTDB.archaea.txt;
done
```

**Step 3.** Clean the domtblout data

``` r
# Read in domtlout files
GTDB_domtable <- read_table("DATA/hmm.domtable.GTDB.txt", col_names = FALSE)
GTDB_domtable_archaea <- read_table("DATA/hmm.domtable.GTDB.archaea.txt", col_names = FALSE)

fix_domtblout_data <- function(X) {
   X %>%
  rename("genome" = X1,
         "gene_callers_id" = X2,
         "target_accession" = X3,
         "gene_length" = X4,
         "hmm_name" = X5,
         "hmm_id" = X6,
         "hmm_length" = X7,
         "evalue" = X8,
         "bitscore" = X9,
         "bias" = X10,
         "match_num" = X11,
         "num_matches" = X12,
         "dom_c_evalue" = X13,
         "dom_i_evalue" = X14,
         "dom_bitscore" = X15,
         "dom_bias" = X16,
         "hmm_start" = X17,
         "hmm_stop" = X18,
         "gene_start" = X19,
         "gene_stop" = X20,
         "env_to" = X21,
         "env_from" = X22,
         "mean_post_prob" = X23,
         "description" = X24) %>%
  mutate(model_coverage = (hmm_stop - hmm_start)/hmm_length) %>%
  mutate(gene_coverage = (gene_stop - gene_start)/gene_length) %>%
  select(genome, gene_callers_id, hmm_name, model_coverage, gene_coverage) %>%
  arrange(desc(model_coverage))  
}

GTDB_domtable_ed <- fix_domtblout_data(GTDB_domtable)
GTDB_domtable_archaea_ed <- fix_domtblout_data(GTDB_domtable_archaea)

# Save for later
GTDB_domtable_ed %>% write_tsv("DATA/GTDB_domtable_ed.tsv")
GTDB_domtable_archaea_ed %>% write_tsv("DATA/GTDB_domtable_archaea_ed.tsv")
```

**Step 4.** Plot the distribution of model and gene alignment coverages

Here we will plot the distribution of HMM alignment coverages to all the
ORF hits across the genome datasets

``` r
GTDB_domtable_ed <- read_tsv("DATA/GTDB_domtable_ed.tsv")
GTDB_domtable_archaea_ed <- read_tsv("DATA/GTDB_domtable_archaea_ed.tsv")

plot_boxplot_HMM_alignment_coverage_distribution <- function(X, HMM_source) {
  
  ###
  # X <- GTDB_domtable_archaea_ed
  # HMM_source <- "Archaea_76"
  ###
  plot_model_cov <- X %>%
  group_by(hmm_name) %>%
  mutate(avg = mean(model_coverage),
         median = median(model_coverage)) %>%
  ggplot(aes(x=model_coverage, y=fct_reorder(hmm_name, median))) +
  geom_boxplot(outlier.size = 0.5,
               outlier.shape = 1,
               outlier.alpha = 0.5) +
  geom_vline(xintercept = 0.8, color = "blue") +
  geom_vline(xintercept = 0.95, color = "red") +
  theme_light() +
  theme(legend.position = "none",
        axis.ticks.x=element_blank(),
        legend.title = element_blank(),
        axis.text=element_text(size = 6),
        text=element_text(family="Helvetica Neue Condensed Bold"),
        plot.title = element_text(hjust = 0.5)
        ) +
  xlab("HMM model alignment coverage") +
  ylab(str_c(HMM_source, "model ordered by mean coverage", sep = " ")) +
  ggtitle(str_c("hmmsearch model coverage for", HMM_source, sep = " "))
  
  plot_gene_cov <- X %>%
  group_by(hmm_name) %>%
  mutate(avg = mean(gene_coverage),
         median = median(gene_coverage)) %>%
  ggplot(aes(x=gene_coverage, y=fct_reorder(hmm_name, median))) +
  geom_boxplot(outlier.size = 0.5,
               outlier.shape = 1,
               outlier.alpha = 0.5) +
  geom_vline(xintercept = 0.8, color = "blue") +
  geom_vline(xintercept = 0.95, color = "red") +
  theme_light() +
  theme(legend.position = "none",
        axis.ticks.x=element_blank(),
        legend.title = element_blank(),
        axis.text=element_text(size = 6),
        text=element_text(family="Helvetica Neue Condensed Bold"),
        plot.title = element_text(hjust = 0.5)
        ) +
  xlab("gene alignment coverage") +
  ylab(str_c(HMM_source, "gene ordered by mean coverage", sep = " ")) +
  ggtitle(str_c("hmmsearch gene-coverage for", HMM_source, sep = " "))
  
  ridges_model_coverage <- X %>%
    group_by(hmm_name) %>%
    mutate(avg = mean(model_coverage),
           median = median(model_coverage)) %>%
    ggplot(aes(x=model_coverage, y=fct_reorder(hmm_name, median))) +
      geom_density_ridges(rel_min_height = 0.01,
                          # jittered_points = TRUE,
                          # position = position_points_jitter(width = 0.05, height = 0),
                          # point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7
                          ) +
      geom_vline(xintercept = 0.8, color = "blue") +
    geom_vline(xintercept = 0.95, color = "red") +
    theme_light() +
    theme(legend.position = "none",
          axis.ticks.x=element_blank(),
          legend.title = element_blank(),
          axis.text=element_text(size = 6),
          text=element_text(family="Helvetica Neue Condensed Bold"),
          plot.title = element_text(hjust = 0.5)
          ) +
    xlab("HMM model alignment coverage") +
    ylab(str_c(HMM_source, "model ordered by mean coverage", sep = " ")) +
    ggtitle(str_c("hmmsearch model-coverage for", HMM_source, sep = " "))
  
  ridges_gene_coverage <- X %>%
    group_by(hmm_name) %>%
    mutate(avg = mean(gene_coverage),
           median = median(gene_coverage)) %>%
    ggplot(aes(x=gene_coverage, y=fct_reorder(hmm_name, median))) +
      geom_density_ridges(rel_min_height = 0.01) +
      geom_vline(xintercept = 0.8, color = "blue") +
    geom_vline(xintercept = 0.95, color = "red") +
    theme_light() +
    theme(legend.position = "none",
          axis.ticks.x=element_blank(),
          legend.title = element_blank(),
          axis.text=element_text(size = 6),
          text=element_text(family="Helvetica Neue Condensed Bold"),
          plot.title = element_text(hjust = 0.5)
          ) +
    xlab("gene alignment coverage") +
    ylab(str_c(HMM_source, "gene ordered by mean coverage", sep = " ")) +
    ggtitle(str_c("hmmsearch gene-coverage for", HMM_source, sep = " "))
  
  ridges_model_coverage <- X %>%
    group_by(hmm_name) %>%
    mutate(avg = mean(model_coverage),
           median = median(model_coverage)) %>%
    ggplot(aes(x=model_coverage, y=fct_reorder(hmm_name, median))) +
      geom_density_ridges(rel_min_height = 0.01) +
      geom_vline(xintercept = 0.8, color = "blue") +
    geom_vline(xintercept = 0.95, color = "red") +
    theme_light() +
    theme(legend.position = "none",
          axis.ticks.x=element_blank(),
          legend.title = element_blank(),
          axis.text=element_text(size = 6),
          text=element_text(family="Helvetica Neue Condensed Bold"),
          plot.title = element_text(hjust = 0.5)
          ) +
    xlab("HMM model alignment coverage") +
    ylab(str_c(HMM_source, "model ordered by mean coverage", sep = " ")) +
    ggtitle(str_c("hmmsearch model alignment coverage for", HMM_source, sep = " "))
  
  list(plot_model_cov = plot_model_cov,
       plot_gene_cov = plot_gene_cov,
       ridges_model_coverage = ridges_model_coverage, 
       ridges_gene_coverage = ridges_gene_coverage)
}

plot_model_cov_bacteria_71 <- plot_boxplot_HMM_alignment_coverage_distribution(GTDB_domtable_ed, HMM_source = "Bacteria_71")
plot_model_cov_archaea_76 <- plot_boxplot_HMM_alignment_coverage_distribution(GTDB_domtable_archaea_ed, HMM_source = "Archaea_76")

ggarrange(plot_model_cov_archaea_76$plot_model_cov, plot_model_cov_bacteria_71$plot_model_cov, plot_model_cov_archaea_76$plot_gene_cov, plot_model_cov_bacteria_71$plot_gene_cov, ncol = 2, nrow = 2)
```

![](REPRODUCIBLE_WORKFLOW_files/figure-commonmark/unnamed-chunk-6-1.png)

{% include IMAGE path="images/Supplemental_1_reformatted.png" width=50 caption="Distribution of HMM alignment coverage of Bacteria_71 and Archaea_76 SCG HMM collections searched against GTDB RefSeq Archaea and Bacteria representative genomes." %}

## Explore SCG HMM copy number across GTDB r95 RefSeq

To confirm and identify ribosomal proteins that come in single-copy per
genomes, we searched `Bacteria_71` and `Archaea_76` HMM collections
across GTDB with 80% model alignment coverage (threshold identified in
the previous analysis).

**Step 1.** Re-run EcoPhylo with 80% HMM alignment coverage cut-off

``` bash
# Bacteria
# get default config
anvi-run-workflow -w ecophylo --get-default-config ecophylo_config_GTDB_reps.json

# Make log dir
mkdir -p ECOPHYLO_WORKFLOW_GTDB_REPS_80//00_LOGS

# Run workflow
anvi-run-workflow -w ecophylo -c ecophylo_config_GTDB_reps_80.json --skip-dry-run --additional-params --cluster \"clusterize -j={rule} -o={log} -n={threads} --exclude \'\' -x \" --jobs=120 --resources nodes=120 --latency-wait 100 --keep-going --rerun-incomplete

# Archaea
# cp config from bacteria and modify
cp ecophylo_config_GTDB_reps.json ecophylo_config_GTDB_reps_archaea_80.json

# Make log dir
mkdir -p ECOPHYLO_WORKFLOW_GTDB_REPS_ARCHAEA_80/00_LOGS

# Run workflow
anvi-run-workflow -w ecophylo -c ecophylo_config_GTDB_reps_archaea.json --skip-dry-run --additional-params --cluster \"clusterize -j={rule} -o={log} -n={threads} --exclude \'\' -x \" --jobs=120 --resources nodes=120 --latency-wait 100 --keep-going --rerun-incomplete
```

**Step 2.** Extract matrix of `HMM-hits` across genome datasets

Here we will extract a matrix of `HMM-hits` for the HMM collections
(`Bacteria_71` and `Archaea_76`) over the GTDB representative genomes.

``` bash
cd /project2/meren/PROJECTS/ECOPHYLO/Surface_Ocean

sed 's|\.|_|' GTDB_representatives_refseq_external_genomes_archaea.txt  > GTDB_representatives_refseq_external_genomes_archaea_ed.txt

sed 's|\.|_|' GTDB_representatives_refseq_external_genomes.txt > GTDB_representatives_refseq_external_genomes_ed.txt

# archaea
anvi-script-gen-hmm-hits-matrix-across-genomes -e GTDB_representatives_refseq_external_genomes_archaea_ed.txt \
                                                    --hmm-source Archaea_76 \
                                                  -o GTDB_representatives_refseq_external_genomes_archaea_GENOME_MATRIX.txt

# bacteria
anvi-script-gen-hmm-hits-matrix-across-genomes -e GTDB_representatives_refseq_external_genomes_ed.txt \
                                                    --hmm-source Bacteria_71 \
                                                  -o GTDB_representatives_refseq_external_genomes_bacteria_GENOME_MATRIX.txt
```

**Step 3.** Count copy number of SCGs across genomes

``` r
# load data
SCG_MATRIX_ARCHAEA <- read_tsv("DATA/GTDB_representatives_refseq_external_genomes_archaea_GENOME_MATRIX.txt")
SCG_MATRIX_BACTERIA <- read_tsv("DATA/GTDB_representatives_refseq_external_genomes_bacteria_GENOME_MATRIX.txt")

# make ggplot theme
plot_theme <- theme_light() +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, family = "Helvetica Neue Condensed Bold"),
                      axis.ticks.x = element_blank(),
                      legend.title = element_blank(),
                      axis.text = element_text(size = 12, family = "Helvetica Neue Condensed Bold"),
                      axis.title = element_text(size = 18, family = "Helvetica Neue Condensed Bold"),
                      text = element_text(size = 8,  family = "Helvetica Neue Condensed Bold"),
                      plot.title = element_text(hjust = 0.5, size = 18, family = "Helvetica Neue Condensed Bold"),
                      legend.position = "bottom")

# What percent of genomes does each SCG cover?
plot_num_scg <- function(X, TITLE, XAXIS, YAXIS) {
  
  ###
  # X <- SCG_MATRIX_ARCHAEA
  ###
  
  SCG_MATRIX_long <- X %>% pivot_longer(-genome_or_bin) %>% rename(SCG = "name")
  
  # Get order of bar charts
  percent_order <- X %>% 
    pivot_longer(-genome_or_bin) %>% 
    rename(SCG = "name") %>%
    mutate(presence = case_when(value >= 1 ~ 1,
                              value < 1  ~ 0)) %>%
    group_by(SCG) %>%
    summarize(total = sum(presence)) %>%
    mutate(percent = (total/SCG_MATRIX_long$genome_or_bin %>% unique() %>% length())*100)
  
   
  num_scgs <- SCG_MATRIX_long %>% .$SCG %>% unique() %>% length()
  num_genomes <- SCG_MATRIX_long %>% .$genome_or_bin %>% unique() %>% length()
   
  SCG_MATRIX_long %>%
    group_by(SCG) %>%
    count(value) %>%
    mutate(value = as.character(value),
           once = n[value == "1"]) %>%
    left_join(percent_order) %>%
    ggplot(aes(x =  fct_reorder(SCG, percent), y = n/num_genomes, fill = value)) +
    geom_bar(stat="identity") +
    plot_theme + 
    ylab(XAXIS) +
    xlab(YAXIS) + 
    ggtitle(glue(TITLE, " (n = {num_genomes} genomes)"))
}

SCG_count_archaea <- plot_num_scg(SCG_MATRIX_ARCHAEA, 
             XAXIS = "Percent detection GTDB Archaea genomes",
             YAXIS = "Archaea_76 SCG Pfam models",
             TITLE = "Distribution of SCG hits per GTDB Archaea representative genome")

SCG_count_bacteria <- plot_num_scg(SCG_MATRIX_BACTERIA, 
             XAXIS = "Percent detection GTDB Bacteria genomes",
             YAXIS = "Bacteria_71 Pfam models",
             TITLE = "Distribution of SCG hits per GTDB Bacteria representative genome")

plot_final <- ggarrange(SCG_count_bacteria, 
                        SCG_count_archaea, 
                        nrow = 2,
                        labels = c("A", "B"), font.label = (family = "Helvetica Neue Condensed Bold")) %>%
  annotate_figure(top = text_grob("SCG HMM detection across GTDB r95 representative genomes", size = 24, family = "Helvetica Neue Condensed Bold"))

plot_final
```

![](REPRODUCIBLE_WORKFLOW_files/figure-commonmark/unnamed-chunk-9-1.png)

## Benchmarking EcoPhylo workflow with Ribosomal proteins using CAMI synthetic metagenomes

Before analyzing the oral, ocean, and human gut metagenomic dataset, we
selected ribosomal proteins that contextualized the majority of the
associated genome collections. Additionaly, we examined the assembly
rate of ribosomal proteins in metagenomic assemblies to make sure any
candidate proteins were not over assembled. This was show in figures:
`Supplemental_8.png`, `Supplemental_8.png`, FIXME. Here is a code
outline to reproduce these analyses.

We explored different ribosomal protein clustering thresholds and their
impact on (non-specific read
recruitment)\[https://anvio.org/vocabulary/#non-specific-read-recruitment\]
in the ecophylo workflow. To do this, we used
[CAMI](https://doi.org/10.1038/s41592-022-01431-4) and performed a grid
search across ribosomal protein clustering thresholds from 95%-100%.
Here is how we did it with the
[CAMI](https://doi.org/10.1038/s41592-022-01431-4) marine dataset.

**Step 1.** Identify ribosomal proteins that contextualize the majority
of the genome collection by running the ecophylo workflow over the
genome collection and calculating what percentage of genomes have which
SCG.

**Step 2.** Make config files for ribosomal protein clustering grid
search from 95%-100%

``` bash
# Get default ecophylo config file
anvi-run-workflow -w ecophylo --get-default-config ecophylo.json

# Edit config file to perform grid search
DATASET="MARINE"
for DIR_SCG in RP_L22 RP_L2 RP_L5 RP_15 RP_S8
do
  for i in 95 96 97 98 99 1
  do
    percent_id=$(echo "${i}" | sed 's/0\.//')
    python /SCRIPTS/paramaterize_clustering_threshold.py --input-config-json ecophylo.json \
                                                         --output-config-json ecophylo_"${DIR_SCG}"_"${percent_id}"_"${DATASET}".json \
                                                         --min-seq-id $i \
                                                         --cov-mode 1 \
                                                         --hmm-list hmm_list_"${DIR_SCG}".txt \
                                                         --SCG "${DIR_SCG}" \
                                                         --cami-dataset "${DATASET}";
  done
done
```

**Step 3.** Run ecophylo workflow over grid search

``` bash
DATASET="MARINE"
for DIR_SCG in RP_L22 RP_L2 RP_L5 RP_15 RP_S8
do
  mkdir -p ECOPHYLO_WORKFLOW_95_PERCENT_SIM_"${DIR_SCG}"_"${DATASET}"/00_LOGS
  anvi-run-workflow -w ecophylo -c ecophylo_"${DIR_SCG}"_95_"${DATASET}".json --skip-dry-run --additional-params --cluster \"clusterize -j={rule} -o={log} -n={threads} --exclude \'\' -x \" --jobs=120 --resources nodes=120 --latency-wait 100 --keep-going --rerun-incomplete
done
```

**Step 4.** Collect non-specific read recruitment data and plot results

``` bash
# NOTE: to speed this up, we ran the first workflow at 95% clustering then cp over the directory structure so that the following workflows only have to repeat the clustering step
DATASET="Marine"
for DIR_SCG in RP_L22
do
  mkdir -p ECOPHYLO_WORKFLOW_95_PERCENT_SIM_"${DIR_SCG}"_"${DATASET}"/00_LOGS
  anvi-run-workflow -w ecophylo -c ecophylo_"${DIR_SCG}"_95_"${DATASET}".json --skip-dry-run --additional-params --cluster \"clusterize -j={rule} -o={log} -n={threads} --exclude \'\' -x \" --jobs=120 --resources nodes=120 --latency-wait 100 --keep-going --rerun-incomplete

# cp first directory over to new clustering threshold
DATASET="MARINE"
for DIR_SCG in RP_L22 RP_L2 RP_L5 RP_15 RP_S8
do
  for i in 96 97 98 99 1
  do
      mkdir -p ECOPHYLO_WORKFLOW_"${i}"_PERCENT_SIM_"${DIR_SCG}"_"${DATASET}"/00_LOGS
      cp -r ECOPHYLO_WORKFLOW_95_PERCENT_SIM_"${DIR_SCG}"_"${DATASET}"/01_REFERENCE_PROTEIN_DATA ECOPHYLO_WORKFLOW_"${i}"_PERCENT_SIM_"${DIR_SCG}"_"${DATASET}"
  done
done

# Run workflow for the rest of the clustering thresholds
DATASET="Marine"
for DIR_SCG in RP_L22 RP_L2 RP_L5 RP_15 RP_S8
do
  for i in 96 97 98 99 1
  do
      echo -e "ECOPHYLO_WORKFLOW_"${i}"_PERCENT_SIM_"${DIR_SCG}"_"${DATASET}""
      mkdir -p ECOPHYLO_WORKFLOW_"${i}"_PERCENT_SIM_"${DIR_SCG}"_"${DATASET}"/00_LOGS
      anvi-run-workflow -w ecophylo -c ecophylo_"${DIR_SCG}"_"${i}"_"${DATASET}".json --additional-params --cluster \"clusterize -j={rule} -o={log} -n={threads} --exclude \'\' -x \" --jobs=120 --resources nodes=120 --latency-wait 100 --keep-going --rerun-incomplete
  done
done


# Collect mismapping data

MISMAPPING_FILE_NAME="mismapping_MARINE.tsv"
DATASET="MARINE"
echo -e "SCG\tclustering_threshold\tsample\tTotal_reads_mapped\tMAQ_lessthan_2\tnum_mismapping" > "$MISMAPPING_FILE_NAME"

# Loop through the desired directories and collect missmapping datat
for DIR_SCG in RP_L22 RP_L2 RP_L5 RP_15 RP_S8
do
  for i in 95 96 97 98 99 1
  do
    RP=$(echo "$DIR_SCG" | sed 's|RP_||')
    MAPPING_PATH="/project2/meren/PEOPLE/mschechter/CAMI/${DATASET}/ECOPHYLO_WORKFLOW_${i}_PERCENT_SIM_${DIR_SCG}_${DATASET}/METAGENOMICS_WORKFLOW/04_MAPPING/Ribosomal_${RP}"

    # Loop through the BAM files in the directory
    for sample in "$MAPPING_PATH"/*bam
    do
      SAMPLE_NAME=$(basename "$sample")
      NUM_MISMAPPING=$(samtools view "$sample" | grep XS:i | cut -f12-13 | sed 's/..:i://g' | awk '$1==$2' | wc -l) # magic command to collect missmapping data
      TOTAL_NUM_READS_MAPPED=$(samtools view "$sample" | wc -l)
      MAPQ_lessthan_2=$(samtools view "$sample" | awk '$5<2' | wc -l)
      echo -e "$DIR_SCG\t$i\t$SAMPLE_NAME\t$TOTAL_NUM_READS_MAPPED\t$MAPQ_lessthan_2\t$NUM_MISMAPPING" >> "$MISMAPPING_FILE_NAME"
    done
  done
done
```

Percent of missmapping reads

``` r
plot_theme_2 <- theme_light() +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, family = "Helvetica Neue Condensed Bold"),
                      axis.ticks.x = element_blank(),
                      legend.title = element_blank(),
                      axis.text = element_text(size = 8, family = "Helvetica Neue Condensed Bold"),
                      axis.title = element_text(size = 8, family = "Helvetica Neue Condensed Bold"),
                      text = element_text(size = 8,  family = "Helvetica Neue Condensed Bold"),
                      plot.title = element_text(hjust = 0.5, size = 8, family = "Helvetica Neue Condensed Bold"),
                      legend.position = "bottom")

plot_missmapping <- function(PATH, TITLE){
  
  ####
  # PATH = "/DATA/mismapping_MARINE.tsv"
  # TITLE = "Marine"
  ####
  mismapping_MARINE <- read_tsv(PATH)

  mismapping_MARINE_ed <- mismapping_MARINE %>%
    mutate(percent_mismapping = (num_mismapping/Total_reads_mapped),
           clustering_threshold = as_factor(clustering_threshold))
  
  mismapping_MARINE_ed$clustering_threshold <- factor(mismapping_MARINE_ed$clustering_threshold, levels = c('95','96', '97', '98', '99', '1'), ordered = TRUE)
  
  mismapping_MARINE_ed %>%
    ggplot(aes(x = clustering_threshold, y = percent_mismapping, fill = SCG)) +
      geom_boxplot() +
      scale_y_continuous(labels = scales::percent) +
      # coord_cartesian(ylim = c(0, 100)) +
      plot_theme_2 +
      # guides(shape = guide_legend(override.aes = list(size = 0.5))) +
      ggtitle(TITLE)
}

# Plot missmapping
MARINE_missmapping <- plot_missmapping("DATA/mismapping_MARINE.tsv", TITLE = "Marine")
STRAIN_MADNESS_missmapping <- plot_missmapping("DATA/mismapping_STRAIN_MADNESS.tsv", TITLE = "STRAIN MADNESS")
PLANT_ASSOCIATED_missmapping <- plot_missmapping("DATA/mismapping_PLANT_ASSOCIATED.tsv", TITLE = "PLANT ASSOCIATED")

missmapping_p <- ggarrange(MARINE_missmapping, STRAIN_MADNESS_missmapping, PLANT_ASSOCIATED_missmapping,
          ncol = 3, nrow = 1, labels = c("A", "B", "C"))

missmapping_p
```

![](REPRODUCIBLE_WORKFLOW_files/figure-commonmark/unnamed-chunk-13-1.png)

## Human oral microbiome genome recovery with ribosomal protein rpL19

To explore genome recovery from the human oral cavity, we contextualized
a genomic collection with a set of metagenomes from tongue and plaque
samples. The first step of this exploration was to identify ribosomal
proteins that occur in the majority of the genome collection in
single-copy.

### Find top ribosomal proteins to genomic collection

Here we tabulate the copy number of ribosomal proteins across the
genomic dataset and identified that rpL19, rpS15, and rpS2 contextualize
the majority of the genomic collection.

1.  Run EcoPhylo workflow to annotate HOMD with Bacteria_71 SCG
    collection and extract HMM-hits

``` bash
# Get default config
anvi-run-workflow -w ecophylo --get-default-config ecophylo_config_survey_Bacteria_71.json

# Run workflow
mkdir -p ECOPHYLO_WORKFLOW/00_LOGS

anvi-run-workflow -w ecophylo -c ecophylo_config_survey_Bacteria_71.json --skip-dry-run --additional-params --cluster \"clusterize -j={rule} -o={log} -n={threads} -x \" --jobs=120 --resources nodes=120 --latency-wait 100 --keep-going --rerun-incomplete
```

Find SCGs that detect the majority of the genomic dataset

``` bash
# Across metagenome dataset
anvi-script-gen-hmm-hits-matrix-across-genomes -e metagenomes.txt --hmm-source Bacteria_71 -o Bacteria_71_Oral_METAGENOME_MATRIX.txt

# Across Genome dataset
anvi-script-gen-hmm-hits-matrix-across-genomes -e external_genomes_MAGs_HOMD.txt --hmm-source Bacteria_71 -o Bacteria_71_Oral_GENOME_MATRIX.txt

# JUST MAGs
anvi-script-gen-hmm-hits-matrix-across-genomes -e external-genomes.txt --hmm-source Bacteria_71 -o Bacteria_71_Oral_GENOME_MATRIX_just_MAGs.txt
```

2.  Plot SCG detection across metagenomic assemblies

Plot distribution of ribosomal proteins across metagenomic assemblies.

``` r
Bacteria_71_METAGENOME_MATRIX <- read_tsv("DATA/ORAL/Bacteria_71_Oral_METAGENOME_MATRIX.txt")

plot_metagenome_SCGs_frequencies <- plot_metagenome_SCGs(Bacteria_71_METAGENOME_MATRIX, 
                                                         YLAB = "Single-copy core genes (anvi'o Bacteria_71)",
                                                         XLAB = "Total SCG counts",
                                                         TITLE = "oral metagenomes")
```

3.  Plot SCG detection across genome dataset

``` r
Bacteria_71_GENOME_MATRIX <- read_tsv("DATA/ORAL/Bacteria_71_Oral_GENOME_MATRIX.txt")

# Make list of MAGs that are 50% complete, 10% redundant
MAG_metadata <- read_excel("DATA/ORAL/13059_2020_2195_MOESM2_ESM.xlsx", sheet = "(b) Non_redundant_DB")
list_50_compete <- MAG_metadata %>% 
  filter(percent_completion > 50, percent_redundancy < 10) %>%
  mutate(name = gsub("-", "_", MAGs_Id)) %>%
  .$name

MAGs_to_remove <- MAG_metadata %>% 
  filter(percent_completion <= 50) %>%
  mutate(name = gsub("-", "_", MAGs_Id)) %>%
  .$name

# Separate out MAGs and HOMD genomes
Bacteria_71_GENOME_MATRIX_MAGs <- Bacteria_71_GENOME_MATRIX %>% filter(grepl("ORAL_", genome_or_bin))
Bacteria_71_GENOME_MATRIX_MAGs_HQ <- Bacteria_71_GENOME_MATRIX %>% filter(grepl("ORAL_", genome_or_bin), grepl("MAG", genome_or_bin))
Bacteria_71_GENOME_MATRIX_HOMD <- Bacteria_71_GENOME_MATRIX %>% filter(!grepl("ORAL_", genome_or_bin))
Bacteria_71_GENOME_MATRIX_HOMD_MAGs_50 <- Bacteria_71_GENOME_MATRIX %>% filter(!genome_or_bin %in% MAGs_to_remove) # 50% complete MAGs
Bacteria_71_GENOME_MATRIX_MAGs_50 <- Bacteria_71_GENOME_MATRIX %>% filter(genome_or_bin %in% list_50_compete) # 50% complete MAGs

# How much of MAGs does the RP cover?
Bacteria_71_GENOME_MATRIX %>% 
  pivot_longer(-genome_or_bin) %>%
  rename(SCG = "name") %>%
  filter(SCG == "Ribosomal_L20") %>%
  filter(grepl("ORAL_", genome_or_bin)) %>%
    group_by(SCG) %>%
    count(value) %>% 
    mutate(percent_occ = n/sum(n)*100) %>%
    mutate(value = as.character(value),
           once = n[value == "1"])
```

    # A tibble: 3 × 5
    # Groups:   SCG [1]
      SCG           value     n percent_occ  once
      <chr>         <chr> <int>       <dbl> <int>
    1 Ribosomal_L20 0       470      59.5     319
    2 Ribosomal_L20 1       319      40.4     319
    3 Ribosomal_L20 2         1       0.127   319

``` r
# How much of Genomes does the RP cover?
Bacteria_71_GENOME_MATRIX %>% 
  pivot_longer(-genome_or_bin) %>%
  rename(SCG = "name") %>%
  filter(SCG == "Ribosomal_L20") %>%
  filter(!grepl("ORAL_", genome_or_bin)) %>%
    group_by(SCG) %>%
    count(value) %>% 
    mutate(percent_occ = n/sum(n)*100) %>%
    mutate(value = as.character(value),
           once = n[value == "1"])
```

    # A tibble: 3 × 5
    # Groups:   SCG [1]
      SCG           value     n percent_occ  once
      <chr>         <chr> <int>       <dbl> <int>
    1 Ribosomal_L20 0        67      0.778   8544
    2 Ribosomal_L20 1      8544     99.2     8544
    3 Ribosomal_L20 2         4      0.0464  8544

``` r
plot_bacteria_SCG_frequency_per_genomes <- plot_genome_SCGs(Bacteria_71_GENOME_MATRIX, 
                                                        YLAB = "Single-copy core genes (anvi'o Bacteria_71)",
                                                        XLAB = "Percent of genome collection detected",
                                                        TITLE = glue("Single-copy core gene percent frequencies across oral genome dataset (HOMD + MAGs) \n(n = {Bacteria_71_GENOME_MATRIX$genome_or_bin %>% length()})"))

plot_bacteria_SCG_frequency_per_MAGs <- plot_genome_SCGs(Bacteria_71_GENOME_MATRIX_MAGs, 
                                                        YLAB = "Single-copy core genes (anvi'o Bacteria_71)",
                                                        XLAB = "Percent of genome collection detected",
                                                        TITLE = glue("Single-copy core gene percent frequencies across Shaiber et al., 2020 NR MAGs \n(n = {Bacteria_71_GENOME_MATRIX_MAGs$genome_or_bin %>% length()})"))

plot_bacteria_SCG_frequency_per_MAGs_complete <- plot_genome_SCGs(Bacteria_71_GENOME_MATRIX_MAGs_HQ, 
                                                        YLAB = "Single-copy core genes (anvi'o Bacteria_71)",
                                                        XLAB = "Percent of genome collection detected",
                                                        TITLE = glue("Single-copy core gene percent frequencies across Shaiber et al., 2020 MAGs (completion > 70%) \n(n = {Bacteria_71_GENOME_MATRIX_MAGs_HQ$genome_or_bin %>% length()})"))

plot_bacteria_SCG_frequency_per_HOMD <- plot_genome_SCGs(Bacteria_71_GENOME_MATRIX_HOMD, 
                                                        YLAB = "Single-copy core genes (anvi'o Bacteria_71)",
                                                        XLAB = "Percent of genome collection detected",
                                                        TITLE = glue("Single-copy core gene percent frequencies across HOMD \n(n = {Bacteria_71_GENOME_MATRIX_HOMD$genome_or_bin %>% length()})"))

plot_bacteria_SCG_frequency_per_HOMD_MAGs_50 <- plot_genome_SCGs(Bacteria_71_GENOME_MATRIX_HOMD_MAGs_50, 
                                                        YLAB = "Single-copy core genes (anvi'o Bacteria_71)",
                                                        XLAB = "Percent of genome collection detected",
                                                        TITLE = glue("Single-copy core gene percent frequencies across oral genome dataset (HOMD + MAGs) \n(completion > 50%, redundancy < 10%; n = {Bacteria_71_GENOME_MATRIX_HOMD_MAGs_50$genome_or_bin %>% length()})"))

plot_bacteria_SCG_frequency_per__MAGs_50 <- plot_genome_SCGs(Bacteria_71_GENOME_MATRIX_MAGs_50, 
                                                        YLAB = "Single-copy core genes (anvi'o Bacteria_71)",
                                                        XLAB = "Percent of genome collection detected",
                                                        TITLE = glue("Single-copy core gene percent frequencies across Shaiber et al., 2020 MAGs \n(completion > 50%, redundancy < 10%; n = {Bacteria_71_GENOME_MATRIX_MAGs_50$genome_or_bin %>% length()})"))


supp_fig_8 <- ggarrange(plot_bacteria_SCG_frequency_per_HOMD_MAGs_50, plot_metagenome_SCGs_frequencies, ncol = 2)

supp_fig_8
```

![](REPRODUCIBLE_WORKFLOW_files/figure-commonmark/unnamed-chunk-17-1.png)

### Run the workflow with rpL19

``` bash
mkdir -p ECOPHYLO_WORKFLOW_RP_L19_FINAL/00_LOGS

anvi-run-workflow -w ecophylo -c ecophylo_config_SCGs_RP_L19.json --skip-dry-run --additional-params --cluster \"clusterize -j={rule} -o={log} -n={threads} -x \" --jobs=120 --resources nodes=120 --latency-wait 100 --keep-going --rerun-incomplete
```

#### Remove genomes that are not detected in any samples.

Some genomes in the genomic collections we analyzed were not detected in
the metagenomes. Ecophylo offers an efficient way to detect if a
populations is detected in a metagenome. We used this feature to quickly
filter out genomes that were not detected in the data. Subsquently, we
re-ran the workflow with only genomes that were detected.

``` bash
anvi-export-table PROFILE.db --table detection_splits -o detection_splits.txt
```

Here we will only work with splits that were detected with 90% coverage
in at least one metagenomes. We will then get the cluster members for
these splits and select the genomes that were detected. Finally we will
filter out the genomes that are not detected in any sample.

``` r
detect_and_write_genomes <- function(DIR_PATH, SCG, DETECTION_VALUE) {
  
  ###
  # DIR_PATH <- "DATA/ORAL/ECOPHYLO_WORKFLOW_RP_L19_FINAL"
  # SCG <- "Ribosomal_L19"
  # DETECTION_VALUE <- 0.9
  ###

  detection_splits <- read_tsv(path(DIR_PATH, "detection_splits.txt"))
  mmseqs_NR_cluster <- read_tsv(path(DIR_PATH, glue("{SCG}-mmseqs_NR_cluster.tsv")), col_names = c("representative", "cluster_member"))
  
  splits_0_detection <- filter_splits_by_detection(detection_splits, DETECTION_VALUE, SCG)

  genomes_not_detected <- mmseqs_NR_cluster %>%
    separate(representative, into = c("representative", "representative_split"), sep = glue("_{SCG}_")) %>%
    separate(cluster_member, into = c("cluster_member", "cluster_member_split"), sep = glue("_{SCG}_")) %>%
    select(representative, cluster_member) %>%
    filter(representative %in% splits_0_detection$representative_name) %>% # Filter for splits that were not detected in any sample
    select(cluster_member) %>%
    filter(grepl("^GCA_", cluster_member)) # Filter for only HOMD genomes
}

DIR_PATH <- "DATA/ORAL/ECOPHYLO_WORKFLOW_RP_L19_FINAL_cov_mode"
SCG <- "Ribosomal_L19"

genomes_NOT_detected_Ribosomal_L19 <- detect_and_write_genomes(DIR_PATH, SCG, DETECTION_VALUE = 0.9)

external_genomes_MAGs_HOMD <- read_tsv("DATA/ORAL/external_genomes_MAGs_HOMD.txt") 

external_genomes_MAGs_HOMD %>%
  filter(!name %in% genomes_NOT_detected_Ribosomal_L19$cluster_member) %>%
  write_tsv(path(glue("DATA/ORAL/external_genomes_MAGs_HOMD_{SCG}_detected.txt")), col_names = TRUE)
```

#### Re-run Ecophylo with only detected genomes

``` bash
# cp config and edit paths
cp ecophylo_config_SCGs_RP_L19_cov_mode.json ecophylo_config_SCGs_RP_L19_cov_mode_detected.json

mkdir -p ECOPHYLO_WORKFLOW_RP_L19_FINAL_cov_mode_detected/00_LOGS

cp -r ECOPHYLO_WORKFLOW_RP_L19_FINAL_cov_mode/01_REFERENCE_PROTEIN_DATA ECOPHYLO_WORKFLOW_RP_L19_FINAL_cov_mode_detected/01_REFERENCE_PROTEIN_DATA

anvi-run-workflow -w ecophylo -c ecophylo_config_SCGs_RP_L19_cov_mode_detected.json --skip-dry-run --additional-params --cluster \"clusterize -j={rule} -o={log} -n={threads} -x \" --jobs=40 --resources nodes=40 --latency-wait 100 --keep-going --rerun-incomplete
```

### Make miscellaneous data

Ecophylo provides some basic miscilaneous infoirmation about each of the
ribosomal proteins detected in your data. However, anvi’o make it easy
to enrich the interface with extra data to get more insights. Here is an
example of leveraging input and various output files from the workflow
to create more metadata:

#### Import genome-types

Get default collection to make metadata

``` bash
cd DATA/ORAL/ECOPHYLO_WORKFLOW_RP_L19_FINAL_cov_mode_detected

anvi-export-collection -p PROFILE.db -C DEFAULT
```

Add genome type metadata

``` r
SCG <- "Ribosomal_L19"
SCG_suffix <- "L19"
DIR_PATH <- "DATA/ORAL/ECOPHYLO_WORKFLOW_RP_L19_FINAL_cov_mode_detected"
DEFAULT_collection <- read_tsv(path(DIR_PATH, "collection-DEFAULT.txt"), col_names = c("representative", "collection")) %>%
  mutate(representative = str_replace(pattern = "_split_00001",replacement = "", representative))
detected_mmseqs_NR_cluster <- read_tsv(path(DIR_PATH, glue("{SCG}-mmseqs_NR_cluster.tsv")), col_names = c("representative", "cluster_member")) %>%
  filter(representative %in% DEFAULT_collection$representative)
genomes <- read_tsv(glue("DATA/ORAL/external_genomes_MAGs_HOMD_{SCG}_detected.txt"))

# Make list of MAGs that are 50% complete, 10% redundant
MAG_metadata <- read_excel("DATA/ORAL/13059_2020_2195_MOESM2_ESM.xlsx", sheet = "(b) Non_redundant_DB")
list_50_compete <- MAG_metadata %>% 
  filter(percent_completion > 50, percent_redundancy < 10) %>%
  mutate(name = gsub("-", "_", MAGs_Id)) %>%
  .$name

make_genome_type_file <- function(genomes, mmseqs_NR_cluster, SCG) {
  
  ####
  # genomes <- genomes
  # mmseqs_NR_cluster <- detected_mmseqs_NR_cluster
  # SCG = "Ribosomal_L19"
  ####
  
  REFG <- genomes %>% filter(grepl("GCA_", name)) 
  MAG <- genomes %>% filter(!grepl("GCA_", name))
  
  genome_type_metadata <- mmseqs_NR_cluster %>%
    separate(cluster_member, into = c("name", "number"), sep = glue("_{SCG}_"), remove = FALSE) %>%
    dplyr::select(-number) %>% 
    mutate(genome_type = case_when(name %in% MAG$name  ~ "MAG",
                                   name %in% REFG$name ~ "REFG",
                                   TRUE ~ "METAG")) %>%
    dplyr::select(representative, genome_type) %>% 
    distinct() %>%
    pivot_wider(id_cols = "representative", names_prefix = "genome_type.", names_from = "genome_type", values_from = "genome_type") %>%
    mutate_at(vars(starts_with("genome_type")), ~ ifelse(!is.na(.x), TRUE, FALSE)) %>%
    mutate(representative = str_c(representative, "_split_00001"))
  
  # Find MAGs with completion > 70% (they are labeled with "_MAG")
  genome_type_metadata_HQ_MAG <- mmseqs_NR_cluster %>%
    separate(cluster_member, into = c("name", "number"), sep = glue("_{SCG}_"), remove = FALSE) %>%
    dplyr::select(-number) %>% 
    # mutate(HQ_MAG = case_when(grepl("_MAG", name)           ~ TRUE,
    #                           TRUE                          ~ FALSE)) %>%
    mutate(HQ_MAG = case_when(name %in% list_50_compete           ~ TRUE,
                              TRUE                                ~ FALSE)) %>%
    dplyr::select(representative, HQ_MAG) %>% 
    distinct() %>%
    pivot_wider(id_cols = "representative", names_prefix = "HQ_MAG.", names_from = "HQ_MAG", values_from = "HQ_MAG") %>%
    filter(HQ_MAG.TRUE == TRUE) %>%
    select(-HQ_MAG.FALSE) %>%
    mutate(representative = str_c(representative, "_split_00001"))
  
  genome_type_metadata %>% 
    mutate("HQ_MAG.TRUE" = case_when(representative %in% genome_type_metadata_HQ_MAG$representative ~ TRUE,
                                     TRUE                                                           ~ FALSE),
           contains_genome = case_when(HQ_MAG.TRUE == TRUE | genome_type.REFG == TRUE ~ TRUE,
                                      TRUE                                            ~ FALSE),
           metagenome_ONLY = case_when(genome_type.METAG == TRUE & genome_type.REFG == FALSE & HQ_MAG.TRUE == FALSE ~ TRUE,
                                       TRUE                                                                         ~ FALSE))
}

genome_type_metadata <- make_genome_type_file(genomes, mmseqs_NR_cluster = detected_mmseqs_NR_cluster, SCG = "Ribosomal_L19")

# Save the misc data
genome_type_metadata %>%
  write_tsv(path(DIR_PATH,"genome_types.txt"))
```

Import genome types misc data into anvi’o

``` bash
cd DATA/ORAL/ECOPHYLO_WORKFLOW_RP_L19_FINAL_cov_mode_detected

anvi-import-misc-data genome_types.txt \
                              -p PROFILE.db \
                              --target-data-table items \
                              --just-do-it
```

#### Re-calculate phylogenetic tree with IQtree

``` bash
SCG="Ribosomal_L19"

iqtree2 -s "${SCG}"_aligned_trimmed_filtered.fa -m WAG -B 1000 -pre "${SCG}"_trimmed_filtered_IQTREE_ultrafast_bootstrap -T 20
```

Add “\_split_0001” string to each tree leaf so it can be imported back
into anvi’o

``` r
# Import tree
SCG <- "L19"
DIR_PATH <- glue("DATA/ORAL/ECOPHYLO_WORKFLOW_RP_{SCG}_FINAL_cov_mode_detected")

add_split_string_to_tree(IN_PATH = path(DIR_PATH, glue("Ribosomal_{SCG}_trimmed_filtered_IQTREE_ultrafast_bootstrap.contree")),
                         OUT_PATH = path(DIR_PATH, glue("Ribosomal_{SCG}_trimmed_filtered_IQTREE_ultrafast_bootstrap_ed.nwk")))
```

    DATA/ORAL/ECOPHYLO_WORKFLOW_RP_L19_FINAL_cov_mode_detected/Ribosomal_L19_trimmed_filtered_IQTREE_ultrafast_bootstrap_ed.nwk

Import reformatted tree back into anvi’o as a new items-order

``` bash
cd DATA/ORAL/ECOPHYLO_WORKFLOW_RP_L19_FINAL_cov_mode_detected

SCG="Ribosomal_L19"
anvi-import-items-order -p PROFILE.db \
                        -i "${SCG}"_trimmed_filtered_IQTREE_ultrafast_bootstrap_ed.nwk \
                        --name IQTree
```

Visualize final tree

``` bash
anvi-interactive -p DATA/ORAL/ECOPHYLO_WORKFLOW_RP_L19_FINAL_cov_mode_detected/PROFILE.db \
                 -c DATA/ORAL/ECOPHYLO_WORKFLOW_RP_L19_FINAL_cov_mode_detected/CONTIGS.db
```

### Calculating genome recovery rates

#### anvi-estimate-scg-taxonomy

One way we calculate genome recovery rates of different genome
aquisistion methods in this manuscript is by calculating the proportion
of detected taxa by a certain method. To do this, the ecophylo workflow
leverages
[anvi-estimate-scg-taxonomy](https://anvio.org/help/main/programs/anvi-estimate-scg-taxonomy/)
to annotate ribsomal proteins with taxonomy from GTDB.

``` bash
SCG="S15"
cd DATA/ORAL/ECOPHYLO_WORKFLOW_RP_"${SCG}"_FINAL_cov_mode_detected

anvi-export-misc-data -p PROFILE.db --target-data-table items -o items.tsv
```

``` r
SCG <- "L19"
DIR_PATH <- glue("DATA/ORAL/ECOPHYLO_WORKFLOW_RP_{SCG}_FINAL_cov_mode_detected")
item_additional_data_L19 <- read_tsv(path(DIR_PATH, "items.tsv"))

find_genome_recovery_rate <- function(taxonomic_rank, item_additional_data) {
  
  #####
  # WDIR="DATA/ORAL/"
  # DIR=path(WDIR, "ECOPHYLO_WORKFLOW_RP_L19_FINAL_detected")
  # SCG="Ribosomal_L19"
  # 
  # item_additional_data <- read_tsv(path(DIR, "items.tsv"))
  # taxonomic_rank = "t_class"
  ######

  # These taxa had x < 50% MAG recovery rate
  name_taxa_rank <- taxonomic_rank
  Genome_recovery <- item_additional_data %>%
   group_by(!!sym(taxonomic_rank)) %>%
    summarize(total_detected = n(),
              num_MAGs = sum(HQ_MAG.TRUE == TRUE), # we are only counting MAGs with completion > 70%
              MAG_recovery_rate = num_MAGs/total_detected,
              num_REFGs = sum(genome_type.REFG == TRUE),
              REF_recovery_rate = num_REFGs/total_detected,
              num_Genomes = sum(metagenome_ONLY == FALSE),
              total_genome_recovery_rate = num_Genomes/total_detected) %>%
    dplyr::select(!!sym(taxonomic_rank), num_MAGs, num_REFGs, num_Genomes, total_detected, MAG_recovery_rate, REF_recovery_rate, total_genome_recovery_rate) %>%
    mutate(rank = name_taxa_rank) %>%
    rename(taxa_name = !!sym(taxonomic_rank))
  
  return(Genome_recovery)
}

RP_L19_genome_recovery_rates_SCG_taxonomy_list <- purrr::map(taxonomic_ranks_GTDB, ~find_genome_recovery_rate(item_additional_data_L19, taxonomic_rank = .x))

RP_L19_genome_recovery_rates_SCG_taxonomy <- bind_rows(RP_L19_genome_recovery_rates_SCG_taxonomy_list)

RP_L19_genome_recovery_rates_SCG_taxonomy
```

    # A tibble: 343 × 9
       taxa_name       num_M…¹ num_R…² num_G…³ total…⁴ MAG_r…⁵ REF_r…⁶ total…⁷ rank 
       <chr>             <int>   <int>   <int>   <int>   <dbl>   <dbl>   <dbl> <chr>
     1 Bacteria            163      59     199     272   0.599   0.217   0.732 t_do…
     2 <NA>                  4       0       4       5   0.8     0       0.8   t_do…
     3 Actinomycetota       21      11      27      39   0.538   0.282   0.692 t_ph…
     4 Bacillota            31      18      42      59   0.525   0.305   0.712 t_ph…
     5 Bacteroidota         46      10      50      58   0.793   0.172   0.862 t_ph…
     6 Campylobactero…       7       4      10      12   0.583   0.333   0.833 t_ph…
     7 Fusobacteriota        5       4       8      12   0.417   0.333   0.667 t_ph…
     8 Patescibacteria      35       6      39      52   0.673   0.115   0.75  t_ph…
     9 Pseudomonadota       16       6      21      34   0.471   0.176   0.618 t_ph…
    10 Spirochaetota         2       0       2       6   0.333   0       0.333 t_ph…
    # … with 333 more rows, and abbreviated variable names ¹​num_MAGs, ²​num_REFGs,
    #   ³​num_Genomes, ⁴​total_detected, ⁵​MAG_recovery_rate, ⁶​REF_recovery_rate,
    #   ⁷​total_genome_recovery_rate

#### Taxonomic binning

Bin taxa of interest in the interactive interface and export the
collections

``` bash
cd DATA/ORAL/ECOPHYLO_WORKFLOW_RP_L19_FINAL_cov_mode_detected

anvi-export-collection -p PROFILE.db --collection-name taxa --output-file-prefix taxa
anvi-export-collection -p PROFILE.db --collection-name all_Bacteria --output-file-prefix Bacteria
```

``` r
# read in data
taxa_collection <- read_tsv("DATA/ORAL/ECOPHYLO_WORKFLOW_RP_L19_FINAL_cov_mode_detected/taxa.txt", col_names = c("items", "bin"))
bacteria_collection <- read_tsv("DATA/ORAL/ECOPHYLO_WORKFLOW_RP_L19_FINAL_cov_mode_detected/Bacteria.txt", col_names = c("items", "bin"))

# Calculate genome recovery rates with taxonomic binning
SCG <- "L19"
DIR_PATH <- glue("DATA/ORAL/ECOPHYLO_WORKFLOW_RP_{SCG}_FINAL_cov_mode_detected")
item_additional_data_L19 <- read_tsv(path(DIR_PATH, "items.tsv"))

# Get genome recovery for bins
RP_L19_taxonomic_binning <- item_additional_data_L19 %>%
  left_join(taxa_collection) %>%
  group_by(bin) %>%
  summarize(total_detected = n(),
            num_MAGs = sum(HQ_MAG.TRUE == TRUE),
            MAG_recovery_rate = num_MAGs/total_detected,
            num_IGs = sum(genome_type.REFG == TRUE),
            IG_recovery_rate = num_IGs/total_detected,
            num_Genomes = sum(metagenome_ONLY == FALSE),
            total_genome_recovery_rate = num_Genomes/total_detected) %>%
  drop_na()

# Get genome recovery for ALL bacteria
RP_L19_taxonomic_binning_bacteria <- item_additional_data_L19 %>%
  left_join(bacteria_collection) %>%
  group_by(bin) %>%
  summarize(total_detected = n(),
            num_MAGs = sum(HQ_MAG.TRUE == TRUE),
            MAG_recovery_rate = num_MAGs/total_detected,
            num_IGs = sum(genome_type.REFG == TRUE),
            IG_recovery_rate = num_IGs/total_detected,
            num_Genomes = sum(metagenome_ONLY == FALSE),
            total_genome_recovery_rate = num_Genomes/total_detected)

RP_L19_taxonomic_binning %>% 
  bind_rows(RP_L19_taxonomic_binning_bacteria) %>%
  write_tsv("TABLES/Supplemental_table_ii.tsv")

Supplemental_table_ii <- read_tsv("TABLES/Supplemental_table_ii.tsv")
```

## Global surface ocean with rpL14

rpl14 rpS8 rpS11

``` bash
# make new config
anvi-run-workflow -w ecophylo --get-default-config ecophylo_config_SCGs_surface_ocean_RP_L14_50.json

# Make log dir
mkdir -p ECOPHYLO_WORKFLOW_survey_SCGs_surface_ocean_L14_50/00_LOGS

# Run whole workflow
anvi-run-workflow -w ecophylo -c ecophylo_config_SCGs_surface_ocean_RP_L14_50.json --skip-dry-run --additional-params --cluster \"clusterize -j={rule} -o={log} -n={threads} -x \" --jobs=120 --resources nodes=120 --latency-wait 120 --keep-going --rerun-incomplete
```

### Find top ribosomal proteins to genomic collection

Same as above.

### Remove undetected genomes

Here we will export the anvi’o detection table to identify genome in the
collection that were not identified in the metagenomes.

``` bash
cd ~/Google_drive_uchicago/00_MANUSCRIPTS/2024_Schechter_EcoPhylo/DATA/OCEAN/ECOPHYLO_WORKFLOW_survey_SCGs_surface_ocean_L14_50_cov_mode

anvi-export-table PROFILE.db --table detection_splits -o detection_splits.txt
```

Rewrite `external-genomes.txt` with only detected genomes.

``` r
detect_and_write_genomes <- function(WORKING_DIR_PATH, COLLECTION, DIR_PATH, SCG, SCG_abbrev, DETECTION_VALUE){
  
  ###
  # WORKING_DIR_PATH = "DATA/OCEAN"
  # COLLECTION = 70
  # DIR_PATH = glue("DATA/OCEAN/ECOPHYLO_WORKFLOW_survey_SCGs_surface_ocean_RP_L14_{COLLECTION}")
  # SCG = "Ribosomal_L14"
  # SCG_abbrev = "RP_L14"
  # DETECTION_VALUE = 0.9
  ###
  
  # load input data
  DIR_PATH = glue("DATA/OCEAN/ECOPHYLO_WORKFLOW_survey_SCGs_surface_ocean_{SCG_abbrev}_{COLLECTION}")
  detection_splits <- read_tsv(path(DIR_PATH, "detection_splits.txt"))
  mmseqs_NR_cluster <- read_tsv(path(DIR_PATH, glue("{SCG}-mmseqs_NR_cluster.tsv")), col_names = c("representative", "cluster_member"))
  external_genomes <- read_tsv(path(WORKING_DIR_PATH, glue("PAOLI_ET_AL_MAGs_greaterthan_30m_samples_{COLLECTION}_external_genomes.txt")))
  misc <- read_tsv(path(DIR_PATH, glue("{SCG}_misc.tsv")))
  all_completeness <- read_tsv("DATA/OCEAN/all_completeness.txt")
  
  # Here we will filter for splits that contain genomes in their clusters AND are not detected in the metagenomic dataset
  splits_0_detection <- filter_splits_by_detection(detection_splits, DETECTION_VALUE, SCG)
  
  # only filter for genomes
  splits_0_detection_genomes_only <- splits_0_detection %>%
    left_join(misc %>% rename(split = split_name)) %>%
    filter(genomic_seq_in_cluster == "yes") # filter for only clusters that contain genomes
  
  # Here we will get all of the cluster members too
  # Get cluster members of splits with 0 detection across dataset
  mmseqs_NR_cluster_ed <- mmseqs_NR_cluster %>%
    separate(representative, into = c("representative_name", "representative_split"), sep = glue("_{SCG}_"), remove = FALSE) %>%
    separate(cluster_member, into = c("cluster_member_name", "cluster_member_split"), sep = glue("_{SCG}_"), remove = FALSE) %>%
    select(-representative_split, -cluster_member_split) 
  
  cluster_members_to_filter <- mmseqs_NR_cluster_ed %>%
    filter(representative_name %in% splits_0_detection_genomes_only$representative_name) %>%
    left_join(all_completeness %>% rename(cluster_member_name = Genome)) %>%
    filter(!is.na(Project)) %>% 
    select(representative, cluster_member_name, Project, genome_type, `Anvio Completion`, `Anvio Redundancy`)
  
  # Remove all undetected genomes from external-genomes.txt and export
  external_genomes_filtered <- external_genomes %>%
    filter(!name %in% cluster_members_to_filter$cluster_member_name)
  
  # export filtered external-genomes.txt
  print(path(WORKING_DIR_PATH, glue("PAOLI_ET_AL_MAGs_greaterthan_30m_samples_{COLLECTION}_external_genomes_{SCG}_detection_filtered.txt")))
  external_genomes_filtered %>% write_tsv(path(WORKING_DIR_PATH, glue("PAOLI_ET_AL_MAGs_greaterthan_30m_samples_{COLLECTION}_external_genomes_{SCG}_detection_filtered.txt")))
  
  return(cluster_members_to_filter)
}

RP_L14_50_filtered_genomes <- detect_and_write_genomes(WORKING_DIR_PATH = "DATA/OCEAN", COLLECTION = 50, DIR_PATH = glue("DATA/OCEAN/ECOPHYLO_WORKFLOW_survey_SCGs_surface_ocean_RP_L14_{COLLECTION}_cov_mode"), SCG = "Ribosomal_L14", SCG_abbrev = "RP_L14", DETECTION_VALUE = 0.9)
```

    DATA/OCEAN/PAOLI_ET_AL_MAGs_greaterthan_30m_samples_50_external_genomes_Ribosomal_L14_detection_filtered.txt

``` r
# How many genome types were removed?
RP_L14_50_filtered_genomes %>%
  count(Project, genome_type)
```

    # A tibble: 8 × 3
      Project                 genome_type     n
      <chr>                   <chr>       <int>
    1 BGEO                    METAG           1
    2 Delmont_prochlorococcus REFG            8
    3 Delmont_SAR11           REFG            3
    4 GORG                    SAGS           18
    5 MARD                    REFG         1347
    6 MARD                    SAGS           87
    7 Perez                   SAGS            6
    8 TARA                    METAG         449

### Re-run with filtered external-genomes.txt

``` bash
# make new config
cp ecophylo_config_SCGs_surface_ocean_RP_L14_50_cov_mode.json ecophylo_config_SCGs_surface_ocean_RP_L14_50_cov_mode_detected.json

# Make log dir
mkdir -p ECOPHYLO_WORKFLOW_survey_SCGs_surface_ocean_L14_50_cov_mode_detected/00_LOGS

anvi-run-workflow -w ecophylo -c ecophylo_config_SCGs_surface_ocean_RP_L14_50_cov_mode_detected.json --skip-dry-run --additional-params --cluster \"clusterize -j={rule} -o={log} -n={threads} -x \" --jobs=120 --resources nodes=120 --latency-wait 120 --keep-going --rerun-incomplete
```

### Make miscellaneous data

Get default collection to make metadata

``` bash
anvi-export-collection -p PROFILE.db -C DEFAULT
  
anvi-export-misc-data -p PROFILE.db --target-data-table items -o items_metadata.txt
```

Get genome-type metadata

``` r
DIR <- "DATA/OCEAN/ECOPHYLO_WORKFLOW_survey_SCGs_surface_ocean_L14_50_cov_mode_detected"
Ribosomal_L14_mmseqs_NR_cluster <- read_tsv(path(DIR, "Ribosomal_L14-mmseqs_NR_cluster.tsv"), col_names = c("representative", "cluster_member"))
DEFAULT_collection <- read_tsv(path(DIR, "collection-DEFAULT.txt"), col_names = c("representative", "collection")) %>%
 mutate(representative = str_replace(pattern = "_split_00001",replacement = "", representative))
items_metadata <- read_tsv(path(DIR, "items_metadata.txt"))

Ribosomal_L14_mmseqs_NR_cluster <- Ribosomal_L14_mmseqs_NR_cluster %>%
  filter(representative %in% DEFAULT_collection$representative)

SAR11_external_genomes <- read_tsv("DATA/OCEAN/SAR11_external_genomes.txt")
prochlorococcus_external_genomes <- read_tsv("DATA/OCEAN/prochlorococcus_external_genomes.txt")
MARTINEZ_PEREZ_SAGS_external_genomes <- read_tsv("DATA/OCEAN/MARTINEZ_PEREZ_SAGS_external_genomes.txt")
PAOLI_ET_AL_MAGs_greaterthan_30m_samples_50_external_genomes <- read_tsv("DATA/OCEAN/PAOLI_ET_AL_MAGs_greaterthan_30m_samples_50_external_genomes.txt")

paoli_MAGs_ONLY <- PAOLI_ET_AL_MAGs_greaterthan_30m_samples_50_external_genomes %>%
  separate(name, into = c("project", "biosample", "genome_type", "code"), remove = FALSE) %>%
  filter(genome_type == "METAG")

all_genomes <- c(PAOLI_ET_AL_MAGs_greaterthan_30m_samples_50_external_genomes$name, SAR11_external_genomes$name, prochlorococcus_external_genomes$name, MARTINEZ_PEREZ_SAGS_external_genomes$name)


metagenome_only <- Ribosomal_L14_mmseqs_NR_cluster %>% # add in layer to show clusters that only contain sequences from the assembly
  separate(cluster_member, into = c("name", "number"), sep = "_Ribosomal_L14_", remove = FALSE) %>%
  select(-number) %>% 
  mutate(genome_type = case_when(name %in% paoli_MAGs_ONLY$name                      ~ "MAG",
                                 grepl("_SAGS_", name)                               ~ "SAG",
                                 grepl("REFG", name)                                 ~ "REFG",
                                 name %in% prochlorococcus_external_genomes$name     ~ "REFG",
                                 name %in% MARTINEZ_PEREZ_SAGS_external_genomes$name ~ "SAG",
                                 name %in% SAR11_external_genomes$name               ~ "REFG",
                                 # name %in% Delmont_et_al_2021_external_genomes$name ~ "MAG",
                                 TRUE                                                ~ "METAG")) %>%
  select(representative, genome_type) %>% 
  distinct() %>%
  group_by(representative)  %>%
  filter(all(genome_type == "METAG"))

genome_types <- Ribosomal_L14_mmseqs_NR_cluster %>%
  separate(cluster_member, into = c("name", "number"), sep = "_Ribosomal_L14_", remove = FALSE) %>%
  select(-number) %>% 
  mutate(genome_type = case_when(name %in% paoli_MAGs_ONLY$name                      ~ "MAG",
                                 grepl("_SAGS_", name)                               ~ "SAG",
                                 grepl("REFG", name)                                 ~ "REFG",
                                 name %in% prochlorococcus_external_genomes$name     ~ "REFG",
                                 name %in% MARTINEZ_PEREZ_SAGS_external_genomes$name ~ "SAG",
                                 name %in% SAR11_external_genomes$name               ~ "REFG",
                                 # name %in% Delmont_et_al_2021_external_genomes$name ~ "MAG",
                                 TRUE                                                ~ "METAG")) %>%
  select(representative, genome_type) %>% 
  distinct() %>%
  pivot_wider(id_cols = "representative", names_prefix = "genome_type.", names_from = "genome_type", values_from = "genome_type") %>%
  mutate_at(vars(starts_with("genome_type")), ~ ifelse(!is.na(.x), TRUE, FALSE)) %>%
  mutate(items = str_c(representative, "_split_00001")) 

genome_types_1 <- genome_types %>%
  mutate(genome_type.METAG_only = case_when(representative %in% metagenome_only$representative ~ TRUE,
                                            TRUE                                               ~ FALSE))

items_metadata %>%
  left_join(genome_types_1) %>%
  select(items, genomic_seq_in_cluster, genome_type.METAG_only, genome_type.MAG, genome_type.SAG, genome_type.REFG) %>%
  write_tsv(path(DIR, "genome_types.txt"))
```

Import the metadata

``` bash
# import the genome_types
anvi-import-misc-data genome_types.txt -p PROFILE.db --target-data-table items

# color the ocean samples by project
anvi-export-state -p PROFILE.db -s default -o default_states.json

python SCRIPTS/color_samples.py ../sample_colors.tsv default_states.json

anvi-import-state -p PROFILE.db -s default_states_new.json -n default
```

### Manually curate EcoPhylo treetree

In our study, we calculated phylogenetic trees of large collections of
ribosomal proteins (in the surface ocean microbiome, thousands!) which
included sequences from all domains of life, including plastids and
mitochondria. Due to the high-throughput nature of the workflow and the
chance of recruiting assembly artifacts, we manually inspected and
removed unusually long branches and recalculated the tree. Additionally,
in some scenarios, we removed the entire mitochondria signal to improve
the topology of the tree. Here are the basic steps to do this:

Step 1. Make collection of bad branches

- FIXME: we can add a picture with examples here

Step 2. Export collection and remove those sequences from the ribosomal
protein fasta file

``` bash
# Export default collection and collection of bad branches
anvi-export-collection -p PROFILE.db -C DEFAULT
anvi-export-collection -C bad_branches -p PROFILE.db --output-file-prefix bad_branches

# Clean fasta headers
grep -v Eukaryotes bad_branches.txt | cut -f 1 | sed 's|_split_00001||' > bad_branches_headers.txt

SCG="Ribosomal_S11"

anvi-script-reformat-fasta $SCG-AA_subset.fa --exclude-ids bad_branches_headers.txt -o $SCG-AA_subset_remove_long_seqs.fa

ALIGNMENT_PREFIX=""${SCG}"-AA_subset_remove_long_seqs_aligned_maxiters_2"

# Align with two iterations
clusterize "muscle -in $SCG-AA_subset_remove_long_seqs.fa -out "${ALIGNMENT_PREFIX}".faa -maxiters 2" -n 15 -o 00_LOGS/align.log

# Trim
clusterize "trimal -in "${ALIGNMENT_PREFIX}".faa -out "${ALIGNMENT_PREFIX}"_trimmed.faa -gappyout" -n 5 -o 00_LOGS/trim.log

# Calculate tree
clusterize "FastTree "${ALIGNMENT_PREFIX}"_trimmed_filtered.faa > "${ALIGNMENT_PREFIX}"_trimmed_filtered_FastTree.nwk" -n 15 -o 00_LOGS/FastTree.log
```

Step 3. Use `anvi-split` to remove bad brances from the EcoPhylo
interface

``` bash
SCG="Ribosomal_S11"

grep -v -f bad_branches_headers.txt collection-DEFAULT.txt | sed 's|EVERYTHING|EVERYTHING_curated|' > my_bins.txt

anvi-import-collection my_bins.txt -C curated -p PROFILE.db -c "${SCG}"-contigs.db
                        
anvi-split -C curated --bin-id EVERYTHING_curated -p PROFILE.db -c "${SCG}"-contigs.db --output-dir "${SCG}"_curated
```

Step 4. Add the string “\_split_00001” to each tree leaf so we can
import it back into the interface

``` r
add_split_string_to_tree <- function(IN_PATH, OUT_PATH) {
  
  # Import tree
  tree <- read.tree(IN_PATH)
  
  # Create DF with tree tipe metadata
  tree_tip_metadata <- tree$tip.label %>% 
    as_tibble() %>%
    rename(tip_label = value) %>%
    mutate(tip_label = str_c(tip_label, "_split_00001"))
  
  
  tree$tip.label <- tree_tip_metadata$tip_label
  
  # Write DF 
  print(OUT_PATH)
  write.tree(tree, file = OUT_PATH)
}


add_split_string_to_tree(IN_PATH = glue("{SCG}-AA_subset_remove_long_seqs_aligned_maxiters_2_trimmed_filtered_FastTree.nwk"),
                         OUT_PATH = glue("{SCG}-AA_subset_remove_long_seqs_aligned_maxiters_2_trimmed_filtered_FastTree_ed.nwk"))
```

Step 5. Import new tree and visualize

``` bash
TREE_NAME="FastTree_curated"
anvi-import-items-order -p "${SCG}"_curated/EVERYTHING_curated/PROFILE.db \
                        -i "${SCG}"-AA_subset_remove_long_seqs_aligned_maxiters_2_trimmed_filtered_FastTree_ed.nwk \
                        --name $TREE_NAME

anvi-interactive -p "${SCG}"_curated/EVERYTHING_curated/PROFILE.db \
                 -c "${SCG}"_curated/EVERYTHING_curated/CONTIGS.db
```

### Find taxa with low genome recovery rates

Notes:

- anvi-scg-taxonomy metadata to calculate general genome recovery rates
  with the caveat that a few sequences are NA’s for taxonomy
- make taxonomy bins from FastTree taxa bin version - taxa of interest
  which will pass along taxonomic annotation to NA’s
- use bins to remove chloroplast signal
- calculate genome recovery for taxonomic bins
- export and look for other taxa with low genome recovery rate with
  anvi-scg-taxonomy metadata

``` r
find_genome_recovery_rate <- function(taxonomic_rank, item_additional_data) {
  
  #####
  # WDIR="DATA/OCEAN/"
  # DIR=path(WDIR, "ECOPHYLO_WORKFLOW_survey_SCGs_surface_ocean_L14_70_filtered")
  # SCG="Ribosomal_L14"
  # 
  # item_additional_data <- read_tsv(path(DIR, glue("{SCG}_curated/EVERYTHING_curated/items.txt")))
  # taxonomic_rank = "t_class"
  ######

  # These taxa had x < 50% MAG recovery rate
  name_taxa_rank <- taxonomic_rank
  Genome_recovery <- item_additional_data %>%
   group_by(!!sym(taxonomic_rank)) %>%
    summarize(total_detected = n(),
              num_genomes = sum(genomic_seq_in_cluster == "yes"),
              genome_recovery_rate = num_genomes/total_detected,
              num_MAGs = sum(genome_type.MAG == TRUE),
              MAG_recovery_rate = num_MAGs/total_detected,
              num_SAGs = sum(genome_type.SAG == TRUE),
              SAG_recovery_rate = num_SAGs/total_detected,
              num_REFGs = sum(genome_type.REFG == TRUE),
              REF_recovery_rate = num_REFGs/total_detected) %>%
    dplyr::select(!!sym(taxonomic_rank), num_MAGs, num_SAGs, num_REFGs, num_genomes, total_detected, genome_recovery_rate, MAG_recovery_rate, SAG_recovery_rate, REF_recovery_rate) %>%
    # arrange(genome_recovery_rate) %>%
    mutate(rank = name_taxa_rank) %>%
    rename(taxa_name = !!sym(taxonomic_rank))
  
  return(Genome_recovery)
}

WDIR="DATA/OCEAN/"
DIR=path(WDIR, "ECOPHYLO_WORKFLOW_survey_SCGs_surface_ocean_L14_50_cov_mode_detected")
SCG="Ribosomal_L14" 

items <- read_tsv(path(DIR, glue("{SCG}_curated/EVERYTHING_curated/items.txt")))
collection_taxa <- read_tsv(path(DIR, glue("{SCG}_curated/EVERYTHING_curated/collection-taxa-IQtree.txt")), col_names = c("items", "bin"))
Bacteria_Archaea <- read_tsv(path(DIR, glue("{SCG}_curated/EVERYTHING_curated/Bacteria_Archaea.txt")), col_names = c("items", "bin"))

chloroplasts <- Bacteria_Archaea %>% filter(bin == "Chloroplasts") %>% .$items
eukaryotes <- Bacteria_Archaea %>% filter(bin == "Eukaryotes") %>% .$items
items_no_chloroplasts <- items %>% filter(!items %in% chloroplasts, !items %in% eukaryotes) # remove chloroplasts and eukaryotes

# anvi-scg-taxonomy genome recovery rates
RP_L14_genome_recovery_rates_SCG_taxonomy_list <- purrr::map(taxonomic_ranks_GTDB, ~find_genome_recovery_rate(items_no_chloroplasts, taxonomic_rank = .x))

RP_L14_genome_recovery_rates_SCG_taxonomy <- bind_rows(RP_L14_genome_recovery_rates_SCG_taxonomy_list)

num_NA <- RP_L14_genome_recovery_rates_SCG_taxonomy %>% filter(rank == "t_domain" & is.na(taxa_name)) %>% .$total_detected
total_archaea <- RP_L14_genome_recovery_rates_SCG_taxonomy %>% filter(rank == "t_domain" & taxa_name == "Archaea") %>% .$total_detected
total_bacteria <- RP_L14_genome_recovery_rates_SCG_taxonomy %>% filter(rank == "t_domain" & taxa_name == "Bacteria") %>% .$total_detected

percent_NA_SCG_taxonomy <- num_NA / (total_archaea + total_bacteria)
```

``` r
# Genome recovery rate by bins
bin_genome_recovery_rate <- items_no_chloroplasts %>%
  left_join(collection_taxa) %>% 
  select(items, bin, genome_type.METAG_only, genome_type.MAG, genome_type.SAG, genome_type.REFG, genomic_seq_in_cluster) %>% 
  group_by(bin) %>%
  summarize(total_detected = n(),
            num_genomes = sum(genomic_seq_in_cluster == "yes"),
            genome_recovery_rate = num_genomes/total_detected,
            num_MAGs = sum(genome_type.MAG == TRUE),
            MAG_recovery_rate = num_MAGs/total_detected,
            num_SAGs = sum(genome_type.SAG == TRUE),
            SAG_recovery_rate = num_SAGs/total_detected,
            num_REFGs = sum(genome_type.REFG == TRUE),
            REF_recovery_rate = num_REFGs/total_detected) %>%
  select(bin, num_MAGs, num_SAGs, num_REFGs, num_genomes, total_detected, genome_recovery_rate, MAG_recovery_rate, SAG_recovery_rate, REF_recovery_rate) %>%
  drop_na()

# Genome recovery rates only bacteria and archaea
RP_L14_bin_archaea_bacteria <- items %>% 
  left_join(Bacteria_Archaea) %>% 
  group_by(bin) %>%
  summarize(total_detected = n(),
            num_genomes = sum(genomic_seq_in_cluster == "yes"),
            genome_recovery_rate = num_genomes/total_detected,
            num_MAGs = sum(genome_type.MAG == TRUE),
            MAG_recovery_rate = num_MAGs/total_detected,
            num_SAGs = sum(genome_type.SAG == TRUE),
            SAG_recovery_rate = num_SAGs/total_detected,
            num_REFGs = sum(genome_type.REFG == TRUE),
            REF_recovery_rate = num_REFGs/total_detected) %>%
  select(bin, num_MAGs, num_SAGs, num_REFGs, num_genomes, total_detected, genome_recovery_rate, MAG_recovery_rate, SAG_recovery_rate, REF_recovery_rate)

# Phylogenetic binning results
bin_genome_recovery_rate %>% 
  filter(bin != "Archaea",
         bin != "Chloroplasts") %>%
  bind_rows(RP_L14_bin_archaea_bacteria) %>%
  write_tsv("TABLES/Supplemental_table_5.tsv")
```

## Hadza human gut microbiome with rpS15

### Run EcoPhylo

``` bash
mkdir -p ECOPHYLO_WORKFLOW_RP_S15_internal_cov_mode/00_LOGS

anvi-run-workflow -w ecophylo -c ecophylo_config_RP_S15_internal_cov_mode.json --skip-dry-run --additional-params --cluster \"clusterize -j={rule} -o={log} -n={threads} -x --exclude \'\' \" --jobs=160 --resources nodes=160 --latency-wait 100 --keep-going --rerun-incomplete
```

### Re-calculate phylogenetic tree with IQtree

``` bash
iqtree2 -s Ribosomal_S15_renamed.faa -m WAG -B 1000 -pre Ribosomal_S15_renamed_IQTREE_ultrafast_bootstrap -T 20
```

``` r
# Import tree

DIR="DATA/GUT/ECOPHYLO_WORKFLOW_RP_S15_internal_cov_mode/"
SCG="Ribosomal_S15" 

add_split_string_to_tree(IN_PATH = path(DIR, glue("{SCG}_renamed_IQTREE_ultrafast_bootstrap.contree")),
                         OUT_PATH = path(DIR, glue("{SCG}_renamed_IQTREE_ultrafast_bootstrap_ed.contree")))
```

    DATA/GUT/ECOPHYLO_WORKFLOW_RP_S15_internal_cov_mode/Ribosomal_S15_renamed_IQTREE_ultrafast_bootstrap_ed.contree

``` bash
anvi-import-items-order -p PROFILE.db \
                        -i Ribosomal_S15_renamed_IQTREE_ultrafast_bootstrap.contree \
                        --name IQtree
```

### Import sample categorical data

``` bash
SCG="S15"
cd DATA/GUT/ECOPHYLO_WORKFLOW_RP_"${SCG}"_internal_cov_mode

anvi-import-misc-data ../layers_additional_data.txt \
                        -p PROFILE.db \
                        --target-data-table layers

# color the ocean samples by project
anvi-export-state -p PROFILE.db -s default -o default_states.json

python ../../../SCRIPTS/color_samples.py ../sample_colors.tsv default_states.json

anvi-import-state -p PROFILE.db -s default_states_new.json -n default

anvi-interactive
```

### Find taxa with low genome recovery rate

``` bash
SCG="S15"
cd DATA/GUT/ECOPHYLO_WORKFLOW_RP_"${SCG}"_internal_cov_mode
  
anvi-export-misc-data -p PROFILE.db --target-data-table items -o items.tsv

anvi-export-collection -p PROFILE.db --collection-name bacteria_archaea --output-file-prefix bacteria_archaea_collection

anvi-export-collection -p PROFILE.db --collection-name taxa_IQtree --output-file-prefix low_high_MAG_recovery_rates
```

#### anvi-scg-taxonomy

``` r
SCG <- "S15"
DIR_PATH <- glue("DATA/GUT/ECOPHYLO_WORKFLOW_RP_{SCG}_internal_cov_mode")
item_additional_data_S15 <- read_tsv(path(DIR_PATH, "items.tsv"))

find_MAG_recovery_rate <- function(item_additional_data, taxonomic_rank) {
  
  ###
  # item_additional_data <- item_additional_data_S16
  # taxonomic_rank = "t_class"
  ###

  # These taxa had x < 50% MAG recovery rate
  name_taxa_rank <- taxonomic_rank
  MAG_recovery <- item_additional_data %>%
   group_by(!!sym(taxonomic_rank)) %>%
    summarize(total_detected = n(),
              num_MAGs = sum(genomic_seq_in_cluster == "yes"),
              MAG_recovery_rate = num_MAGs/total_detected) %>%
    rename(taxa_name = !!sym(taxonomic_rank)) %>%
    select(taxa_name, num_MAGs, total_detected, MAG_recovery_rate) %>%
    arrange(desc(total_detected)) %>%
    mutate(rank = name_taxa_rank)
  
  return(MAG_recovery)
}

RP_S15_genome_recovery_rates_SCG_taxonomy_list <- purrr::map(taxonomic_ranks_GTDB, ~find_MAG_recovery_rate(item_additional_data_S15, taxonomic_rank = .x))

RP_S15_genome_recovery_rates_SCG_taxonomy <- bind_rows(RP_S15_genome_recovery_rates_SCG_taxonomy_list) %>% mutate(SCG = "Ribosomal_S15")

RP_S15_genome_recovery_rates_SCG_taxonomy
```

    # A tibble: 2,323 × 6
       taxa_name         num_MAGs total_detected MAG_recovery_rate rank     SCG     
       <chr>                <int>          <int>             <dbl> <chr>    <chr>   
     1 Bacteria              1569           2219            0.707  t_domain Ribosom…
     2 <NA>                    13            137            0.0949 t_domain Ribosom…
     3 Archaea                  6              7            0.857  t_domain Ribosom…
     4 Bacillota             1225           1721            0.712  t_phylum Ribosom…
     5 <NA>                    13            137            0.0949 t_phylum Ribosom…
     6 Bacteroidota            96            121            0.793  t_phylum Ribosom…
     7 Pseudomonadota          66            111            0.595  t_phylum Ribosom…
     8 Actinomycetota          46             95            0.484  t_phylum Ribosom…
     9 Cyanobacteriota         47             62            0.758  t_phylum Ribosom…
    10 Verrucomicrobiota       24             29            0.828  t_phylum Ribosom…
    # … with 2,313 more rows

#### Taxonomic binning

``` r
# bacteria and archaea only
bacteria_archaea_collection <- read_tsv("DATA/GUT/ECOPHYLO_WORKFLOW_RP_S15_internal_cov_mode/bacteria_archaea_collection.txt", col_names = c("items", "bin"))

RP_S15_Bacteria_Archaea_recovery <- item_additional_data_S15 %>% 
  left_join(bacteria_archaea_collection) %>%
  group_by(bin) %>%
  summarize(total_detected = n(),
            num_MAGs = sum(genomic_seq_in_cluster == "yes"),
            MAG_recovery_rate = num_MAGs/total_detected) 

# phylogenetic bins
SCG <- "S15"
DIR_PATH <- glue("DATA/GUT/ECOPHYLO_WORKFLOW_RP_{SCG}_internal_cov_mode")
taxa_S15 <- read_tsv(path(DIR_PATH, "low_high_MAG_recovery_rates.txt"), col_names = c("items", "bin"))
item_additional_data_S15 <- read_tsv(path(DIR_PATH, "items.tsv"))

taxonomic_binning_S15 <- item_additional_data_S15 %>% 
  left_join(taxa_S15) %>%
  group_by(bin) %>%
  summarize(total_detected = n(),
            num_MAGs = sum(genomic_seq_in_cluster == "yes"),
            MAG_recovery_rate = num_MAGs/total_detected) 

RP_S15_Bacteria_Archaea_recovery %>% 
  filter(bin == "Bacteria") %>% 
  bind_rows(taxonomic_binning_S15) %>% 
  write_tsv("TABLES/Supplemental_table_ee.tsv")

taxonomic_binning_S15 <- read_tsv("TABLES/Supplemental_table_ee.tsv")
```
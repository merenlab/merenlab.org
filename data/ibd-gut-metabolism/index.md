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

Here is a list of links for quick access to the data described in our manuscript and on this page:

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

## Downloading the datapack for this reproducible workflow

We've pre-packaged a lot of the data and scripts that you need for this workflow into a DATAPACK. To download it, run the following in your terminal:

```bash
FIXME download
cd FIXME
```

We suggest working from within this datapack. Most of the commands and paths referenced below will assume that your working directory is the uppermost level of the datapack directory structure.

We also wanted to provide access to the main set of metagenome assemblies that we analyzed in the paper. This is a rather large dataset, so we separated it from the rest of the datapack to avoid overburdening your storage system unnecessarily. If you want to access these assemblies, you can download them into your datapack by running the following commands:

{:.warning}
The size of this metagenome dataset is **96 GB** (the archive alone is ~33 GB). Please make sure you have enough space to download it!

```bash
FIXME download
```

Here is a quick overview of the datapack structure:

```
VESELI_2023_DATAPACK/
  |
  |- 00_FASTQ_FILES                       ## this is where you can download metagenome FASTQ files
  |- 01_ALL_METAGENOME_DBS                ## this is where you will generate contigs databases for the gut metagenome assemblies
  |- VESELI_ET_AL_METAGENOME_CONTIGS_DBS  ## this holds the contigs databases we generated for the subset of samples used in the primary analyses of this paper
  |- TABLES                               ## this holds the important data tables
  |- SCRIPTS                              ## this contains scripts that you can run to reproduce some of the work described below
  |- MISC                                 ## this holds miscellanous files
  |- OUTPUT                               ## you will be generating some output in this directory
```

## Computational environment details

The bulk of analyses in this study were done using anvi'o version 7.1-dev (that is, the development version of anvi'o following the stable release v7.1). You can use anvi'o version 8.0 (once it is released) to reproduce our results, as all of the relevant code has been included as part of that stable release.

The only relevant difference between v7.1-dev and v8.0 (with respect to reproducing our results) is the default KEGG snapshot, which is newer in v8.0 than the version we used for the analyses in this paper. The choice of KEGG version affects the results of `anvi-run-kegg-kofams` and `anvi-estimate-metabolism`. In order to use the same version we did, you should run the following code to download the appropriate snapshot onto your computer into the directory `./KEGG_2020-12-23/` (you can change that path if you want):

```bash
anvi-setup-kegg-kofams --kegg-snapshot v2020-12-23 \
                       --kegg-data-dir KEGG_2020-12-23 
```

Whenever KEGG-related programs are used, you can make them use the appropriate KEGG version by adding `--kegg-data-dir KEGG_2020-12-23` (replacing that path with wherever you decided to store the KEGG data on your computer).


## Obtaining our initial dataset of public fecal metagenomes

This section covers the steps for acquiring and processing our initial set of publicly-available gut metagenomes. We downloaded, assembled and annotated 2,893 samples from 13 different studies. We wanted a large number of samples from various sources in order to evaluate our metabolic competency hypothesis across a wide diversity of cohorts from different geographical locations, age groups, hospital systems, and degress of healthiness. Note that this extensive dataset was later filtered to remove samples with low-sequencing depth (as described in the next section), and as a result, not all of these samples were utilized for the main analyses in our study. However, you can access the full list of the 2,893 samples that we considered in sheet (c) of [Supplementary Table 1](https://doi.org/10.6084/m9.figshare.22679080). 

{:.warning}
This section is computationally intensive and requires a lot of storage resources. If you want to reproduce this section, you should make sure that your high-performance computing system is prepared to shoulder the burden. Note that the large dataset covered here is only relevant to a few of the analyses described later, so there may not even be a need for you to go through this section at all. If you are only interested in reproducing the main analyses of the paper, the datapack at [FIXME LINK]() provides the final contigs databases for the relevant subset of samples, so you can skip this part. :)

### Criteria for sample selection and sample groups
We sought to obtain a large number of fecal metagenomes from healthy individuals and from individuals with IBD. We used the following criteria to search for studies offering such samples:

1) The study provided publicly-available shotgun metagenomes of fecal matter. These were usually stool metagenomes, but we also accepted luminal aspirate samples from the ileal pouch of patients that have undergone a colectomy with an ileal pouch anal anastomosis (from [Vineis et al. 2016](https://doi.org/10.1128/mBio.01713-16)), since these samples also represent fecal matter.
2) The study sampled from people living in industrialized countries. These countries have a higher incidence of IBD and share dietary and lifestyle tendencies that result in a similar gut microbiome composition when compared to developing countries. We found one study that sampled from both industrialized and developing countries ([Rampelli et al. 2015](https://doi.org/10.1016/j.cub.2015.04.055)), and from this study we used only those samples from industrialized areas.
3) The study included samples from people with IBD and/or they included samples from people without gastrointestinal (GI) disease or inflammation. The latter were often healthy controls from studies of diseases besides IBD, or from studies of treatments such as dietary interventions and antibiotics, and in these cases we included only the control samples in our dataset.
4) The study provided clear metadata for differentiating between case and control samples, so that we could accurately assign samples to the appropriate group.

These criteria led us to 13 studies of the human gut microbiome, which are summarized in sheet (a) of [Supplementary Table 1](https://doi.org/10.6084/m9.figshare.22679080). We almost certainly did not find all possible studies that fit our requirements, but we were sufficiently satisfied with the large number of samples and the breadth of human diversity encompassed by these studies, so we stopped there.

Each of the samples had to be assigned to a group based upon the general health status of the sample donor. In order to do this, we had to make some decisions about what to include (or not) in our characterization of "healthy" people. Not being clinical GI experts, we did our best with the metadata and cohort descriptions that were provided by each study (though there is room for disagreement here). We decided upon three groups: a healthy group of samples from people without gastrointestinal disease or inflammation, which contained the control samples from most of the studies; an IBD group of samples from people with a confirmed diagnosis of ulcerative colitis (UC), Crohn's disease (CD), or unclassified IBD; and an intermediate 'non-IBD' group of samples from people without a definite IBD diagnosis but who nevertheless may be presenting symptoms of GI distress or inflammation. [Lloyd-Price et al. 2019](https://doi.org/10.1038/s41586-019-1237-9) describes the criteria and justification for this last group quite eloquently:

<blockquote>
Subjects not diagnosed with IBD based on endoscopic and histopathologic findings were classified as ‘non-IBD’ controls, including the aforementioned healthy individuals presenting for routine screening, and those with more benign or non-specific symptoms. This creates a control group that, while not completely ‘healthy’, differs from the IBD cohorts specifically by clinical IBD status.

<div class="blockquote-author">Lloyd-Price et al. 2019</div>
</blockquote>

For the studies characterizing their controls as 'non-IBD' samples (there were 3), we applied the same grouping of their samples within our dataset in order to be consistent. We also decided to include samples from people with colorectal cancer (CRC) (from [Feng et al. 2015](https://doi.org/10.1038/ncomms7528)) in the 'non-IBD' group, since inflmmation in the GI tract (including from IBD) can promote the development of CRC ([Kraus and Arber, 2009](https://doi.org/10.1016/j.coph.2009.06.006)), though no IBD diagnoses were described for the individuals in the CRC study.

We note that we did not exclude samples from individuals with high BMI from the healthy group, because the study that included such samples ([Le Chatelier et al. 2013](https://doi.org/10.1038/nature12506)) had already excluded individuals with GI disease, diabetes, and other such conditions.

### Downloading the public metagenomes used in this study

The SRA accession number of each sample is listed in [Supplementary Table 1c](https://doi.org/10.6084/m9.figshare.22679080). We downloaded the samples from each contributing study individually, over time, using the [NCBI SRA toolkit](https://github.com/ncbi/sra-tools) and particularly the [`fasterq-dump` program](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump) that will download the FASTQ files for each given SRA accession. We then gzipped each read file to save on space (and you will see us refer to these gzipped FASTQ files later in the workflow). 

If you want to download all of the samples we used in this work, we recommend a similar strategy (but please keep in mind that the storage requirements for almost 3,000 metagenomes will be huge). Feel free to reach out to us for help if you need it. For convenience, we've provided a plain-text version of Table 1c in our DATAPACK, which can be accessed at the path `TABLES/00_ALL_SAMPLES_INFO.txt`. The last column of that file provides the SRA accessions that can be used for downloading each sample. We recommend downloading each read file to the paths listed in the `r1` and `r2` columns for each sample (don't forget the gzip step), as this will keep the file organization consistent with what we expect later and will minimize the need to alter file paths. The `00_FASTQ_FILES/` directory referenced in those paths should already exist in the datapack, ready for you to add files to it.

There is one exception to this strategy, and that is the study by [Quince et al. 2015](https://doi.org/10.1038/ajg.2015.357). There are no deposited sequences under [the NCBI BioProject for this study](https://www.ncbi.nlm.nih.gov/bioproject/270985). The SRA accession column for these samples in `TABLES/00_ALL_SAMPLES_INFO.txt` contains NaN values, and these rows should be skipped when using `fasterq-dump` to download samples. We accessed these metagenomes directly from the study authors.

{:.notice}
If you are paying close attention, you might notice that not all of the samples from each contributing study are included in our dataset. A few samples here and there were dropped due to errors during processing (i.e., failed assembly) or missing metadata.

### Metagenome processing: single assemblies and annotations

We used the [anvi'o metagenomic workflow](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#metagenomics-workflow), which makes use of the workflow management tool [snakemake](https://snakemake.readthedocs.io/en/stable/), for high-throughput assembly and annotation of our large dataset. The samples from each contributing study were processed using individual workflow runs with similar configurations. Here are the most important steps in the workflow that directly impact the downstream analyses: 

* quality filtering of sequencing reads using the [Minoche et al. 2011](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-11-r112) guidelines via the [`illumina-utils` package](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0066643), specifically the program `iu-filter-quality-minoche`
* single assembly with [IDBA-UD](https://academic.oup.com/bioinformatics/article/28/11/1420/266973). We used all default parameters except that we set the minimum contig length (`--min_contig`) to be 1000
* generation of an anvi'o contigs database (and gene-calling) for each assembly with `anvi-gen-contigs-database`
* annotation of single-copy core genes with `anvi-run-hmms`
* annotation of KEGG KOfams with `anvi-run-kegg-kofams`

(the workflow has other steps, namely read recruitment of each sample against its assembly and the consolidation of the resulting read mapping data into anvi'o profile databases, but these are not crucial for our downstream analyses in this paper.)

We provide an example configuration file (`MISC/config.json`) in the DATAPACK that can be used for reproducing our assemblies. To run the workflow, you simply create a 3-column `samples.txt` file containing the sample name, path to the R1 file, and path to the R2 file for each sample that you downloaded. An example file is described in our [workflow tutorial](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#samplestxt). In fact, you can derive this file from `TABLES/00_ALL_SAMPLES_INFO.txt` (assuming you have either followed our file naming/organization recommendations or updated those paths to reflect the names/organization you decided upon):

```bash
grep -v Vineis_2016 TABLES/00_ALL_SAMPLES_INFO.txt | cut -f 1-3 > MISC/samples.txt
```

Then, you can start the workflow with the following command (hopefully adapted for use on a high-performance computing cluster):

```bash
anvi-run-workflow -w metagenomics -c MISC/config.json
```

A few notes:
* We renamed the samples from each study to incorporate information such as country of origin (for healthy samples) or host diagnosis (for IBD samples) for better readability and downstream sorting. To match our sample names, your `samples.txt` file should use the same sample names that are described in Supplementary Table 1c (or the first column of `TABLES/00_ALL_SAMPLES_INFO.txt` in the DATAPACK). If you generated the `samples.txt` from that file, you are good to go.
* We used the default snapshot of KEGG data associated with anvi'o v7.1-dev, which can be downloaded onto your computer by running `anvi-setup-kegg-kofams --kegg-snapshot  v2020-12-23`, as described earlier. To exactly replicate the results of this study, the metagenome samples need to be annotated with this KEGG version by changing the `--kegg-data-dir` parameter (in the `anvi_run_kegg_kofams` rule of the config file) to point to this snapshot wherever it located on your computer
* The number of threads used for each rule is set in the config file. We conservatively set this number in the example `MISC/config.json` to be 1 for all rules, but you will certainly want to adjust these to take advantage of the resources of your particular system.
* We set the output directory for contigs databases to be `01_ALL_METAGENOME_DBS` to be consistent with the file organization described in `TABLES/00_ALL_SAMPLES_INFO.txt`

The steps in the workflow described above apply to all of the metagenome samples except for those from [Vineis et al. 2016](https://doi.org/10.1128/mBio.01713-16), which had to be processed differently since the downloaded samples contain merged reads (rather than paired-end reads described in R1 and R2 files, as in the other samples). Since the metagenomics workflow currently works only on paired-end reads, we had to run the assemblies manually. Aside from the lack of workflow, there are only two major differences in the processing of the 96 Vineis et al. samples:

* No additional quality-filtering was run on the downloaded samples, because the merging of the sequencing reads as described in the paper's methods section already included a quality-filtering step
* Single assembly of the merged reads was done with [MEGAHIT](https://academic.oup.com/bioinformatics/article/31/10/1674/177884), using all default parameters except for a minimum contig length of 1000

We wrote a loop to run an individual assembly on each sample from [Vineis et al. 2016](https://doi.org/10.1128/mBio.01713-16). This loop makes use of the sample names and paths as established in the `TABLES/00_ALL_SAMPLES_INFO.txt` file:

```bash
while read name path; do \
  megahit -r $path  \
      --min-contig-len 1000 \
      -t 7 \
      -o ${name}_TMP; \
done < <(grep Vineis_2016 TABLES/00_ALL_SAMPLES_INFO.txt | cut -f 1-2 | tail -n+2)
```

Once the assemblies were done, we extracted the final assembly files from each output directory, renamed them with the sample name, and put them all in one folder:

```bash
mkdir -p VINEIS_ASSEMBLIES
while read name path; do \
  mv ${name}_TMP/final.contigs.fa VINEIS_ASSEMBLIES/${name}.fasta
done < <(grep Vineis_2016 TABLES/00_ALL_SAMPLES_INFO.txt | cut -f 1-2 | tail -n+2)
rm -r *TMP/
```

Then, we were able to leverage the anvi'o contigs workflow to generate the contigs databases for each assembly and run the annotation steps. We've provided the relevant configuration file for this workflow (`MISC/vineis_config.json`) as well as the input file that lists the path to each assembly (`MISC/vineis_fasta.txt`) in the DATAPACK, and this is how you could run it for yourself:

```bash
anvi-run-workflow -w contigs -c MISC/vineis_config.json
```

The same notes about setting the KEGG data version and the number of threads that we detailed for the metagenomics workflow apply to this workflow as well.

### A final note on data organization

If you've elected to reproduce the metagenome download and processing described in this section, you will have ended up with a lot of samples, assemblies, and contigs databases on your computer system. You may have noticed that so far, we've been relying a lot on the organization of files as described in the file `TABLES/00_ALL_SAMPLES_INFO.txt`, which includes the paths to each read file (`r1` and `r2`) and to each contigs database (`contigs_db_path`). We will continue to rely on that organization in subsequent sections, so now is a good time to double-check that all of the file paths are correct.

If you decided to follow a different file naming and organization strategy, that is fine. Where the files are on your computer does not matter for following the remainder of this workflow, as long as you prepare a file for yourself that describes the correct paths to each sample's read files and contigs database, and you use that file in place of the `TABLES/00_ALL_SAMPLES_INFO.txt` file as needed.


## Selecting our final dataset of metagenomes

Before we began analyzing metabolism within the large dataset we compiled, we discovered that we had to reduce our sample set to ensure accurate calculations. One of the critical steps in comparing community-level copy numbers of metabolic pathways between microbial communities of differing richness is **normalization of these data by the size of the community**. Pathway copy numbers will naturally tend to be higher in metagenome assemblies that describe larger communities, so comparing these values is not very meaningful when community sizes are vastly different. The gut microbiomes of people with IBD tend to harbor much less diversity than healthy gut microbiomes, so this is certainly a problem in our case. We therefore came up with a strategy of normalizing the pathway copy numbers calculated for a given sample with the number of microbial populations represented within that metagenome assembly.

However, this normalization only works if we can accurately estimate that number of populations in each sample. Yet this is not the case for samples of low sequencing depth, as you will see later in this section - to a certain extent, the estimated number of populations is correlated with sample sequencing depth. We interpreted this to mean that low-depth samples fail to capture enough sequences from low-abundance populations, thereby skewing the population estimates in these assemblies. Therefore, we tried to mitigate the issue by filtering out low-depth samples and keeping only samples with high enough sequencing depth for our downstream analyses.

This section will cover how we estimate the number of populations in each metagenome assembly using single-copy core genes, our analysis of its relationship with sequencing depth, and the removal of low-depth samples to establish our final set of samples for analysis.

### Estimating number of populations per sample

To estimate how many microbial populations are represented in a metagenome assembly, we can rely on the fact that all microbial genomes (with few exceptions) contain exactly one copy of each gene in a special set of essential genes called single-copy core genes (SCGs). These include ribosomal proteins and other housekeeping genes. Anvi'o ships with a few generic sets of SCGs that is each specific to a domain of microbial life (Bacteria, Archaea, and Protista) and these genes are annotated using the program `anvi-run-hmms` (which, you might recall, we ran earlier as part of our metagenome processing workflows).

Since we expect to find one copy of each SCG in each microbial population, we can count the total number of copies of an SCG in a metagenome assembly and use that as the number of populations. However, using just one SCG for estimation would be error-prone due to missing SCG annotations in incomplete data (or the occasional duplication within one genome). Instead, we can use all of the SCGs for a given domain to make the estimate more robust to noise. The **mode of the number of SCGs** in the assembly gives us our estimate of the number of populations, in this case. The sketch below illustrates this process - each SCG is annotated in the metagenome assembly as indicated by the black boxes on the highlighted sequences on the bottom (shown for the first 3 SCGs only), these annotations are tallied (top histogram), and then the mode of the counts is computed:

[![Estimating the number of microbial populations using the mode of single-copy core gene annotations](images/estimate_pops_with_SCGs.png)](images/estimate_pops_with_SCGs.png){:.center-img .width-50}

The anvi'o codebase includes a `NumGenomesEstimator` class that does exactly this: takes the mode of the number of copies of the SCGs for a particular domain, giving you an estimate of the number of bacteria, archaea, and protists in a given metagenome assembly. Those values can be added together to obtain the total number of microbial populations in the sample (for gut metagenomes, usually only bacterial species are found).

We wrote a script that runs this estimation on each of our 2,893 samples. You can run it using the following command, and it will produce a table in the output directory (`OUTPUT/NUM_POPULATIONS_ALL_SAMPLES.txt`).

```bash
python SCRIPTS/estimate_num_genomes.py
```

If you didn't download all the metagenome samples, or you don't feel like running this part, you can access the estimates in the `num_populations` column of the `TABLES/00_ALL_SAMPLES_INFO.txt` file.

### Determining sequencing depth

We used a BASH loop to count the number of reads in each metagenome sample. It counts the number of lines in each (gzipped) FASTQ file and divides by 4 (the number of lines per read) to get the number of reads, which is then stored in a tab-delimited file. R1 files and R2 are counted separately. For the Vineis et al. samples, the count of merged reads per sample is stored in the R1 column, and the R2 column is 0 (because there is only one FASTQ file for each of these samples).

```bash
echo -e "name\tnum_reads_r1\tnum_reads_r2" > OUTPUT/metagenome_num_reads.txt
while read name r1 r2 db diagnosis study remainder; do \
  numr1=$(echo $(zcat $r1 | wc -l) / 4 | bc); \
  if [ "$study" = "Vineis_2016" ]; then \
    numr2=0; \
  else numr2=$(echo $(zcat $r2 | wc -l) / 4 | bc); \
  fi;
  echo -e "$name\t$numr1\t$numr2" >> OUTPUT/metagenome_num_reads.txt; \
done < <(tail -n+2 TABLES/00_ALL_SAMPLES_INFO.txt)
```

Since the number of files to process is so large, this will take quite a while to run. We've already stored the resulting read counts in the `TABLES/00_ALL_SAMPLES_INFO.txt` file in case you don't want to wait.

### Generating Supplementary Figure 1
Supplementary Figure 1 is a scatterplot demonstrating the correlation betwen sequencing depth and the estimated number of populations in a metagenome assembly. The code to plot this figure can be found in the R script at `SCRIPTS/plot_figures.R` in the DATAPACK. Here is the relevant code taken from the script (note: this code snippet does not include some required setup, like loading packages and setting some global variables, and will not run on its own).

```r
#### SUPP FIG 1 - SEQUENCING DEPTH SCATTERPLOT ####
all_metagenomes = read.table(file=paste(data_dir, "00_ALL_SAMPLES_INFO.txt", sep=""), 
                             header = TRUE, sep = "\t")
gplt = all_metagenomes %>%
  ggplot(aes(r1_num_reads/1e6, num_populations, color=group)) +
  geom_point(alpha=0.5, size=2) +
  labs(y="Estimated Number of Genomes", x="Number of Sequencing Reads (R1) * 10^6", subtitle="ALL Samples: Sequencing Depth vs # Genomes") +
  geom_vline(xintercept = 25) +
  scale_color_manual(values = c(HEALTHY_color, IBD_color, NONIBD_color)) +
  scale_fill_manual(values = c(HEALTHY_color, IBD_color, NONIBD_color)) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) +
  stat_cor(data=subset(all_metagenomes, r1_num_reads < 25000000), method='spearman', label.x.npc = 'left', show.legend = FALSE) + 
  geom_smooth(data=subset(all_metagenomes, r1_num_reads < 25000000), method='lm', aes(fill=group), alpha=.2) + 
  stat_cor(data=subset(all_metagenomes, r1_num_reads >= 25000000), method='spearman', label.x.npc = 'center', show.legend = FALSE) +
  geom_smooth(data=subset(all_metagenomes, r1_num_reads >= 25000000), method='lm', aes(fill=group), alpha=.2)
gplt
```

This produces the following plot (which we cleaned up Inkscape to produce the final polished figure for the manuscript):

[![Supplementary Figure 1. Scatterplot of sequencing depth vs estimated number of populations.](images/supp_fig_1.png)](images/supp_fig_1.png){:.center-img .width-50}

As you can see, the correlation between a sample's sequencing depth and the estimated number of microbial populations it contains is fairly strong, particularly for lower depth samples. That correlation starts to weaken at higher sequencing depths. We selected our sequencing depth threshold to be 25 million reads as a compromise between the need for accurate estimates of microbiome richness and the need for enough samples for a robust and powerful analysis.

The script at `SCRIPTS/plot_figures.R` contains code for most of the other figures in the manuscript (those that were generated from data). Not all of the code and figures will be highlighted in this workflow, but you can always find it in that file.

### Removal of samples with low sequencing depth

Running the following script will subset the samples with >= 25 million sequencing reads. It will generate a new table at `TABLES/00_SUBSET_SAMPLES_INFO.txt` containing the subset of 408 samples. You'll see some information about the resulting sample groups (and which studies contributed to them) in the output of the program.

```bash
python SCRIPTS/subset_metagenome_samples.py
```

### Final set of samples and their contigs DBs

In the remainder of the analyses described in our manuscript, we utilized the subset of 408 samples with high sequencing depth described in `TABLES/00_SUBSET_SAMPLES_INFO.txt`. Since 408 is a much more reasonable number than 2,893 we have provided the contigs databases for our metagenome assemblies of these samples in [this datapack](FIXME LINK). If you elect to download this datapack via the instructions at the start of this workflow, you will find the assemblies in the directory called `VESELI_ET_AL_METAGENOME_CONTIGS_DBS/`, and the table already contains the path to these samples (relative to the uppermost directory of the datapack). So, even if you didn't download and process all of the metagenome samples as described in the first section of this workflow, you can still continue with the subsequent sections.

## Metabolism analyses (metagenomes)

In this section, we will cover the analyses used to determine the set of metabolic pathways that are enriched in the IBD sample group. First, we calculated the copy numbers of all pathways in the KEGG modules database, in each metagenome assembly. Those numbers are made from the combined genes from all microbial populations represented in a given metagenome and thus quantify the community-level metabolic potential. To make them comparable across different gut communities of varying diversity, we normalized each copy number by the estimated number of populations in the same sample to obtain a 'per-population copy number', or PPCN. We then ran a statistical test on each module to determine which pathways were **the most different between the sample groups**, specifically looking for the ones that had **higher PPCN in IBD**.

### Metabolism estimation

We used the program `anvi-estimate-metabolism` to compute copy numbers of KEGG modules in each sample. You can find details about that program and its calculation strategies on [this page](https://anvio.org/help/main/programs/anvi-estimate-metabolism/). To ensure that all annotated KOfams in the metagenome contributed to the calculations, we ran the program in 'genome mode' on each sample, and to do this in a high-throughput manner, we provided the program with an [external genomes file](https://anvio.org/help/main/artifacts/external-genomes/) containing the paths to all the samples' contigs databases at once. Here is the code to both generate that input file (from the table containing the sample information, including paths) and run the metabolism estimation code.

```bash
echo -e "name\tcontigs_db_path" > OUTPUT/METAGENOME_PATHS.txt    # header for the input file
tail -n+2 TABLES/00_SUBSET_SAMPLES_INFO.txt | cut -f 1,4 >> OUTPUT/METAGENOME_PATHS.txt    # copy the sample names and paths over

# run metabolism estimation
anvi-estimate-metabolism -e OUTPUT/METAGENOME_PATHS.txt \
                    -O OUTPUT/METAGENOME_METABOLISM \
                    --kegg-data-dir KEGG_2020-12-23 \
                    --add-copy-number \
                    --matrix-format
```

Note the use of `--kegg-data-dir KEGG_2020-12-23`, which points the program to the version of KEGG data that we used for our analysis (which you should have downloaded when going through the 'Computational environment details' section). If the path to this data on your computer is different, you will have to change that in the command before you run it.

Since we used the `--matrix-format` flag, the output of `anvi-estimate-metabolism` will be a set of matrices, each containing a different statistic summarized across all of the metagenomes and all of the pathways. The one that we want for downstream analysis is the one containing stepwise copy numbers, which can be found at the path `OUTPUT/METAGENOME_METABOLISM-module_stepwise_copy_number-MATRIX.txt`. 

We've included our version of this output file in the DATAPACK. You will find it at `TABLES/METAGENOME_METABOLISM-module_stepwise_copy_number-MATRIX.txt` (and this table is referenced in the scripts below that depend on this data, just in case you don't run this step).

### Normalization of pathway copy numbers to PPCN

We just calculated the pathway copy numbers, and we estimated the number of microbial populations in each metagenome assembly in the previous section. To normalize the former, we divide by the latter. It is that simple :) 

The R script at `SCRIPTS/module_stats_and_medians.R` contains the code for computing PPCN values. It gets the module copy numbers and per-sample population estimates from the tables in the `TABLES/` folder, and generates a long-format data table at `OUTPUT/ALL_MODULES_NORMALIZED_DATA.txt` that includes the PPCN values in the `PPCN` column.

You can see the full code in the script, but here are the most relevent lines for your quick perusal:

```r
## NORMALIZING FUNCTION
normalize_values = function(matrix, normalizing_matrix, normalize_by_col){
  # function that normalizes the 'value' column
  # matrix: long-form dataframe containing the 'value' column and 'sample' column
  # normalizing_matrix: dataframe containing the data with which to normalize. 'sample' column must match samples in matrix
  # normalize_by_col: which column to normalize 'value' with. Must be passed in form normalizing_matrix$colname
  matrix$normalize_with = normalize_by_col[match(matrix$sample, normalizing_matrix$sample)]
  matrix$normalized_value = matrix$value / matrix$normalize_with
  return(matrix)
}

#### COMPUTING PPCN VALUES (COPY NUMBER NORMALIZATION) ####
## LOAD DATA
metagenomes = read.table(file=paste(data_dir, "00_SUBSET_SAMPLES_INFO.txt", sep=""), 
                         header = TRUE, sep = "\t")
module_info = read.table(paste(data_dir, "ALL_MODULES_INFO.txt", sep=""), header = TRUE, sep="\t")

## LOAD COPY NUMBER MATRIX AND NORMALIZE AND GROUP SAMPLES
table_from_file = read.table(file=paste(data_dir, "METAGENOME_METABOLISM-module_stepwise_copy_number-MATRIX.txt", sep=""), 
                             header = TRUE, sep = "\t")
stepwise_matrix = melt(table_from_file)
colnames(stepwise_matrix) = c('module', 'sample', 'value')
normalized_df = normalize_values(stepwise_matrix, normalizing_matrix = metagenomes, normalize_by_col = metagenomes$num_populations)
normalized_df$sample_group = metagenomes$group[match(normalized_df$sample, metagenomes$sample)]
```

Later on, we will also make use of the median per-population copy number for each module in each sample group (for visualizing the data as boxplots). The script also generates these median values into the file `OUTPUT/ALL_MODULES_MEDIAN_PPCN.txt`. If you plan to reproduce our figures, you should run that section of code, too (it is clearly marked in the comments).

### Enrichment analysis for IBD-enriched pathways

For each module, we used a one-sided Wilcoxon Rank-Sum test on its PPCN values to determine whether the module had significantly higher copy number in the IBD sample group compared to the healthy group. Here is the code we used to do that, which is coming from the script at `SCRIPTS/module_stats_and_medians.R`:

```r
#### PER-MODULE ENRICHMENT TEST ####
## SELECT GROUPS OF INTEREST (HEALTHY, IBD)
ibd_and_healthy = median_copy_num %>% filter(sample_group %in% c('HEALTHY', 'IBD'))

## PER-MODULE STAT TEST (IBD vs HEALTHY)
median_per_module = spread(ibd_and_healthy, key = sample_group, value = normalized_value)
median_per_module[, "p_value"] = NA
median_per_module[, "W_stat"] = NA
for (mod in levels(normalized_df$module)){
  ibd_norm_values = normalized_df %>% filter(sample_group == "IBD" & module == mod)
  healthy_norm_values = normalized_df %>% filter(sample_group == "HEALTHY" & module == mod)
  w2 = wilcox.test(ibd_norm_values$normalized_value, healthy_norm_values$normalized_value, alternative = 'greater')
  median_per_module$p_value[median_per_module$module == mod] = w2$p.value
  median_per_module$W_stat[median_per_module$module == mod] = w2$statistic
}

## FDR adjustment of p-value
median_per_module$fdr_adjusted_p_value = p.adjust(median_per_module$p_value, method = "fdr")
```

The resulting p-values were used to filter the pathways for those that were most enriched in the IBD gut microbiome. To be even more conservative in what we considered to be 'enriched' in IBD, we also required the module to have a minimum difference in PPCN values between the two sample groups ('effect size'). We calculated this difference by taking the median of the pathway's PPCN values in each sample group, and subtracting the healthy median from the IBD median.

```r
## DETERMINE SET OF IBD-ENRICHED MODULES
median_per_module$diff = median_per_module$IBD  - median_per_module$HEALTHY
adj_p_value_threshold = 2e-10

# check distribution of effect sizes after p-value filter
over_adj_p = median_per_module %>% filter(fdr_adjusted_p_value <= adj_p_value_threshold)
stripchart(over_adj_p$diff)
mean(over_adj_p$diff)

# set effect size threshold as mean of IBD-HEALTHY medians for modules over the p-value threshold
diff_threshold = mean(over_adj_p$diff)
over_adj_threshold = median_per_module %>% filter(fdr_adjusted_p_value <= adj_p_value_threshold & diff >= diff_threshold)
```

If you keep following the code in the script, you will see that it will generate a file with the statistical test results and a variety of other data for each module, including the enrichment status, at `OUTPUT/ALL_MODULES_MEDIANS_AND_STATS.txt`. It also prints a subset of this information for the IBD-enriched modules to the file at `OUTPUT/IBD_ENRICHED_MODULES.txt`.

One more note - you may notice that we exclude one module from the automatically-generated IBD-enriched set:

```r
## FINAL LIST OF IBD-ENRICHED MODULES
# exclude M00006 because it is the first part of M00004
ibd_enriched = ordered_modules %>% filter(module != "M00006")
```

This is because the module in question, [M00006](https://www.genome.jp/entry/M00006), is the oxidative phase of the pentose phosphate cycle and overlaps completely with the module describing the entirety of the pentose phosphate cycle, [M00004](https://www.genome.jp/entry/M00004). We also considered removing module [M00007](https://www.genome.jp/entry/M00007) because it is the non-oxidative phase of this cycle; however, it has some enzymatic differences with module M00004, so we left it in.

### Generating Figure 2

Figure 2 includes several plots of both unnormalized copy numbers and PPCN values. We won't copy the code to generate these plots here, but just wanted to remind you that you can find it in the script `SCRIPTS/plot_figures.R`. The section of code for each plot was designed to be an independent as possible from the other sections of code, so that you need only run that section to get what you want. Most of them do require the packages, paths, and color variables set up at the top of the script, however.

We used Inkscape to polish up the resulting plots into a publication-ready figure. In particular, we adjusted the background color opacity for the heatmaps in panels B and E to make the signal more visible, which is why you may notice a difference between the generated plots and the plots in the manuscript.

### Summary of metagenome results

So far, we've analyzed gut metagenomes, using per-population copy number as our metric for the typical metabolic capacity of a microbial population in the community described by each metagenome. We obtained a list of metabolic pathways that are enriched in the IBD sample group, that largely overlap with the modules associated with metabolic independence from [our previous study](https://doi.org/10.1186/s13059-023-02924-x), and that are mostly biosynthesis pathways for important cellular metabolites. These 33 IBD-enriched pathways could be important to microbial survival in the depleted communities of the IBD gut environment, and represent a refined list of modules associated with high metabolic independence in this environment. 

We've essentially confirmed [our previous observations](https://doi.org/10.1186/s13059-023-02924-x), but with a much more extensive dataset of publicly-available gut metagenomes than was used in that study. However, will these observations hold up at the genome level? That is going to be the topic of the next section.


## Analyzing a dataset of gut microbial genomes from the GTDB

We wanted to confirm our results at the genome level. The ideal way to do this would be to carefully bin metagenome-assembled genomes (MAGs) from each metagenome and individually analyze their metabolic capacity. However, binning those MAGs would take a really long time, and even if we automated the process, [automatic binning is difficult and not always conclusive](https://merenlab.org/2020/01/02/visualizing-metagenomic-bins/). Instead, we decided to leverage a high-quality set of reference genomes from the [Genome Taxonomy Database (GTDB)](https://gtdb.ecogenomic.org/) ([Parks et al. 2021](https://doi.org/10.1093/nar/gkab776)).

First, we determined which genomes represent typical gut microbes by running read recruitment and analyzing the resulting coverage information. Then, we analyzed the metabolic potential of each genome by calculating the stepwise completeness of each KEGG module. We used the completeness scores of our 33 IBD-enriched pathways to determine whether each genome represented a microbe with high metabolic independence (HMI) - that is, high average completenesss of all these pathways - or low metabolic independence (LMI). Finally, we used read recruitment results from the gut metagenome dataset to analyze the distribution of each group of genomes across healthy individuals and individuals with IBD.

### Genome processing: the anvi'o contigs workflow

You might remember the [contigs workflow](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#contigs-workflow) from the metagenome processing section above. This workflow, implemented using [snakemake](https://snakemake.readthedocs.io/en/stable/), can take a large number of genome sequences and convert each one into a contigs database. It can also run gene annotation from a variety of functional databases on these genomes.

We used the contigs workflow on all of the representative genomes for species clusters in GTDB release 95.0. The two most critical steps for our downstream analyses were:

* the identification of single-copy core genes (SCGs) in each genome using `anvi-run-hmms`
* annotation of KEGG KOfams using `anvi-run-kegg-kofams` (with our `v2020-12-23` KEGG snapshot)

If you want to see what other programs we ran, you can check the workflow configuration file at `MISC/GTDB_config.json`. 

And if you want to run this for yourself, then you need to 1) download all representative genome sequences for GTDB release 95, 2) make a [fasta.txt](https://anvio.org/help/main/artifacts/fasta-txt/) file that includes the path to each genome, 3) update the config file for your computer system (i.e., changing the number of threads for each rule and possibly adding the path to the right version of KEGG data with `--kegg-data-dir`), and 4) running the following command (possibly modified to work with your HPC's scheduler):

```bash
anvi-run-workflow -w contigs -c MISC/GTDB_config.json
```

Just please keep in mind that you will need a large amount of computational resources. The workflow took about 2 months to finish on our HPC (using 80 cores). Just the zipped fasta files for each of the ~31k representative genomes takes up about 33 GB of storage space, and that does not even consider the storage required for auxiliary files and for the files generated during processing.

### Using the EcoPhylo workflow for quick identification of relevant gut microbes

Once we had all of the GTDB representative genomes annotated, we had to figure out which of those species were relevant to the gut environment. We started by restricting the dataset to genomes assigned to the 3 most common phyla found in the human gut: Firmicutes, Proteobacteria, and Bacteroidetes (`Bacteroidota` in GTDB taxonomy). We simply read the metadata file describing the taxonomy of each representative genome and filtered for those 3 phyla:

```bash
for phylum in Firmicutes Proteobacteria Bacteroidota; do \
  grep p__${phylum} bac120_taxonomy_r95.tsv | cut -f 1 | cut -d '_' -f 2,3
done
```

In the above loop, the `bac120_taxonomy_r95.tsv` file is a metadata file that can be downloaded from [the GTDB FTP site here](https://data.gtdb.ecogenomic.org/releases/release95/95.0/). The `grep` command finds the lines with matching phyla, and the `cut` commands extract the genome accession from the line.

That gives us 19,226 genomes to work with, which is a lot. We wanted to identify the gut microbes within this set by mapping sequencing reads from healthy human gut metagenomes - not the ones used in our study, but an external dataset of 150 gut metagenomes from the [Human Microbiome Project (HMP)](https://doi.org/10.1038/nature11209). However, with so many genomes to map against, a read recruitment workflow to the full genome sequences would have taken years to finish (even on our fancy HPC cluster). So instead we leveraged [Matt Schechter's](https://anvio.org/people/mschecht/) [EcoPhylo workflow](https://anvio.org/help/main/artifacts/ecophylo-workflow/) to do a much faster and less computationally-intensive analysis of the distribution of these genomes in the HMP samples. This workflow takes one gene of interest from each genome and each metagenome, clusters them and picks representative sequences with `mmseqs2` ([Steinegger and Söding 2017](https://doi.org/10.1038/nbt.3988)), and uses the representative sequences to rapidly summarize the distribution of each cluster across the metagenomic samples.

#### Downloading the HMP gut metagenomes

We used 150 gut metagenome samples from healthy people that were published in [this paper for the HMP](https://doi.org/10.1038/nature11209). We downloaded these samples from [the HMP data website](https://www.hmpdacc.org/hmp/resources/data_browser.php). The file at `MISC/HMP_metagenomes.txt` gives the accession numbers of the samples that we used, in case you want to download them for yourself.

Before using them for analysis, we ran quality filtering on these samples using the illumina-utils program `iu-filter-quality-minoche ` and made single assemblies out of 100 of the samples (one per individual, since some people provided multiple stool samples in the HMP study). We used the anvi'o metagenomics workflow for this, as described in the 'Metagenome Processing' section above.

#### Picking a gene to use for EcoPhylo

We wanted to use a single-copy core gene for the EcoPhylo workflow because every genome should have one copy of the gene and because clustering SCG sequences can roughly resolve species-level differences (so not too many irrelevant genomes should end up in the same cluster as genomes from gut microbes). To select our gene of interest, we picked the SCG that was most frequently found across all 19,226 of the GTDB genomes and all 100 of the HMP metagenome assemblies.

First, we obtained a matrix of SCG frequencies:
```bash
anvi-estimate-scg-taxonomy --metagenomes MISC/GTDB_genomes_and_metagenomes.txt \
    --report-scg-frequencies OUTPUT/GTDB_SCG_MATRIX.txt \
    -O GTDB_SCG
```

In the command above, the `GTDB_genomes_and_metagenomes.txt` is a [file describing the paths](https://anvio.org/help/main/artifacts/metagenomes/) to each of the GTDB genomes and HMP metagenomes. You will find it in the `MISC` directory of the datapack.

Then, we used the following R code to list the total number of hits to each SCG within this dataset:

```r
library(tidyverse)

SCG_frequencies <- read_tsv("OUTPUT/GTDB_SCG_MATRIX.txt")

# look in all metagenomes and genomes
SCG_frequencies %>%
  pivot_longer(cols = starts_with("Ribosomal")) %>% 
  dplyr::rename(SCG = "name") %>%
  group_by(SCG) %>%
  summarize(total = sum(value)) %>%
  arrange(desc(total))

# look in only the genomes
SCG_frequencies %>%
filter(!grepl("USA", genome)) %>%
pivot_longer(cols = starts_with("Ribosomal")) %>% 
dplyr::rename(SCG = "name") %>%
group_by(SCG) %>%
summarize(total = sum(value)) %>%
arrange(desc(total))
```

The one that was most frequently identified across all the data (genomes and metagenome assemblies) was Ribosomal Protein S2, but the one that was most frequently present in the GTDB genomes was Ribosomal Protein S6 (`Ribosomal_S6`). We decided to use `Ribosomal_S6` for the EcoPhylo workflow so that we could include as many genomes as possible in our analysis.

#### Running the EcoPhylo workflow

To run the workflow, you need a file describing which gene to use, a file describing the paths to all of the GTDB genomes, a file describing the paths to all of the metagenome assemblies, and a file describing the paths to all of the HMP samples (quality-filtered sequencing reads) to be used for read recruitment. We show how to generate the first three files:

```bash
echo -e "name\tsource\tpath" > MISC/hmm_list.txt
echo -e "Ribosomal_S6\tBacteria_71\tINTERNAL" >> MISC/hmm_list.txt

grep -v USA MISC/GTDB_genomes_and_metagenomes.txt > MISC/GTDB_external-genomes.txt

echo -e "name\tcontigs_db_path" > MISC/HMP_metagenome_assemblies.txt
grep USA MISC/GTDB_genomes_and_metagenomes.txt >> MISC/HMP_metagenome_assemblies.txt
```

To generate the last file, make a [samples-txt](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#samplestxt) with the paths to the HMP samples, wherever they are on your computer (we called this file `GTDB_HMP_samples.txt`). You will probably also have to change the paths to the GTDB genomes and HMP assemblies in `GTDB_genomes_and_metagenomes.txt` and all derivative files.

The configuration file for the workflow can be found at `MISC/ecophylo_config.json`. It contains the names of each of the input files from the `MISC/` folder. You can run the workflow using the following command:

```bash
anvi-run-workflow -w ecophylo -c MISC/ecophylo_config.json
```

As always, please remember to change the thread counts appropriately for your system (we made sure the config files in the datapack doesn't assign more than 5 threads to any rule so it doesn't inadvertently overwhelm your computer, but you definitely want to increase those if you can). On our system, the workflow took a few months to finish (with some stops and restarts for troubleshooting), which is still much faster than read recruitment to the full genome sequences would have taken. :)

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
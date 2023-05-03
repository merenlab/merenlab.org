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
* [doi:10.5281/zenodo.7883421](https://doi.org/10.5281/zenodo.7883421): Contigs databases for 338 GTDB genomes.


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

## Downloading the data for this reproducible workflow

We've pre-packaged a lot of the data and scripts that you need for this workflow into a DATAPACK. To download it, run the following in your terminal:

```bash
# download the datapack
FIXME download
# extract it 
tar -xvzf VESELI_2023_DATAPACK.tar.gz

# move inside the datapack
cd VESELI_2023_DATAPACK/
```

We suggest working from within this datapack. Most of the BASH commands in subsequent sections will assume that your starting working directory is the uppermost level of the datapack directory structure. To keep things organized, we'll have you generate a folder to work in for each major section of this workflow.

The following two subsections describe how to download the [contigs databases](https://anvio.org/help/main/artifacts/contigs-db/) that we generated for the metagenomes and genomes analyzed in the paper. They are pretty big, so keep their storage requirements in mind before you download them. You can still follow along with some parts of the workflow even if you don't download them.

### Downloading our contigs databases for metagenome assemblies

We provide access to the main set of metagenome assemblies that we analyzed in the paper [here](FIXME LINK). This is a rather large dataset, so we separated it from the rest of the datapack to avoid overburdening your storage system unnecessarily. If you want to access these assemblies, you can download them into your datapack by running the following commands:

{:.warning}
The size of this metagenome dataset is **96 GB** (the archive alone is ~33 GB). Please make sure you have enough space to download it!

```bash
# download the metagenome data archive
FIXME download
# extract the metagenome data
tar -xvzf VESELI_ET_AL_METAGENOME_CONTIGS_DBS.tar.gz

# generate a table of sample names and paths
echo -e "name\tcontigs_db_path" > METAGENOME_EXTERNAL_GENOMES.txt
while read db; do \
  sample=$(echo $db | sed 's/-contigs.db//'); \
  path=$(ls -d $PWD/VESELI_ET_AL_METAGENOME_CONTIGS_DBS/${db}); \
  echo -e "$sample\t$path" >> METAGENOME_EXTERNAL_GENOMES.txt; \
done < <(ls VESELI_ET_AL_METAGENOME_CONTIGS_DBS/)
```

Once you run the above code, you should see in the datapack a folder called `VESELI_ET_AL_METAGENOME_CONTIGS_DBS` that contains 408 database files, and a file called `METAGENOME_EXTERNAL_GENOMES.txt` that describes the name and absolute path to each sample on your computer. If everything looks good, you can delete the archive to get back some storage space:

```
# clean up the archive
rm VESELI_ET_AL_METAGENOME_CONTIGS_DBS.tar.gz
```

### Downloading our contigs databases for GTDB genomes

Likewise, there is [yet another link](FIXME link) to download the set of GTDB genomes that we analyzed in the paper. It takes up only **2 GB** of space. You can download it by running the following:

```bash
# download the genome data archive
wget https://zenodo.org/record/7883421/files/VESELI_ET_AL_GENOME_CONTIGS_DBS.tar.gz
# extract the genome data
tar -xvzf VESELI_ET_AL_GENOME_CONTIGS_DBS.tar.gz

# generate a table of genome names and paths
echo -e "name\tcontigs_db_path" > GTDB_EXTERNAL_GENOMES.txt
while read db; do \
  acc=$(echo $db | cut -d '.' -f 1); \
  ver=$(echo $db | sed 's/-contigs.db//' | cut -d '.' -f 2); \
  genome="${acc}_${ver}"; \
  path=$(ls -d $PWD/SUBSET_GTDB_CONTIGS_DBS/${db}); \
  echo -e "$genome\t$path" >> GTDB_EXTERNAL_GENOMES.txt; \
done < <(ls SUBSET_GTDB_CONTIGS_DBS/)
```

Once you run the above code, you should see in the datapack a folder called `SUBSET_GTDB_CONTIGS_DBS` that contains 338 database files, and a file called `GTDB_EXTERNAL_GENOMES.txt` that describes the name and absolute path to each genome's database on your computer. If everything looks good, you can delete the archive:

```
# clean up the archive
rm VESELI_ET_AL_GENOME_CONTIGS_DBS.tar.gz
```

### Datapack description

Here is a quick overview of the datapack structure:

```
VESELI_2023_DATAPACK/
  |
  |- VESELI_ET_AL_METAGENOME_CONTIGS_DBS  ## this holds the contigs databases we generated for the subset of samples we analyzed
  |- SUBSET_GTDB_CONTIGS_DBS              ## this holds the contigs databases we generated for the subset of GTDB genomes we analyzed
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

Whenever KEGG-related programs are used, you can make them use the appropriate KEGG version by adding `--kegg-data-dir KEGG_2020-12-23` (replacing that path with wherever you decided to store the KEGG data on your computer). In the workflows and commands on this page, we'll assume you followed the setup command exactly as written and add the directory name `KEGG_2020-12-23`.


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

If you want to download all of the samples we used in this work (keeping in mind that the storage requirements for almost 3,000 metagenomes will be huge), here we show you how to download all of the samples with one script, into the same folder (for better organization and easier compatibility with the later sections of this workflow). But it truly doesn't matter how you decide to organize the samples as long as you can keep track of the paths to each sample on your computer. So if you want to do it differently, go for it (and feel free to reach out to us for help if you need it). 

For convenience, we've provided a plain-text version of Table 1c in our DATAPACK, which can be accessed at the path `TABLES/00_ALL_SAMPLES_INFO.txt`. The last column of that file provides the SRA accessions that can be used for downloading each sample.

There is one exception to this strategy, and that is the study by [Quince et al. 2015](https://doi.org/10.1038/ajg.2015.357). There are no deposited sequences under [the NCBI BioProject for this study](https://www.ncbi.nlm.nih.gov/bioproject/270985). The SRA accession column for these samples in `TABLES/00_ALL_SAMPLES_INFO.txt` contains NaN values, and these rows should be skipped when using `fasterq-dump` to download samples. We accessed these metagenomes directly from the study authors.

Let's make a folder to work in for this section of the workflow:

```bash
mkdir 01_METAGENOME_DOWNLOAD
cd 01_METAGENOME_DOWNLOAD/
```

And, let's make a folder in which you can download the samples:

```bash
mkdir 00_FASTQ_FILES
```

Now, we'll extract a list of the SRA accessions (for all except the Quince et al. samples) from the provided table.

```bash
grep -v Quince_2015 ../TABLES/00_ALL_SAMPLES_INFO.txt | cut -f 9 | tail -n+2 > sra_accessions_to_download.txt
```

There should be 2,824 accessions in that file. In the `SCRIPTS` folder of your datapack, there is a script called `download_sra.sh` that will download each sample with `prefetch`, unpack it with `fasterq-dump` into the `00_FASTQ_FILES` folder that you just created, gzip the resulting FASTQ files, and then delete the intermediate files. You can run it using the following command (we recommend running this on an HPC cluster, with plenty of threads):

```bash
../SCRIPTS/download_sra.sh sra_accessions_to_download.txt
```

Note that the script has no error-checking built in, so if any of the samples fail to download, you will have to manage those yourself. But once it is done, you should have the FASTQ files for all 2,824 samples in the `00_FASTQ_FILES` directory. You can add in the remaining Quince et al. samples however you acquire those (we'd be happy to pass them along).

Once you have all the samples ready, you can generate a file called `METAGENOME_SAMPLES.txt` describing the path to each sample like this:

```bash
echo -e "sample\tr1\tr2" > METAGENOME_SAMPLES.txt
while read sample diag study doi grp num1 num2 nump sra; do \
  r1=$(ls -d $PWD/00_FASTQ_FILES/${sra}_1.fastq.gz) ;\
  r2=$(ls -d $PWD/00_FASTQ_FILES/${sra}_2.fastq.gz) ;\
  echo -e "$sample\t$r1\t$r2" >> METAGENOME_SAMPLES.txt ;\
done < <(tail -n+2 ../TABLES/00_ALL_SAMPLES_INFO.txt | grep -v Vineis)
while read sample diag study doi grp num1 num2 nump sra; do \
  r1=$(ls -d $PWD/00_FASTQ_FILES/${sra}.fastq.gz) ;\
  r2="NaN" ;\
  echo -e "$sample\t$r1\t$r2" >> METAGENOME_SAMPLES.txt ;\
done < <(grep Vineis ../TABLES/00_ALL_SAMPLES_INFO.txt)
```

We will use this file later to make input files for our workflows.

Now that all the sequences are downloaded, let's go back to the top-level datapack directory in preparation for the next section.

```bash
cd ..
```

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

Let's run this workflow in a new folder:

```bash
mkdir 02_METAGENOME_PROCESSING
cd 02_METAGENOME_PROCESSING/
```

We provide an example configuration file (`MISC/metagenomes_config.json`) in the DATAPACK that can be used for reproducing our assemblies. Copy that file over to your current working directory:

```bash
cp ../MISC/metagenomes_config.json .
```

To run the workflow, you simply create a 3-column `samples.txt` file containing the sample name, path to the R1 file, and path to the R2 file for each sample that you downloaded. An example file is described in our [workflow tutorial](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#samplestxt). In fact, you can derive this file from `TABLES/00_ALL_SAMPLES_INFO.txt` (assuming you have either followed our file naming/organization recommendations or updated those paths to reflect the names/organization you decided upon):

```bash
grep -v NaN ../01_METAGENOME_DOWNLOAD/METAGENOME_SAMPLES.txt > samples.txt
```

Then, you can start the workflow with the following command (hopefully adapted for use on a high-performance computing cluster):

```bash
anvi-run-workflow -w metagenomics -c metagenomes_config.json
```

A few notes:
* We renamed the samples from each study to incorporate information such as country of origin (for healthy samples) or host diagnosis (for IBD samples) for better readability and downstream sorting. To match our sample names, your `samples.txt` file should use the same sample names that are described in Supplementary Table 1c (or the first column of `TABLES/00_ALL_SAMPLES_INFO.txt` in the DATAPACK). If you generated the `samples.txt` from that file, you are good to go.
* We used the default snapshot of KEGG data associated with anvi'o v7.1-dev, which can be downloaded onto your computer by running `anvi-setup-kegg-kofams --kegg-snapshot  v2020-12-23`, as described earlier. To exactly replicate the results of this study, the metagenome samples need to be annotated with this KEGG version by changing the `--kegg-data-dir` parameter (in the `anvi_run_kegg_kofams` rule of the config file) to point to this snapshot wherever it located on your computer
* The number of threads used for each rule is set in the config file. We conservatively set this number in the example `MISC/metagenomes_config.json` to be 1 for all rules, but you will certainly want to adjust these to take advantage of the resources of your particular system.

The steps in the workflow described above apply to all of the metagenome samples except for those from [Vineis et al. 2016](https://doi.org/10.1128/mBio.01713-16), which had to be processed differently since the downloaded samples contain merged reads (rather than paired-end reads described in R1 and R2 files, as in the other samples). Since the metagenomics workflow currently works only on paired-end reads, we had to run the assemblies manually. Aside from the lack of workflow, there are only two major differences in the processing of the 96 Vineis et al. samples:

* No additional quality-filtering was run on the downloaded samples, because the merging of the sequencing reads as described in the paper's methods section already included a quality-filtering step
* Single assembly of the merged reads was done with [MEGAHIT](https://academic.oup.com/bioinformatics/article/31/10/1674/177884), using all default parameters except for a minimum contig length of 1000

We wrote a loop to run an individual assembly on each sample from [Vineis et al. 2016](https://doi.org/10.1128/mBio.01713-16). This loop makes use of the sample names and paths as established in the `01_METAGENOME_DOWNLOAD/METAGENOME_SAMPLES.txt` file:

```bash
while read name path; do \
  megahit -r $path  \
      --min-contig-len 1000 \
      -t 7 \
      -o ${name}_TMP; \
done < <(grep NaN ../01_METAGENOME_DOWNLOAD/METAGENOME_SAMPLES.txt | cut -f 1-2)
```

Once the assemblies were done, we extracted the final assembly files from each output directory, renamed them with the sample name, and put them all in one folder:

```bash
mkdir -p VINEIS_ASSEMBLIES
while read name path; do \
  mv ${name}_TMP/final.contigs.fa VINEIS_ASSEMBLIES/${name}.fasta
done < <(grep NaN ../01_METAGENOME_DOWNLOAD/METAGENOME_SAMPLES.txt | cut -f 1-2)
rm -r *TMP/
```

Then, we were able to leverage the anvi'o contigs workflow to generate the contigs databases for each assembly and run the annotation steps. We've provided the relevant configuration file for this workflow (`MISC/vineis_config.json`) as well as the input file that lists the path to each assembly (`MISC/vineis_fasta.txt`) in the DATAPACK, and this is how you could run it for yourself:

```bash
cp ../MISC/vineis_* .
anvi-run-workflow -w contigs -c vineis_config.json
```

The same notes about setting the KEGG data version and the number of threads that we detailed for the metagenomics workflow apply to this workflow as well.

### A final note on data organization

If you've elected to reproduce the metagenome download and processing described in this section, you will have ended up with a lot of samples, assemblies, and contigs databases on your computer system. You will need to access the databases several times later in the workflow, so you should generate an [external genomes file](https://anvio.org/help/main/artifacts/external-genomes/) containing the paths to each database:

```bash
# generate a table of sample names and paths
echo -e "name\tcontigs_db_path" > ALL_METAGENOME_DBS.txt
while read db; do \
  sample=$(echo $db | sed 's/-contigs.db//'); \
  path=$(ls -d $PWD/03_CONTIGS/${db}); \
  echo -e "$sample\t$path" >> ALL_METAGENOME_DBS.txt; \
done < <(ls 03_CONTIGS/)
```

That file should have a line for all 2,893 samples in it (plus a line for the header), if everything went well in the workflows above.

If you decided to follow a different file naming and organization strategy than the one used above, that is fine. Where the files are on your computer does not matter for following the remainder of this workflow, as long as you prepare a file for yourself that describes the correct paths to each sample's read files and contigs database.

We'll stay in this directory for the next section.

## Selecting our final dataset of metagenomes

{:.notice}
You should be inside the `02_METAGENOME_PROCESSING` directory for this section of the workflow.

Before we began analyzing metabolism within the large dataset we compiled, we discovered that we had to reduce our sample set to ensure accurate calculations. One of the critical steps in comparing community-level copy numbers of metabolic pathways between microbial communities of differing richness is **normalization of these data by the size of the community**. Pathway copy numbers will naturally tend to be higher in metagenome assemblies that describe larger communities, so comparing these values is not very meaningful when community sizes are vastly different. The gut microbiomes of people with IBD tend to harbor much less diversity than healthy gut microbiomes, so this is certainly a problem in our case. We therefore came up with a strategy of normalizing the pathway copy numbers calculated for a given sample with the number of microbial populations represented within that metagenome assembly.

However, this normalization only works if we can accurately estimate that number of populations in each sample. Yet this is not the case for samples of low sequencing depth, as you will see later in this section - to a certain extent, the estimated number of populations is correlated with sample sequencing depth. We interpreted this to mean that low-depth samples fail to capture enough sequences from low-abundance populations, thereby skewing the population estimates in these assemblies. Therefore, we tried to mitigate the issue by filtering out low-depth samples and keeping only samples with high enough sequencing depth for our downstream analyses.

This section will cover how we estimate the number of populations in each metagenome assembly using single-copy core genes, our analysis of its relationship with sequencing depth, and the removal of low-depth samples to establish our final set of samples for analysis.

### Estimating number of populations per sample

To estimate how many microbial populations are represented in a metagenome assembly, we can rely on the fact that all microbial genomes (with few exceptions) contain exactly one copy of each gene in a special set of essential genes called single-copy core genes (SCGs). These include ribosomal proteins and other housekeeping genes. Anvi'o ships with a few generic sets of SCGs that is each specific to a domain of microbial life (Bacteria, Archaea, and Protista) and these genes are annotated using the program `anvi-run-hmms` (which, you might recall, we ran earlier as part of our metagenome processing workflows).

Since we expect to find one copy of each SCG in each microbial population, we can count the total number of copies of an SCG in a metagenome assembly and use that as the number of populations. However, using just one SCG for estimation would be error-prone due to missing SCG annotations in incomplete data (or the occasional duplication within one genome). Instead, we can use all of the SCGs for a given domain to make the estimate more robust to noise. The **mode of the number of SCGs** in the assembly gives us our estimate of the number of populations, in this case. The sketch below illustrates this process - each SCG is annotated in the metagenome assembly as indicated by the black boxes on the highlighted sequences on the bottom (shown for the first 3 SCGs only), these annotations are tallied (top histogram), and then the mode of the counts is computed:

[![Estimating the number of microbial populations using the mode of single-copy core gene annotations](images/estimate_pops_with_SCGs.png)](images/estimate_pops_with_SCGs.png){:.center-img .width-50}

The anvi'o codebase includes a `NumGenomesEstimator` class that does exactly this: takes the mode of the number of copies of the SCGs for a particular domain, giving you an estimate of the number of bacteria, archaea, and protists in a given metagenome assembly. Those values can be added together to obtain the total number of microbial populations in the sample (for gut metagenomes, usually only bacterial species are found).

We wrote a script that runs this estimation on each of our 2,893 samples. You can run it using the following command, and it will produce a table called (`NUM_POPULATIONS_ALL_SAMPLES.txt`). Provide the file indicating the location of all databases to the script:

```bash
python ../SCRIPTS/estimate_num_genomes.py ALL_METAGENOME_DBS.txt
```

If you didn't download all the metagenome samples, or you don't feel like running this part, you can access the estimates in the `num_populations` column of the `TABLES/00_ALL_SAMPLES_INFO.txt` file.

### Determining sequencing depth

We used a BASH loop to count the number of reads in each metagenome sample. It counts the number of lines in each (gzipped) FASTQ file and divides by 4 (the number of lines per read) to get the number of reads, which is then stored in a tab-delimited file. R1 files and R2 are counted separately. For the Vineis et al. samples, the count of merged reads per sample is stored in the R1 column, and the R2 column is 0 (because there is only one FASTQ file for each of these samples).

```bash
echo -e "name\tnum_reads_r1\tnum_reads_r2" > METAGENOME_NUM_READS.txt
while read name r1 r2; do \
  numr1=$(echo $(zcat $r1 | wc -l) / 4 | bc); \
  if [ "$r2" = "NaN" ]; then \
    numr2=0; \
  else numr2=$(echo $(zcat $r2 | wc -l) / 4 | bc); \
  fi;
  echo -e "$name\t$numr1\t$numr2" >> METAGENOME_NUM_READS.txt; \
done < <(tail -n+2 ../01_METAGENOME_DOWNLOAD/METAGENOME_SAMPLES.txt)
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
# first, go back to the parent directory
cd ..

python SCRIPTS/subset_metagenome_samples.py
```

### Final set of samples and their contigs DBs

In the remainder of the analyses described in our manuscript, we utilized the subset of 408 samples with high sequencing depth described in `TABLES/00_SUBSET_SAMPLES_INFO.txt`. Since 408 is a much more reasonable number than 2,893 we have provided the contigs databases for our metagenome assemblies of these samples in [this datapack](FIXME LINK). If you elected to download this datapack via the instructions at the start of this workflow, you will find the assemblies in the directory called `VESELI_ET_AL_METAGENOME_CONTIGS_DBS/`, and the table `METAGENOME_EXTERNAL_GENOMES.txt` contains the paths to these files. So, even if you didn't download and process all of the metagenome samples as described in the first section of this workflow, you can still continue with the subsequent sections.


## Metabolism analyses (for metagenomes)

In this section, we will cover the analyses used to determine the set of metabolic pathways that are enriched in the IBD sample group. First, we calculated the copy numbers of all pathways in the KEGG modules database, in each metagenome assembly. Those numbers are made from the combined genes from all microbial populations represented in a given metagenome and thus quantify the community-level metabolic potential. To make them comparable across different gut communities of varying diversity, we normalized each copy number by the estimated number of populations in the same sample to obtain a 'per-population copy number', or PPCN. We then ran a statistical test on each module to determine which pathways were **the most different between the sample groups**, specifically looking for the ones that had **higher PPCN in IBD**.

In this section, we will work from the top-level directory of the datapack, but we'll make a folder in which we can generate some output:

```bash
mkdir 03_METABOLISM_OUTPUT
```

### Metabolism estimation

We used the program `anvi-estimate-metabolism` to compute copy numbers of KEGG modules in each sample. You can find details about that program and its calculation strategies on [this page](https://anvio.org/help/main/programs/anvi-estimate-metabolism/). To ensure that all annotated KOfams in the metagenome contributed to the calculations, we ran the program in 'genome mode' on each sample, and to do this in a high-throughput manner, we provided the program with an [external genomes file](https://anvio.org/help/main/artifacts/external-genomes/) containing the paths to all the samples' contigs databases at once. Here is the code to both generate that input file (from the table containing the sample information, including paths) and run the metabolism estimation code.

```bash
# run metabolism estimation
anvi-estimate-metabolism -e METAGENOME_EXTERNAL_GENOMES.txt \
                    -O 03_METABOLISM_OUTPUT/METAGENOME_METABOLISM \
                    --kegg-data-dir KEGG_2020-12-23 \
                    --add-copy-number \
                    --matrix-format
```

Note the use of `--kegg-data-dir KEGG_2020-12-23`, which points the program to the version of KEGG data that we used for our analysis (which you should have downloaded when going through the 'Computational environment details' section). If the path to this data on your computer is different, you will have to change that in the command before you run it.

Since we used the `--matrix-format` flag, the output of `anvi-estimate-metabolism` will be a set of matrices, each containing a different statistic summarized across all of the metagenomes and all of the pathways. The one that we want for downstream analysis is the one containing stepwise copy numbers, which can be found at the path `03_METABOLISM_OUTPUT/METAGENOME_METABOLISM-module_stepwise_copy_number-MATRIX.txt`. 

We've included our version of this output file in the DATAPACK. You will find it at `TABLES/METAGENOME_METABOLISM-module_stepwise_copy_number-MATRIX.txt` (and this table is referenced in the scripts below that depend on this data, just in case you don't run this step).

### Normalization of pathway copy numbers to PPCN

We just calculated the pathway copy numbers, and we estimated the number of microbial populations in each metagenome assembly in the previous section. To normalize the former, we divide by the latter. It is that simple :) 

The R script at `SCRIPTS/module_stats_and_medians.R` contains the code for computing PPCN values. It gets the module copy numbers and per-sample population estimates from the tables in the `TABLES/` folder, and generates a long-format data table at `03_METABOLISM_OUTPUT/ALL_MODULES_NORMALIZED_DATA.txt` that includes the PPCN values in the `PPCN` column.

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

Later on, we will also make use of the median per-population copy number for each module in each sample group (for visualizing the data as boxplots). The script also generates these median values into the file `03_METABOLISM_OUTPUT/ALL_MODULES_MEDIAN_PPCN.txt`. If you plan to reproduce our figures, you should run that section of code, too (it is clearly marked in the comments).

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

If you keep following the code in the script, you will see that it will generate a file with the statistical test results and a variety of other data for each module, including the enrichment status, at `03_METABOLISM_OUTPUT/ALL_MODULES_MEDIANS_AND_STATS.txt`. It also prints a subset of this information for the IBD-enriched modules to the file at `03_METABOLISM_OUTPUT/IBD_ENRICHED_MODULES.txt`.

One more note - you may notice that we exclude one module from the automatically-generated IBD-enriched set:

```r
## FINAL LIST OF IBD-ENRICHED MODULES
# exclude M00006 because it is the first part of M00004
ibd_enriched = ordered_modules %>% filter(module != "M00006")
```

This is because the module in question, [M00006](https://www.genome.jp/entry/M00006), is the oxidative phase of the pentose phosphate cycle and overlaps completely with the module describing the entirety of the pentose phosphate cycle, [M00004](https://www.genome.jp/entry/M00004). We also considered removing module [M00007](https://www.genome.jp/entry/M00007) because it is the non-oxidative phase of this cycle; however, it has some enzymatic differences with module M00004, so we left it in.

### Computing proportion of shared enzymes in IBD-enriched pathways

You may have noticed a few statements in our paper about how closely intertwined the IBD-enriched metabolic pathways are. They share a high proportion of enzymes and compounds. To get the data backing up that statement (which can be found in Supplementary Table 2 as well), you can run the following script:

```bash
python SCRIPTS/get_num_shared_enzymes_compounds.py KEGG_2020-12-23/MODULES.db
```

Note that it expects as input the path to the modules database (here `KEGG_2020-12-23/MODULES.db`). If you set up that data elsewhere, you'll need to change the path. The script will output a table at `03_METABOLISM_OUTPUT/IBD_ENRICHED_SHARED_ENZYMES_COMPOUNDS.txt`.

### Generating Figure 2

Figure 2 includes several plots of both unnormalized copy numbers and PPCN values. We won't copy the code to generate these plots here, but just wanted to remind you that you can find it in the script `SCRIPTS/plot_figures.R`. The section of code for each plot was designed to be an independent as possible from the other sections of code, so that you need only run that section to get what you want. Most of them do require the packages, paths, and color variables set up at the top of the script, however.

We used Inkscape to polish up the resulting plots into a publication-ready figure. In particular, we adjusted the background color opacity for the heatmaps in panels B and E to make the signal more visible, which is why you may notice a difference between the generated plots and the plots in the manuscript.

### Summary of metagenome results

So far, we've analyzed gut metagenomes, using per-population copy number as our metric for the typical metabolic capacity of a microbial population in the community described by each metagenome. We obtained a list of metabolic pathways that are enriched in the IBD sample group, that largely overlap with the modules associated with metabolic independence from [our previous study](https://doi.org/10.1186/s13059-023-02924-x), and that are mostly biosynthesis pathways for important cellular metabolites. These 33 IBD-enriched pathways could be important to microbial survival in the depleted communities of the IBD gut environment, and represent a refined list of modules associated with high metabolic independence in this environment. 

We've essentially confirmed [our previous observations](https://doi.org/10.1186/s13059-023-02924-x), but with a much more extensive dataset of publicly-available gut metagenomes than was used in that study. However, will these observations hold up at the genome level? That is going to be the topic of the next section.


## Obtaining a dataset of gut microbial genomes from the GTDB

We wanted to confirm our results at the genome level. The ideal way to do this would be to carefully bin metagenome-assembled genomes (MAGs) from each metagenome and individually analyze their metabolic capacity. However, binning those MAGs would take a really long time, and even if we automated the process, [automatic binning is difficult and not always conclusive](https://merenlab.org/2020/01/02/visualizing-metagenomic-bins/). Instead, we decided to leverage a high-quality set of reference genomes from the [Genome Taxonomy Database (GTDB)](https://gtdb.ecogenomic.org/) ([Parks et al. 2021](https://doi.org/10.1093/nar/gkab776)).

In this section, we show how we determined which genomes represent typical gut microbes by running read recruitment and analyzing the resulting coverage information. Those genomes will be the subject of the analyses described in the subsequent section. The datapack at [https://doi.org/10.5281/zenodo.7883421](https://doi.org/10.5281/zenodo.7883421) provides the contigs databases for the final set of gut genomes selected in this section.

This section includes several workflows, so we will start a new folder to keep everything organized:

```bash
mkdir 04_GTDB_PROCESSING
cd 04_GTDB_PROCESSING/
```

### Genome processing: the anvi'o contigs workflow

You might remember the [contigs workflow](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#contigs-workflow) from the metagenome processing section above. This workflow, implemented using [snakemake](https://snakemake.readthedocs.io/en/stable/), can take a large number of genome sequences and convert each one into a contigs database. It can also run gene annotation from a variety of functional databases on these genomes.

We used the contigs workflow on all of the representative genomes for species clusters in GTDB release 95.0. The two most critical steps for our downstream analyses were:

* the identification of single-copy core genes (SCGs) in each genome using `anvi-run-hmms`
* annotation of KEGG KOfams using `anvi-run-kegg-kofams` (with our `v2020-12-23` KEGG snapshot)

If you want to see what other programs we ran, you can check the workflow configuration file at `MISC/GTDB_contigs_config.json`. 

And if you want to run this for yourself, then you need to 1) download all representative genome sequences for GTDB release 95, 2) make a [fasta.txt](https://anvio.org/help/main/artifacts/fasta-txt/) file that includes the path to each genome, 3) update the config file for your computer system (i.e., changing the number of threads for each rule and possibly adding the path to the right version of KEGG data with `--kegg-data-dir`), and 4) running the following (hopefully, you will modify the workflow command to work with your HPC's scheduler):

```bash
# first make a folder in which to run the workflow
mkdir GTDB_CONTIGS_WORKFLOW
cd GTDB_CONTIGS_WORKFLOW/

# copy over the config file
cp ../../MISC/GTDB_contigs_config.json .

# workflow command
anvi-run-workflow -w contigs -c GTDB_contigs_config.json
```

Just please keep in mind that you will need a large amount of computational resources. The workflow took about 2 months to finish on our HPC (using 80 cores). Just the zipped fasta files for each of the ~31k representative genomes takes up about 33 GB of storage space, and that does not even consider the storage required for auxiliary files and for the files generated during processing.

We'll start the next subsection from the parent directory (`04_GTDB_PROCESSING`):

```bash
cd ..
```

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

We used 150 gut metagenome samples from healthy people that were published in [this paper for the HMP](https://doi.org/10.1038/nature11209). We downloaded these samples from [the HMP data website](https://www.hmpdacc.org/hmp/resources/data_browser.php). The file at `MISC/HMP_metagenomes.txt` gives the accession numbers of the samples that we used (in the 3rd column), in case you want to download them for yourself (which you can do by re-using the `download_sra.sh` script).

Before using them for analysis, we ran quality filtering on these samples using the illumina-utils program `iu-filter-quality-minoche ` and made single assemblies out of 100 of the samples (one per individual, since some people provided multiple stool samples in the HMP study). We used the anvi'o metagenomics workflow for this, as described in the 'Metagenome Processing' section above.

#### Picking a gene to use for EcoPhylo

We wanted to use a single-copy core gene, and specifically a ribosomal protein, for the EcoPhylo workflow because every genome should have one copy of the gene and because clustering ribosomal protein sequences can often roughly resolve species-level differences (they are often used for phylogenomics for this reason). To select our gene of interest, we picked the SCG that was most frequently found across all 19,226 of the GTDB genomes and all 100 of the HMP metagenome assemblies.

First, we obtained a matrix of SCG frequencies:
```bash
anvi-estimate-scg-taxonomy --metagenomes ../MISC/GTDB_genomes_and_metagenomes.txt \
    --report-scg-frequencies GTDB_SCG_MATRIX.txt \
    -O GTDB_SCG
```

In the command above, the `GTDB_genomes_and_metagenomes.txt` is a [file describing the paths](https://anvio.org/help/main/artifacts/metagenomes/) to each of the GTDB genomes and HMP assemblies. You will find it in the `MISC` directory of the datapack. This file assumes that the GTDB genome databases are located in the output directory of the contigs workflow you ran in the folder `GTDB_CONTIGS_WORKFLOW`, and that the databases for the HMP samples are located in the current directory (so you'll have to change those paths to wherever those assemblies are located on your own computer).

Then, we used the following R code to list the total number of hits to each SCG within this dataset:

```r
library(tidyverse)

SCG_frequencies <- read_tsv("GTDB_SCG_MATRIX.txt")

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
echo -e "name\tsource\tpath" > hmm_list.txt
echo -e "Ribosomal_S6\tBacteria_71\tINTERNAL" >> hmm_list.txt

grep -v USA ../MISC/GTDB_genomes_and_metagenomes.txt > GTDB_external-genomes.txt

echo -e "name\tcontigs_db_path" > HMP_metagenome_assemblies.txt
grep USA ../MISC/GTDB_genomes_and_metagenomes.txt >> HMP_metagenome_assemblies.txt
```

To generate the last file, make a [samples-txt](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#samplestxt) with the paths to the HMP samples, wherever they are on your computer (we called this file `GTDB_HMP_samples.txt`). You will probably also have to change the paths to the HMP assemblies in `GTDB_genomes_and_metagenomes.txt` and all derivative files (if these are not in your current working directory).

The configuration file for the workflow can be found at `MISC/ecophylo_config.json`. It contains the names of each of the input files from the `MISC/` folder. You can run the workflow using the following commands:

```bash
cp ../MISC/ecophylo_config.json .
anvi-run-workflow -w ecophylo -c ecophylo_config.json
```

As always, please remember to change the thread counts appropriately for your system (we made sure the config files in the datapack doesn't assign more than 5 threads to any rule so it doesn't inadvertently overwhelm your computer, but you definitely want to increase those if you can). On our system, the workflow took a few months to finish (with some stops and restarts for troubleshooting), which is still much faster than read recruitment to the full genome sequences would have taken. :)

### Subsetting gut genomes by Ribosomal Protein S6 detection in the HMP metagenomes

Once the EcoPhylo workflow finished (producing output within the `04_GTDB_PROCESSING/ECOPHYLO_WORKFLOW` directory), we had access to the following data: 

* clusters of similar Ribosomal Protein S6 sequences (>94% nucleotide identity) from the GTDB genomes and HMP metagenome assemblies
* read recruitment information from the HMP metagenome samples to each cluster's representative sequence

We considered a genome to be a gut microbe if its `Ribosomal_S6` sequence belonged to a cluster for which the representative had at least 50% 'detection' in at least 10% of the HMP metagenomes. Detection is the proportion of nucleotide positions in the sequence that are covered by at least one sequencing read, and its a value that anvi'o calculates automatically from read mapping data. We extracted a table of detection values (from each sample) for each cluster's representative from the profile database produced by the workflow:

```bash
anvi-export-table PROFILE.db --table detection_splits -o Ribosomal_S6_detection.txt
```

The profile database is not in the datapack (so you cannot run the code above unless you ran the EcoPhylo workflow), but the resulting detection table is, and you can find it at `TABLES/Ribosomal_S6_detection.txt`.

In the datapack, you will also find a file at `TABLES/Ribosomal_S6-mmseqs_NR_cluster.tsv` which describes the clusters generated by `mmseqs2` during the workflow. The first column in that file includes the representative sequences, and the second column includes each genome belonging to the cluster with the corresponding representative.

To extract a list of the genomes that fit our detection criteria in the HMP metagenomes, we used a Python script, which can be found in the datapack at `SCRIPTS/extract_ecophylo_gut_genomes.py`. It uses the two files just mentioned to find all the representative sequences with sufficient detection and extract all relevant genome accessions from the sequence names in their corresponding clusters. You can run it like so:

```bash
python ../SCRIPTS/extract_ecophylo_gut_genomes.py
```

It will produce a list of 836 genomes at `04_GTDB_PROCESSING/gut_genome_list.txt`.

### Read recruitment from our dataset of gut metagenomes to gut microbes

Our next task was to recruit reads from our dataset of 408 deeply-sequenced gut metagenomes to the 836 GTDB genomes. To do this efficiently, we concatenated all of the genome sequences into the same FASTA file and used the anvi'o [metagenomics workflow in 'References Mode'](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#references-mode) to run the read mapping.

Let's make a new folder for the read recruitment workflow.

```bash
mkdir GTDB_MAPPING_WORKFLOW
cd GTDB_MAPPING_WORKFLOW/
```

#### Reformatting the genome FASTAs
When making a FASTA file containing multple genomes, you need each contig sequence to be uniquely labeled, and you want to be able to easily identify which genome each contig sequence belongs to (in case you need to backtrack later). The easiest way to do that is to include the genome accession in each contig header. We acheived this by first running `anvi-script-reformat-fasta` on each of the individual GTDB genome FASTA files to get a version of the file with the genome accession prefixing each contig header.

The loop below creates an [external genomes file](https://anvio.org/help/main/artifacts/external-genomes/) containing the name and path of each genome in our set of 836. To run this yourself, you will likely need to change the paths and/or filenames to the genome FASTA files, depending on where you stored them on your computer. Note that we used the naming convention of `${acc}.${ver}_genomic.fna.gz` (where `${acc}` is the genome accession and `${ver}` is its version number), and that we keep these files gzipped to reduce the storage requirement.

```bash
echo -e "name\tpath" > GTDB_genomes_for_read_recruitment.txt
while read g; do \
  ver=$(echo $g | cut -d '_' -f 3); \
  acc=$(echo $g |cut -d '_' -f 1,2); \
  echo -e "${g}\t${acc}.${ver}_genomic.fna.gz" >> GTDB_genomes_for_read_recruitment.txt; #may need to change path in this line \
done < ../gut_genome_list.txt
```

Then, to reformat the contig headers for each genome, we ran the following loop, which unzips the original genome FASTA file and runs the reformatting on it to produce a new file with the appropriate headers:

```bash
mkdir -p 01_GTDB_FASTA_REFORMAT
while read genome path; do \
  gunzip $path; \
  newpath="${path%.*}"; # remove .gz extension \
  outpath="01_GTDB_FASTA_REFORMAT/${genome}.fasta"; \
  reportpath="01_GTDB_FASTA_REFORMAT/${genome}_reformat_report.txt"; \
  anvi-script-reformat-fasta $newpath -o $outpath --simplify-names --prefix $genome --seq-type NT -r $reportpath; \
done < <(tail -n+2 GTDB_genomes_for_read_recruitment.txt)
```

Finally, we concatenated all of the reformatted sequences into one big FASTA file:

```bash
cat 01_GTDB_FASTA_REFORMAT/*.fasta >> GTDB_GENOMES.fasta

# make sure the number of contigs match
grep ">" 01_GTDB_FASTA_REFORMAT/*.fasta | wc -l
grep -c ">" GTDB_GENOMES.fasta
# 64280 in both cases 

# remove the reformatted genomes to save space
rm -r 01_GTDB_FASTA_REFORMAT/
```

The resulting file, `GTDB_GENOMES.fasta` will be our reference for the mapping workflow.

#### Running the read recruitment workflow

The metagenomics workflow in 'References Mode' is slightly different from the previous workflows that we have discussed. It takes a [fasta.txt file](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#fastatxt) describing the reference(s) to map against, turns each reference into a contigs database, does read recruitment from each of the samples in the [samples.txt file](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#samplestxt) with [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml), and summarizes the read-mapping results into individual profile databases.

Unfortunately, once again we have to treat the [Vineis et al.](https://doi.org/10.1128/mBio.01713-16) differently, because the workflow only works with paired-end reads and not with merged reads. So, when generating the input files for the workflow (in the `MISC/` folder), we leave out these samples:

```bash
# fasta txt
echo -e "name\tpath" > GTDB_GENOMES_fasta.txt
echo -e "GTDB_GENOMES\tGTDB_GENOMES.fasta" >> GTDB_GENOMES_fasta.txt

# samples txt
echo -e "sample\tr1\tr2" > GTDB_GENOMES_samples.txt
while read samp; do \
  grep $samp ../../01_METAGENOME_DOWNLOAD/METAGENOME_SAMPLES.txt >> GTDB_GENOMES_samples.txt; \
done < <(grep -v Vineis_2016 ../../TABLES/00_SUBSET_SAMPLES_INFO.txt | cut -f 1 | tail -n+2)
```

You'll find the configuration file for this workflow at `MISC/GTDB_GENOMES_mapping_config.json`. If you take a look, you will notice that most of the optional rules (i.e., gene annotation) are turned off (`"run": false`), and that references mode is turned on (`references_mode": true`). Here is how you can run the mapping workflow (after adjusting the config for your system, of course):

```bash
cp ../../MISC/GTDB_GENOMES_mapping_config.json .
anvi-run-workflow -w metagenomics  -c GTDB_GENOMES_mapping_config.json
```

For us, it took a few days to run. We stopped the workflow just before it ran `anvi-merge` (because, as you will see, we run that later to incorporate all the samples, including the Vineis et al. ones). If you let it run that step, it is fine - you simply may have to overwrite the resulting merged profile database later when you run `anvi-merge`.

#### Mapping the Vineis et al. samples separately

To map the merged reads from the 64 [Vineis et al.](https://doi.org/10.1128/mBio.01713-16) samples, we created a script that replicates the steps of the snakemake workflow for one sample with single-read input. You can find it at `SCRIPTS/map_single_reads_to_GTDB.sh`. The script should be run after the previous workflow finishes because it makes use of the `bowtie` index and output directories generated by the previous workflow (so that all the results are stored in one location). If you changed the output directories before you ran that workflow, you should update the paths in this script before you run it.

Here is the code to run the script (note that it uses 4 threads):

```bash
# get a file of sample paths
while read samp; do \
  grep $samp ../../01_METAGENOME_DOWNLOAD/METAGENOME_SAMPLES.txt | cut -f 2 >> vineis_samples.txt; \
done < <(grep Vineis_2016 ../../TABLES/00_SUBSET_SAMPLES_INFO.txt | cut -f 1 | tail -n+2)

while read samp; do \
  name=$(basename $samp | sed 's/.fastq.gz//'); \
  ../../SCRIPTS/map_single_reads_to_GTDB.sh $samp; \
done < vineis_samples.txt
```

Once those mapping jobs are done, you can summarize the read mapping results into a profile database for each sample. You can multithread the `anvi-profile` step by adding the `-T` parameter, if you wish. Note that these paths also make use of the output directories from the previous workflow:

```bash
while read samp; do \
  name=$(basename $samp | sed 's/.fastq.gz//'); \
  anvi-profile -c 03_CONTIGS/GTDB_GENOMES-contigs.db -i 04_MAPPING/GTDB_GENOMES/${name}.bam -o 05_ANVIO_PROFILE/GTDB_GENOMES/${name} -S $name; \
done < vineis_samples.txt
```

#### Putting it all together with `anvi-merge`

Once all samples have been mapped and profiled, you can merge all of the mapping results into one big profile database containing all samples:

```bash
anvi-merge -c 03_CONTIGS/GTDB_GENOMES-contigs.db -o 06_MERGED 05_ANVIO_PROFILE/GTDB_GENOMES/*/PROFILE.db
```

The resulting database will hold all of the coverage and detection statistics for these genomes across the gut metagenome dataset.

#### Summarizing the read mapping data

Once the coverage and detection data has been nicely calculated and stored in the merged profile database, we extracted the data into tabular text files for downstream processing. First, we created a collection matching each contig in the big FASTA file (which contains all 836 genomes) to its original GTDB genome. We imported that collection into the database, and then used the program `anvi-summarize` to summarize the coverage information across all contigs in a given genome.  Here is the code to do that:

```bash
grep "^>" GTDB_GENOMES.fasta | sed 's/^>//g' > contigs.txt                       # extract all contig headers
while read contig; do echo $contig | cut -d '_' -f 1-2 >> bins.txt; done < contigs.txt  # extract genome name from contig headers
paste contigs.txt bins.txt > GTDB_GENOMES_collection.txt
rm contigs.txt bins.txt

anvi-import-collection -c 03_CONTIGS/GTDB_GENOMES-contigs.db -p 06_MERGED/PROFILE.db -C GTDB_GENOMES --contigs-mode GTDB_GENOMES_collection.txt
anvi-summarize -c 03_CONTIGS/GTDB_GENOMES-contigs.db -p 06_MERGED/PROFILE.db -C GTDB_GENOMES"
```

The program produces a folder of various data tables, one of which is a matrix of detection of each genome in each gut metagenome. You will find this table in the DATAPACK at `TABLES/GTDB_GENOMES_detection.txt`.

### Subsetting gut genomes by detection in our sample groups

One issue with filtering for gut microbes using read recruitment to just one gene is that some non-gut microbes can slip into the set. For some groups of microbes, the Ribosomal Protein S6 gene is similar enough across different populations that they all end up in the same gene cluster, and if the cluster representative has high-enough detection in the HMP metagenomes, all of those genomes will be included in our list even if only a few of them are actually gut microbes. We noticed this problem when we started to look at the read recruitment data produced in the previous subsection - a lot of the genomes were undetected in the healthy and IBD gut metagenomes. When we investigated further, we saw that many of the undetected genomes were coming from two very large `Ribosomal_S6` gene clusters: one cluster of 180 `Enterobacteriaceae` genomes and one cluster of 163 `Streptococcus` genomes.

To mitigate this issue, we decided to take one more filtering step and remove any genomes that were irrelevant to our metagenome dataset based upon their low detection across those samples. We required the genomes to have at least 50% detection (of the entire genome sequence) in at least 2% of the healthy and IBD samples (which translates to at least 7 out of the 331 samples in those two groups). 

You can find the script we used to subset the genomes at `SCRIPTS/subset_gut_genomes_by_detection.py`. Here is how you run it:

```bash
# first, go back to the parent directory 04_GTDB_PROCESSING/
cd ..

python ../SCRIPTS/subset_gut_genomes_by_detection.py
```

It will generate 2 files: a shortened list of GTDB genomes that pass the filter (`genomes_detected_0.02_of_samples.txt`), and a subset of the detection matrix for just this set of genomes (`genomes_detected_0.02_of_samples-detection.txt`). 

### Final set of GTDB genomes and their contigs dbs

Ultimately, we ended up with a relatively small group of 338 gut microbial genomes from the GTDB that we used for downstream analyses. If you elected to download these genomes from [https://doi.org/10.5281/zenodo.7883421](https://doi.org/10.5281/zenodo.7883421) via the instructions at the top of this page, you will find their contigs databases in the `SUBSET_GTDB_CONTIGS_DBS` folder of the datapack. The file at `TABLES/01_GTDB_GENOMES_INFO.txt` describes each genome's taxonomy and also includes much of the data that we will generate in the next section.

Go back to the top-level datapack directory before you start the next section.

```bash
cd ..
```


## Metabolism and distribution analyses (for genomes)

This section covers the genome-level analyses that we ran on the set of 338 gut microbes that was just established. We analyzed the metabolic potential of each genome by calculating the stepwise completeness of each KEGG module. We then used the completeness scores of our 33 IBD-enriched pathways to determine whether each genome represented a microbe with high metabolic independence (HMI) - that is, high average completenesss of all these pathways - or not. Finally, we used our read recruitment results from the gut metagenome dataset to analyze the distribution of each group of genomes across healthy individuals and individuals with IBD. These investigations culminated in Figure 3 of the manuscript.

Just like in the metagenome metabolism section, we will work from the top-level directory of the datapack, but we'll make a folder in which we can generate some output:

```bash
mkdir 05_GTDB_ANALYSES
```

### Metabolism estimation in gut genomes

The file `GTDB_EXTERNAL_GENOMES.txt` holds all the paths to the genome databases in the `SUBSET_GTDB_CONTIGS_DBS` folder (if you downloaded it at the start of the workflow). To estimate metabolism on all of these genomes and produce matrix-formatted output, run the following:

```bash
anvi-estimate-metabolism -e GTDB_EXTERNAL_GENOMES.txt \
                  --matrix-format \
                  -O 05_GTDB_ANALYSES/GTDB_METABOLISM
```

The output file that you want to look at is the matrix of stepwise completeness scores, at `05_GTDB_ANALYSES/GTDB_METABOLISM-module_stepwise_completeness-MATRIX.txt`. It contains the completeness score for each KEGG module in each genome (completeness for an individual genome is analogous to PPCN for a metagenome). Later, we will use these data to make the completeness heatmap in Figure 3. Note that you can find this matrix in Supplementary Table 3b.

### The 'HMI score': labeling genomes by level of metabolic independence

In order to ask questions like "how many metabolically-independent genomes are there in IBD metagenomes?", we needed a way to determine whether a genome encodes high metabolic independence (HMI) or not. There is more than one way to be metabolically-independent, of course. In our [previous study of colonization after FMT](https://doi.org/10.1186/s13059-023-02924-x), we described one way. But in the current study, we care about the pathways that matter to microbial survival in the IBD gut environment. Therefore, our HMI genomes should have high completeness scores for the 33 IBD-enriched modules that we determined in a previous section.

How high should those scores be? That's a good question, and it doesn't have an obvious 'correct' answer. We looked at the distribution of completeness scores from these 33 modules, and didn't see an obvious threshold separating genomes into clear 'high' and 'low' groups. But we had to draw the line somewhere, and we settled on the requirement that the 33 IBD-enriched pathways should be, **on average, 80% complete** in order for a genome to be labeled as 'HMI'. We calculated an 'HMI score' for each genome by adding up the completeness scores of each of these modules. Since there are 33 IBD-enriched pathways each with a theoretical maximum completeness of 1.0, the maximum possible HMI score is 33. Applying our criteria leads to a threshold HMI score of 0.8 * 33 = 26.4, so **any genome with a score above 26.4 was considered to be an 'HMI' genome**, and all the others were labeled as 'non-HMI'.

We used the program `anvi-script-estimate-metabolic-independence` (which you can read about [here](https://anvio.org/help/main/programs/anvi-script-estimate-metabolic-independence/)) to score each genome using the stepwise completeness scores of the 33 IBD-enriched pathways. In the code below, we run that program on each genome and parse the resulting output files to produce a table of scores and labels for each genome at `05_GTDB_ANALYSES/HMI_scores.txt`.

```bash
# get just the module accession numbers as input to the program
tail -n+2 03_METABOLISM_OUTPUT/IBD_ENRICHED_MODULES.txt | cut -f 1 > 05_GTDB_ANALYSES/IBD_MODULE_ACC.txt

# compute the HMI score for each genome
mkdir 05_GTDB_ANALYSES/ESTIMATE_HMI
while read g path; do \
  anvi-script-estimate-metabolic-independence -c $path \
      --use-stepwise-completeness \
      --threshold 26.4 \
      --module-list 05_GTDB_ANALYSES/IBD_MODULE_ACC.txt \
      --kegg-data-dir KEGG_2020-12-23 > 05_GTDB_ANALYSES/ESTIMATE_HMI/${g}.txt 2>&1; \
done < <(tail -n+2 GTDB_EXTERNAL_GENOMES.txt)

# clean up the modules files that were generated
rm G*_modules.txt

# parse the output files into a table
echo -e "genome\tHMI_score\tclassification" > 05_GTDB_ANALYSES/HMI_scores.txt
while read g path; do \
  score=$(grep "Genome score" 05_GTDB_ANALYSES/ESTIMATE_HMI/${g}.txt | cut -d ' ' -f 4) ;\
  label=$(grep "Metabolic independence" 05_GTDB_ANALYSES/ESTIMATE_HMI/${g}.txt | cut -d ' ' -f 4) ;\
  if [[ "$label" = "High" ]]; then \
    class="HMI"; \
  else class="non-HMI"; \
  fi; \
  echo -e "${g}\t${score}\t${class}" >> 05_GTDB_ANALYSES/HMI_scores.txt; \
done < <(tail -n+2 GTDB_EXTERNAL_GENOMES.txt)
```

You can also find the scores and classifications in the `TABLES/01_GTDB_GENOMES_INFO.txt` file.

### Detection calculations

We wanted to quantify the distribution of each genome (and each class of genomes) across the healthy and IBD sample groups. To do this, we summarized the detection data from the mapping workflow by calculating the following proportions:

* a per-genome ratio of detection in the _healthy_ sample group, calculated by dividing the number of healthy samples in which the genome was reliably (> 0.5x) detected by the total number of healthy samples (n = 229)
* a per-genome ratio of detection in the _IBD_ sample group, calculated by dividing the number of IBD samples in which the genome was reliably (> 0.5x) detected by the total number of IBD samples (n = 101)
* a _per-sample_ ratio of HMI genome detection, calculated by dividing the number of HMI genomes reliably (> 0.5x) detected in the sample by the total number of genomes (of any class) reliably (> 0.5x) detected in that sample

The script at `SCRIPTS/summarize_GTDB_detection.py` calculates these proportions, using the detection matrix for the subset of 338 GTDB genomes that was generated earlier (at `04_GTDB_PROCESSING/genomes_detected_0.02_of_samples-detection.txt`). If you don't have that input file, you should go back to the previous section and generate it. Here is how to run the script for computing the proportions:

```
python SCRIPTS/summarize_GTDB_detection.py
```

It will generate two output files in the `05_GTDB_ANALYSES/` folder: `genome_proportions.txt`, which contains the two per-genome proportions, and `sample_proportions.txt` which contains the per-sample proportions. Note that you can also find these data in Supplementary Table 3e and 3f.

### Percent abundance calculations

To calculate the percent abundance of each genome in each sample, we used the following equation: 
`% abundance = sum(num reads mapping to contigs in genome) / total reads in sample`

The total number of reads in each sample was already calculated above (and can be found in the `r1_num_reads` and `r2_num_reads` columns of `TABLES/00_SUBSET_SAMPLES_INFO.txt`). To calculate the number of reads mapping to each contig in each genome, we ran `samtools idxstats` on the indexed BAM files from the mapping workflow that we ran in the `04_GTDB_PROCESSING/GTDB_MAPPING_WORKFLOW/` folder. To replicate that (if you ran the mapping workflow), you can run the following code:

```bash
# go to the place where the workflow output lives
cd 04_GTDB_PROCESSING/GTDB_MAPPING_WORKFLOW
# make a folder to store the idxstats output
mkdir 04_MAPPING/IDXSTATS

# this is for the 344 non-pouchitis samples
while read samp r1 r2; do \
  samtools idxstats 04_MAPPING/GTDB_GENOMES/${samp}.bam > 04_MAPPING/IDXSTATS/${samp}.idxstats; \
done < <(tail -n+2 GTDB_GENOMES_samples.txt)
# this is for the 64 pouchitis samples
while read samp; do \
  name=$(basename $samp | sed 's/.fastq.gz//'); \
  samtools idxstats 04_MAPPING/GTDB_GENOMES/${samp}.bam > 04_MAPPING/IDXSTATS/${samp}.idxstats; \
done < vineis_samples.txt
```

This will generate a table for each sample in the `04_GTDB_PROCESSING/GTDB_MAPPING_WORKFLOW/04_MAPPING/IDXSTATS/` folder that describes each contig from the `GTDB_GENOMES.fasta` file, the contig length, the number of reads mapping to it, and the number of unmapped reads (which is always 0 because the mapping workflow parameters ignore any unmapped reads).

Once that is done, you can run the following script to generate the percent abundance values, including the percent abundance ratio that we plot in Figure 3. The output files will be `percent_abundances.txt` and `percent_abundance_averages.txt`.

```
python ../../SCRIPTS/get_percent_abundance.py
```

After this, you should go back to the top-level directory.

```bash
cd ../..
```

### Generating the genome phylogeny

To generate the phylogenetic tree that organizes the genomes in Figure 3, we used a set of 10 ribosomal proteins that were annotated in at least 90% (n = 304) of the 338 GTDB genomes. We aligned each set of gene sequences, concatenated the alignments, removed any positions with too many gaps using [`trimal`](http://trimal.cgenomics.org/trimal), and built the tree with [IQ-Tree](http://www.iqtree.org/). Here is the code to do that:

```bash
# get aligned and concatenated ribo protein seqs
anvi-get-sequences-for-hmm-hits -e GTDB_EXTERNAL_GENOMES.txt \
                  --return-best-hit \
                  --concatenate-genes \
                  --get-aa-sequences \
                  --hmm-sources Bacteria_71 \
                  -o 05_GTDB_ANALYSES/genomes_hmm_hits_aligned_concatenated.txt \
                  --gene-names Ribosomal_S6,Ribosomal_S16,Ribosomal_L19,Ribosomal_L27,Ribosomal_S15,Ribosomal_S20p,Ribosomal_L13,Ribosomal_L21p,Ribosomal_L20,Ribosomal_L9_C

# remove positions with too many gaps
trimal -in 05_GTDB_ANALYSES/genomes_hmm_hits_aligned_concatenated.txt \
       -out 05_GTDB_ANALYSES/genomes_hmm_hits_aligned_concatenated_GAPS_REMOVED.fa \
       -gt 0.50

# make tree 
iqtree -s  05_GTDB_ANALYSES/genomes_hmm_hits_aligned_concatenated_GAPS_REMOVED.fa \
  -nt 10 -m WAG -B 1000 --prefix  05_GTDB_ANALYSES/GTDB_tree
```

This will produce a tree file at `05_GTDB_ANALYSES/GTDB_tree.contree`, which you can use in the next subsection to visualize all the data that we have generated in this section so far. You can also find our tree file at `MISC/GTDB_tree.contree`.

### Generating Figure 3


## Machine learning analyses

### Classifying gut metagenomes as IBD vs Healthy

### Generating Figure 4

### Classifying antibiotic time-series metagenomes from Palleja et al

### Generating Figure 5


## Supplementary Analyses

### Exploring annotation efficiency (SF 2)

### Additional comparisons of metabolic pathways (SF 3)

### Examining cohort effect (SF 4)

### Testing classifier generalizability
---
layout: page
title: A reproducible workflow for Veseli et al, XXXX
modified: 2026-09-17
excerpt: "A bioinformatics workflow for developing a distribution of coverage metric."
comments: true
authors: [iva]
---

Here is a list of links for quick access to the data described in our manuscript and on this page:

* TBD
* TBD

{:.notice}
If you have any questions, notice an issue, and/or are unable to find an important piece of information here, please feel free to leave a comment down below, send an e-mail to [us](/people/), or get in touch with us through Discord:

{% include _join-anvio-discord.html %}

## Study description

TBD.

## Downloading the data for this reproducible workflow

TBD.

## Computational environment details

TBD.

## Obtaining publicly-available ocean metagenomes

Our group works a lot with publicly-available marine metagenomes in various projects, so it was a natural choice to reuse some of these datasets for developing and testing the discov metric. Since we would need to subsample some metagenomes to mimic lower sequencing depths and we would test the application of discov for studying biogeography, our criteria for selecting which datasets to use was as follows:
- they had to include at least some deeply sequenced samples to allow confident assessment of genome presence-absence and subsampling
- they had to collectively span multiple ocean regions (temperate & polar, open ocean & coastal), depths and size fractions
- they had to include accurate and complete metadata for sample location (latitude, longitude), depth, and size fraction
- they had to include paired-end samples for consistent read-mapping behavior across all samples

We settled on the following set of sequencing datasets from various ocean sampling efforts:

|**`Project`**|**`Name`**|**`Num Samples (paired-end)`**|**`References`**|
|:--|:--|:--|:--|
|PRJEB1787|Tara Oceans (prokaryotic size fraction)| 136 | [Sunagawa et al 2015](https://www.science.org/doi/10.1126/science.1261359)
|PRJEB9740|Tara Arctic (prokaryotic size fraction)| 41 | [Sunagawa et al 2015](https://www.science.org/doi/10.1126/science.1261359)
|PRJEB4352|Tara Oceans (protist size fraction)| 805 | [Carradec et al 2018](https://www.nature.com/articles/s41467-017-02342-1)
|PRJEB1788|Tara Oceans (prokaryotes and large DNA viruses size fraction)| 63 | [Sunagawa et al 2015](https://www.science.org/doi/10.1126/science.1261359)??
|PRJEB4419|Tara Oceans (viral size fraction)| 90 | [Brum et al 2015](https://doi.org/10.1126/science.1261498), [Roux et al 2016](https://doi.org/10.1038/nature19366), [Gregory et al 2019](https://doi.org/10.1016/j.cell.2019.03.040)
|PRJEB9742|Tara Arctic (viral size fraction)| 41 | [Brum et al 2015](https://doi.org/10.1126/science.1261498), [Roux et al 2016](https://doi.org/10.1038/nature19366), [Gregory et al 2019](https://doi.org/10.1016/j.cell.2019.03.040)
|PRJEB8682|Ocean Sampling Day (OSD) 2014| 150 | [Kopf et al 2015](https://doi.org/10.1186/s13742-015-0066-5), https://doi.org/10.1594/PANGAEA.854419
|PRJEB40760|OSD 2018| 52 | https://marineinfo.org/en/doc/dataset/7916
|PRJEB40764|OSD 2019| 45 | https://marineinfo.org/en/doc/dataset/7917
|Malaspina_Acinas_et_al (various BioProjects)|Malaspina| 58 | [Acinas et al 2021](https://www.nature.com/articles/s42003-021-02112-2), [Duarte et al 2015](https://doi.org/10.1002/lob.10008)
|PRJEB52452|Malaspina MProfile| 76 | [Sanchéz and Coutinho et al 2024](https://www.nature.com/articles/s41597-024-02974-1)
|PRJEB83083|Antarctic Circumnavigation Expedition (ACE)| 218 | [Faure et al 2026](https://doi.org/10.1038/s41467-026-69584-w)
| TOTAL | | **1,775** | |

### Download & QC 
We used a standardized procedure to download and process each dataset (individually) on our high-performance computing cluster. Here it is:

1. Obtain a list of NCBI SRA run accessions (WGS samples only) associated with the BioProject, as well as a table mapping BioSamples to their component runs (you can usually get this information directly from SRA metadata files for the BioProject)
2. Use the sra-download workflow in anvi'o to download raw FASTQ files for each run accession
3. Use the metagenomics workflow in anvi'o to combine all runs belonging to the same BioSample and perform quality control with the Illumina utils program `iu-filter-quality-minoche`

<details markdown="1"><summary>Show/Hide  Workflow commands and configuration files </summary>

Here are the minimal workflow commands (without the specific flags needed for running them on our HPC):
```bash
# step 2
anvi-run-workflow -w sra_download -c download_config.json -A --keep-going
# step 3
anvi-run-workflow -w metagenomics -c QC_config.json -A --until gzip_fastqs --keep-going
```
TBD Example config files for the workflows are available in the datapack at TBD.
</details>

With this, we obtained a folder of QC'ed, paired-end FASTQ files (one per BioSample) for each BioProject. We then put the absolute paths to all FASTQ files into a [samples-txt table](https://anvio.org/help/main/artifacts/samples-txt/) that we could use for easy access to the samples from anywhere on our cluster.

### Metadata

For most of the datasets, we were able to obtain the metadata we needed directly from the NCBI SRA Run Selector, where we downloaded the per-sequencing-run metadata for all runs in a given BioProject. For OSD 2014, we combined its SRA metadata with metadata from Pangaea (https://doi.pangaea.de/10.1594/PANGAEA.854419), and for the two Malaspina cruises, we got the metadata from supplementary tables of the corresponding papers (Supp. Table 1 from [Acinas et al 2021](https://www.nature.com/articles/s42003-021-02112-2) and Supp. Table 2 from [Sanchéz and Coutinho et al 2024](https://www.nature.com/articles/s41597-024-02974-1)).

We then standardized and combined a subset of the metadata columns across all projects into one tab-delimited table (Supp. Table TBD of our paper), which is available in the datapack at TBD. Using `pandas`, we specifically extracted size fraction (upper and lower thresholds), depth, and location information (using 'Latitude_Start' and 'Longitude_Start' values in SRA tables when both starting and ending coordinates were available). We converted missing entries of various forms ('NA', 'not provided', etc) to blanks, standardized the column names, and dropped duplicate rows (which occurred for BioSamples with more than one sequencing run). As each initial metadata table was formatted slightly differently, the metadata processing code was correspondingly different across projects, but here is a representative example for BioProject PRJEB40760 (OSD 2018):

```python
import pandas as pd
df = pd.read_csv("OSD_2018_SRA_metadata.txt", sep="\t")
cols_of_interest=['BioSample', 'BioProject', 'Latitude_Start', 'longitude_start', 'Depth', 'sample_size-fraction_lower-threshold', 'sample_size-fraction_upper-threshold']
sub = df[cols_of_interest]
# convert to NA when applicable
sub.replace('not provided', '', inplace=True)
sub.replace('no prefiltration', '', inplace=True)
sub.replace('no size limit', '', inplace=True)
sub.replace('No prefiltration', '', inplace=True)
rename_cols = {'Latitude_Start':'latitude', 'longitude_start':'longitude', 'Depth':'depth', 'sample_size-fraction_lower-threshold':'size_fraction_lower_threshold', 'sample_size-fraction_upper-threshold':'size_fraction_upper_threshold'}
sub.rename(columns=rename_cols, inplace=True)
sub.to_csv("PRJEB40760_SRA_metadata.txt", sep="\t", index=False)
```

Once each project's metadata file had the same columns and format, we combined them into one with a simple BASH loop:
```bash
head -n 1 PRJEB9742_SRA_metadata.txt > sample_metadata.txt; 
for f in *metadata*.txt; do tail -n+2 $f >> sample_metadata.txt; done
```

We later added a column indicating the human-readable name of the sequencing effort:

```python
import pandas as pd

METADATA_FILE="00_DATA/sample_metadata.txt"
df = pd.read_csv(METADATA_FILE, sep="\t", index_col=0)

# any BioProject not in this list is associated with the original Malaspina cruise
bioproject_to_effort_mapping = {
"PRJEB1787": "Tara Oceans",
"PRJEB9740": "Tara Arctic",
"PRJEB4352": "Tara Oceans",
"PRJEB1788": "Tara Oceans",
"PRJEB4419": "Tara Oceans",
"PRJEB9742": "Tara Arctic",
"PRJEB8682": "OSD 2014",
"PRJEB40760": "OSD 2018",
"PRJEB40764": "OSD 2019",
"PRJEB52452": "Malaspina MProfile",
"PRJEB83083": "ACE",
}

df['sampling_effort'] = df['BioProject'].map(bioproject_to_effort_mapping).fillna('Malaspina')
df.to_csv(METADATA_FILE, sep="\t")
```

And we dropped any samples that were missing latitude/longitude metadata. This affected ~23 samples from the Tara Oceans project PRJEB4352 (protist size fraction), which in fact were missing _all_ core metadata fields and not just the coordinates. There were originally 828 samples from this project, but removing the 23 with missing metadata brought the number down to 805.

```python
import pandas as pd

METADATA_FILE="00_DATA/sample_metadata.txt"
df = pd.read_csv(METADATA_FILE, sep="\t", index_col=0)

na_samples = df[df.latitude.isna()].index.to_list()
df.drop(labels=na_samples, inplace=True)
df.to_csv(METADATA_FILE, sep="\t")
```

### Sequencing Depth

We used the number of paired sequencing reads to quantify the sequencing depth of each BioSample. Counting the number of reads in a FASTQ file can take a while, but some bioinformatics tools report the number of reads as part of their output -- we luckily had `bowtie2` logs available from previous read recruitment analyses with many of these samples. When possible, we extracted the number of reads from the `bowtie2` logs, and for samples we hadn't mapped yet, we ran `seqfu count` on the R1 FASTQ files to get a table of counts (R1 and R2 counts were the same because the earlier QC step removed any unpaired reads).

The datapack includes a script to extract the count information at TBD. You can modify the variables at the top of the script to give it access to 1) a folder of `bowtie2` logs and 2) the output of `seqfu count`, then run it like this:
```bash
bash 01_SCRIPTS/get_read_counts.sh
```

We later added these read count data into the metadata file:
```python
import pandas as pd

METADATA_FILE="00_DATA/sample_metadata.txt"
READS_FILE="00_DATA/sample_read_counts.txt"
df = pd.read_csv(METADATA_FILE, sep="\t", index_col=0)
reads = pd.read_csv(READS_FILE, sep="\t", index_col=0)

df = df.join(reads)
df.to_csv(METADATA_FILE, sep="\t")
```

## Obtaining OMDB v2 species cluster reference genomes

TBD.

## Generating a test dataset of manually-verified present/absent genomes

For discov to work, it needs to quantitatively distinguish between read recruitment patterns on _truly present_ genomes and those on _truly absent_ genomes. We therefore needed a set of known 'present' genomes and known 'absent' genomes to compare as we developed the metric. To identify these genomes, we reasoned that we could look at read recruitment data from extremely deeply sequenced metagenomes, use detection and coverage thresholds to find a subset of genomes that are likely to be present or absent, and then visualize the read recruitment data for the subset to manually verify which are truly present and which are truly absent. We required deeply-sequenced metagenomes for this to confidently assess presence-absence; for instance, to ensure genomes of low-abundance populations do not appear absent due to lack of sufficient sequencing coverage.

We selected the top 25 most deeply-sequenced metagenome samples to be our test set:

```python
import pandas as pd

METADATA_FILE="00_DATA/sample_metadata.txt"
df = pd.read_csv(METADATA_FILE, sep="\t", index_col=0)

df.sort_values("num_read_post_QC", ascending=False).iloc[0:25,].to_csv("00_DATA/test_samples_metadata.txt", sep="\t")
```

We then mapped these 25 samples to all 32,022 reference genomes from OMDB v2 using the [anvi'o metagenomics workflow](https://anvio.org/help/main/workflows/metagenomics/#references-mode) in 'references mode'. First, we combined all the genome sequences into a single reference FASTA file called `ALL_OMD_REFERENCE_GENOMES.fa` and ran the workflow until the `bowtie2-build` step to generate the index files using a high-memory node on our HPC. Afterwards, we restarted the workflow (on multiple regular nodes) and ran it until the `anvi-init-bam` step to produce sorted, indexed BAM files. Here are the minimal workflow commands:

```bash
# initial workflow run just to generate bowtie2 index (requires high memory)
anvi-run-workflow -w metagenomics -c mapping_config.json -A --until bowtie_build
# follow-up workflow run for the read recruitment and processing steps
anvi-run-workflow -w metagenomics -c mapping_config.json -A --until anvi_init_bam --keep-going
```
You can find the config file for the mapping in the datapack at TBD.

When the read recruitment workflow finished, we generated a [contigs database](https://anvio.org/help/main/artifacts/contigs-db/) from the reference FASTA file (this also required a high-memory node), and then we created a [collection-txt file](https://anvio.org/help/main/artifacts/collection-txt/) that mapped each contig name to its genome of origin. We used both of these files to obtain per-genome detection and coverage statistics by profiling the BAM files with [anvi-profile-blitz](https://anvio.org/help/main/programs/anvi-profile-blitz/):

```bash
anvi-gen-contigs-database -f ALL_OMD_REFERENCE_GENOMES.fa \
                          -n ALL_OMD_REFERENCE_GENOMES \
						  -T 100 \
						  -o ALL_OMD_REFERENCE_GENOMES-CONTIGS.db

# make a collection txt file
grep "^>" ALL_OMD_REFERENCE_GENOMES.fa | sed 's/>//' > contigs
cat contigs | rev | cut -d'_' -f2- | rev > bins
paste contigs bins > ALL_OMD_REFERENCE_collection.txt
rm contigs bins

anvi-profile-blitz 04_MAPPING/ALL_OMD_REFERENCE_GENOMES/*.bam \
                   -c ALL_OMD_REFERENCE_GENOMES-CONTIGS.db \
                   -C ALL_OMD_REFERENCE_collection.txt \
                   -o omdb_genomes-profile_blitz.txt
```

We used the data in `omdb_genomes-profile_blitz.txt` to select a subset of genomes for manual inspection. To do so, we ran the following script (available in the datapack):

```bash
bash 01_SCRIPTS/extract_candidate_genomes.sh
```

Briefly, this script does the following for each of our 25 deeply-sequenced metagenomes:
- extracts up to 10 candidate 'present' genomes from `omdb_genomes-profile_blitz.txt` that have **> 0.85 detection** in the sample
- extracts up to 10 candidate 'absent' genomes from `omdb_genomes-profile_blitz.txt` that have **<= 0.15 detection and > 2 mean coverage** in the sample. The coverage requirement ensures that the mapped reads are not entirely uniformly distributed across the genome, which could be the case for genomes that are present but in extremely low abundance
- extracts the contig sequences of each candidate genome from `ALL_OMD_REFERENCE_GENOMES.fa`
- generates a [contigs database](https://anvio.org/help/main/artifacts/contigs-db/) containing all candidate genome sequences
- subsets the sample's BAM file to create a smaller one containing only the candidate genome contigs
- runs `anvi-profile` on the subset BAM file to get a [profile database](https://anvio.org/help/main/artifacts/single-profile-db/) for visualizing their read recruitment data
- creates a [collection-txt file](https://anvio.org/help/main/artifacts/collection-txt/) mapping contig name to genome of origin and imports this collection into the profile database
- annotates ribosomal RNA genes and single-copy core genes in these genomes with [`anvi-run-hmms`](https://anvio.org/help/main/programs/anvi-run-hmms/)

A few samples did not yield enough candidate genomes using these initial detection/coverage parameters: SAMEA2620786 only had 4 'absent' candidates (needed 1 more), SAMEA2656986 only had 3 'present' candidates (needed 2 more), and SAMEA2657075 did not have any 'present' candidates (needed 5). For these samples, we relaxed the parameters of the extraction script as follows and ran it again to capture more candidates:
- extract up to 10 candidate 'present' genomes that have **> 0.75 detection** in the sample
- extract up to 10 candidate 'absent' genomes that have **<= 0.20 detection and > 1 mean coverage**
For these samples, we first validated as many genomes as we could from the first-pass before adding in candidates from the second pass.

Using the databases created by the extraction script, we were able to visualize the read recruitment data to all candidate genomes using `anvi-interactive`. We manually examined the coverage patterns across contigs in a given genome using the inspection page (equivalent to running `anvi-inspect` from the command line on different splits of contigs in the genome), paying attention to how the mapped reads were distributed. Specifically, we assessed: (1) how _spread_ out the reads were across the genome and (2) how _even_ the coverage depth was, roughly quantifying these attributes on a scale from 1 to 4. Here is what the scale meant for _spread_:
1. Poor spread. Reads map to relatively few locations on the genome. Often they only map to short contigs, which may represent highly-shared sequences that assemble poorly or contaminating sequences within a metagenome-assembled genome. Within a given contig, reads map so sparsely that there are multiple-kilobases long stretches of sequence with no coverage. Often mapped reads are concentrated within one gene, or a gene neighborhood.
2. Okay spread. Reads map to several locations in the genome, often including some longer contigs. However, the distribution of reads within contigs is highly irregular -- sometimes entire contigs or sections of contigs are fully covered (potentially a shared prophage or genomic island) while everything around it is entirely uncovered, and sometimes reads map only in a few isolated peaks. If reads map regularly within a given contig, there is plenty of distance between them (at least a few kilobases), and this pattern is not consistent across many contigs.
3. Good spread. Reads map to many locations across the genome and most of the contigs (especially longer ones) have at least some coverage. Within a contig, reads map regularly and not very sparsely (e.g., with some reads mapping every 2-4kb), and this pattern is consistent across many contigs. When contigs are entirely uncovered, they are often short. Subsections of longer contigs can be missing coverage, but this often corresponds to clear gene neighborhoods of functionally-related genes or mobile genetic elements.
4. Great spread. Reads map to almost everywhere in the genome, including all longer contigs. There are few gaps, and these gaps correspond to clear gene neighborhoods or mobile genetic elements.

And here is the scale for _evenness_:
1. Very uneven. There are sharp peaks of coverage, and the magnitude of the coverage depth is usually vastly different at distinct locations in the genome. Often there are drastic examples of non-specific read recruitment (i.e., seeing a peak at ~5,000x coverage over one gene when coverage elsewhere is sparse and shallow).
2. Relatively uneven. There are some sharp peaks and some wider stretches of coverage at the same level of depth. The magnitude of coverage depth can fluctuate between a few orders of magnitude but is not always drastically different across distinct regions of the genome.
3. Relatively even. There is a clear median depth of coverage but across different bases the coverage can fluctuate quite a lot; it appears rather 'spiky'. Sometimes there are longer regions of the genome where the coverage depth is quite distinct from its usual pattern. Exceptions can sometimes be explained by finding functionally-related genes (e.g., transport systems, pilus genes, defense mechanisms) or genes related to mobile genetic elements (e.g., transposases, prophage genes).
4. Very even. Almost all regions of coverage have roughly the same coverage depth (e.g., fluctuating around a clear median value with depths within the same order of magnitude, like 5-10x or 30-60x or 100-150x). Occasional spikes or drops in coverage over certain regions (individual genes or gene neighborhoods) are clearly outliers and coverage depth eventually returns back to the same level it was before.

We combined all the evidence to manually confirm whether each candidate genome was 'present' or 'absent' in the sample. We stopped once we had manually verified 5 'present' genomes and 5 'absent' genomes per sample.

{:.warning}
The candidate genomes were not randomly selected from the [anvi-profile-blitz](https://anvio.org/help/main/programs/anvi-profile-blitz/) output -- we simply took the first 10 that fit our threshold criteria. That led to some overlap in the candidates across different samples, especially for the absent genomes (since it is much more likely for a given genome to be absent than present in a sample). 'Absent' candidates often exhibited the same read recruitment patterns across different samples, which makes sense given that sequences prone to recruiting non-specific reads will tend to recruit similar reads from all samples. This led to some redundant examples in our test set.

For 'present' candidates, the read recruitment distribution often followed similar patterns and clearly indicated presence of the genome. Our initial high detection threshold for selecting candidates ensured that the reads were usually spread across most contigs in the genome (most had spread assessed at 4, some at 3), though coverage depth could sometimes vary substantially (evenness of 2-4). For these candidates, we started from one end of the genome and inspected each sequential split until an overall pattern was identified, and then we spot-checked other random regions of the genome to ensure the read recruitment patterns remained consistent. We also spot-checked splits with distinct- or unusual-looking mean coverage values in the main interface page. We paid attention to genes with missing coverage and genes with unusual coverage depth (compared to the rest of the genome) to ensure that these could generally be attributed to potential biological variation between the reference genome and an environmental population, such as differential functional content, presence of mobile genetic elements, or non-specific read recruitment. All 'present' candidates we validated were truly present in the metagenome. In one case (`TARA_SAMEA4398432_MAG_00000068` in sample `SAMEA4398432`), the candidate genome appeared to contain sequences from two distinct microbial populations, both of which were present in the metagenome sample at different abundance levels (there were two different 'modes' of coverage depth) -- however, since all sequences within this potentially-contaminated genome were clearly present in the sample, we kept it in our validated set.

For 'absent' candidates, there was much more variation in the read recruitment patterns and these needed to be scrutinized more carefully. Many contigs in these genomes had no coverage at all; we did not visually inspect most of these but we paid attention to the contig lengths to get a sense of how long the uncovered regions were. We then inspected many contigs with nonzero mean coverage to identify the overall mapping pattern, paying attention to how many peaks of coverage there were, how deep the coverage was at each peak, how long the gaps were between regions with coverage, and how the peaks were distributed within a contig and relative to the locations of genes. If there was a relatively consistent mapping pattern across many contigs of the genome, we didn't exhaustively inspect each contig but rather spot-checked contigs in other regions of the genome to ensure they had the same pattern. Several of these candidate genomes were not clearly absent -- if a genome had reads mapping somewhat regularly with small but mostly even coverage depth and this pattern was consistent across many contigs, we threw out that candidate. For genomes that were verifiably absent, we tried to prioritize including genomes with different 'archetypes' of read recruitment patterns in our test set -- that is, including some absent genomes with higher spread and some absent genomes with higher evenness so the test set was representative of these different patterns. Ultimately, most verified absent genomes had a spread of 1-2 and evenness of 1-2, though evenness was sometimes quite high (3-4) in genomes where very few reads mapped without many clear peaks of non-specific read recruitment.

## Identifying mathematical ways to quantify coverage 'spread' and 'evenness'

We used Claude.AI (Sonnet and Opus models) to help us come up with metrics for quantifying the degree to which mapped reads are spread across the full length of a reference sequence, and the degree to which coverage depth was uniform. We used prompts such as the following:

> I am developing a metric for evaluating the distribution of read recruitment across a nucleotide sequence, and I'd eventually like your help to create a mathematical formula for it. Before we get there, I want to extensively discuss the context, needs, and assumptions of the metric with you.
>
> First, to describe the goals of the metric: we want to be able to robustly detect the presence or absence of a microbial population in the environment using read recruitment of metagenomic reads to its genome, even at low sequencing depths. Currently we use a threshold on 'detection', the proportion of bases in a sequence that have at least one read mapped to it. However, detection can be very low for low-abundance populations that are present in the metagenome sample; for instance, when they have <1x coverage in the sample. We don't want to miss those by using an arbitrarily high detection threshold.
>
> Instead, we want to rely on the geographic distribution of the metagenomic coverage across the genome to determine whether it is truly present or not. The assumption here is: if a population is present even at very low (<1x) coverage, the few reads that do map will be spread out roughly evenly or randomly across the reference genome sequence.
>
> The counter-example to this is the case of non-specific read recruitment, in which the population is not truly present in the sample, but because its genome contains some regions that are widely shared across other microbes, those regions pick up reads originating from other microbial populations. These regions will typically appear as isolated spikes in the read recruitment data.
>
> Do you have any clarifying questions so far? Do these assumptions make sense to you?

and the following prompt describing an early version of the discov calculation to inspire alternative options:

> My planned calculation includes considering the zero-coverage regions as 'gap' regions. Given an array of per-base coverages across a contig, this is what I would do:
>1) first, I would walk over the contig to identify each covered region (bases with coverage >0) and each gap region (bases with coverage == 0). I would store the length of these regions as well as the average coverage within the region (which would be 0 for the gaps, allowing us to easily distinguish the gap regions later).
>2) next, I would do a pre-filtering step to remove regions that appear to have non-specific read recruitment based on the 'expected' depth of coverage. The assumption here is that non-specific read recruitment would have higher levels of coverage than the rest of the genome, so we can remove high depth outliers (for instance, regions of coverage where the depth is >2 std deviations away from the overall mean). In the extreme case where a genome has mostly non-specific read recruitment and little to no population-specific mapping, this is a bit more difficult because the non-specific regions will determine the overall mean coverage -- however, in these cases, the disparity between the overall mean coverage and the low detection value can help us set a global coverage threshold for filtering regions out. For instance, if we see regions of 10X coverage when the detection is 0.1, then we know those regions are likely non-specific even if they make up most of the length of the covered regions.
>3) then, I would use the lengths or locations of the remaining coverage/gap regions to quantify the uniformity of the coverage distribution. Ideas to do this include computing entropy as a measure of randomness (for instance, of the gap lengths), or doing a chi-squared test for a discrete uniform distribution, or simply computing variance.
>4) finally, I would ideally want to normalize the uniformity measure by the length of the contig somehow.
>
> The end result should be one number that quantifies the distribution of coverage across the length of the contig. There are two situations that would result in the metric having a high value: (1) the overall detection is high (almost everything has coverage), in which case there are few gaps to begin with as coverage is distributed throughout the contig. or (2) overall detection is low (even very low), but the coverage is evenly distributed throughout the contig.
>
> I would welcome alternatives for doing the calculation; for instance, if there is a way to work directly with the per-nucleotide coverage array rather than doing (1) to collapse into gap vs covered regions, or if there is a way to skip the prefiltering described in (2) while still being able to avoid the undue influence of non-specific read recruitment.
>
> What are your thoughts or questions on this part?

Claude.AI proposed a number of metrics, and we followed-up via discussions with the model on the calculation, interpretation, and limits of each metric. After lots of back-and-forth, discarding some ideas and coming up with new ones along the way, these discussions enabled us to select a subset for testing. Here are a few examples of things we asked the model:

> Could you please explain to me the differences between CV versus Gini coefficient (and how these are calculated)?
and
> What is the lower limit of the number of regions between which the spatial evenness metric (CV or Gini) breaks down? Does the number of regions have to be >2 or will the calculation be more robust at higher numbers?
and
> I'm curious about the Gini coefficient's limitations as a measure of inequality. Specifically, how does it behave when the number of observations is small? How many observations is 'too few' to reliably compute it?

Much later in the study, once we determined that we would need to incorporate the regular detection metric into version 1.0 of the discov score, we again turned to Claude with the following prompt:

> I'm working with read recruitment data from metagenomes mapped to microbial genomes, and am developing a 'distribution of coverage' (discov) score that quantifies (1) spread of coverage across the genomic reference and (2) evenness of coverage depth at all bases with nonzero coverage. The idea is to use this metric to identify when genomes are present in a sample, even if these genomes have low abundance. These are the current formulations of discov: 'linear' (DisCov = αS + (1-α)E) and 'geometric' (DisCov = S^α * E^(1-α)), where S is the spread metric, E is the evenness metric and α is a weighting factor between 0 and 1. Since S and E also range from 0-1, discov scores range from 0-1. I can tell you more about how S and E are computed, if necessary.
>
> Now, I want to add detection into this formula. Detection is the fraction of bases in the reference sequence that have at least 1 read mapping to them. When detection is high (often we use >0.25 as a threshold), we are fairly confident that the genome is present in the sample. But when detection is low, the genome could still be present yet at low abundance. So what I want is to add detection into the formula with an inverse weight, such that when detection is high, it can boost the value of the score but when detection is low, it doesn't impact the score much or at all. I guess I want it to be additive only. And another requirement -- it would be best for the score to retain a range of 0-1 as this is most interpretable.
>
> My question to you is, what are the mathematical options for incorporating detection into these discov formulas in an additive, inversely weighted fashion while maintaining the overall 0-1 range?

A full report of our exchanges with Claude is out of the scope of this reproducible workflow, but anyone looking for more details should feel free to contact Iva about this.
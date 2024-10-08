---
layout: post
authors: [ahenoch]
title: "Leveling up pangenomics with interactive pangenome graphs"
excerpt: "How to create, understand and appreciate pangenomic graph structures."
modified: 2024-07-10
tags: []
categories: [anvio]
comments: true
redirect_from:
---

{% include _project-anvio-version.html %}

{:.notice}
{% include _fixthispage.html source="_posts/anvio/2024-07-10-pangenome-graphs.md" %}

{% capture images %}{{site.url}}/images/anvio/2024-07-10-pangenome-graphs{% endcapture %}

With the workflow on anvio pangenome graphs you can,

* **Identify and analyse** regions of hypervariability in your genomes
* **Estimate similarity** between the genomes based on graph structures instead of e.g. genome similarity
* **Interactively visualize** your pangenome graph and the impact of single genomes on the consensus
* **Highlight** the distribution of paralogous gene clusters through their synteny based sub clusters
* **Summarize** distribution and frequency of common graph motives into a nice table for downstream analysis

{:.notice}
You can use the anvi'o workflow for pangenome graphs even if you haven't done any metagenomic work with anvi'o. All you need is an anvi'o installation, and a FASTA file for each of your genomes.

## Introduction

The anvi'o pangenomic workflow described here will walk you through the following steps:

* Using the anvi'o pangenomics workflow to generate an anvi'o {% include ARTIFACT name="genomes-storage-db" %} and  {% include ARTIFACT name="pan-db" %} using the program {% include PROGRAM name="anvi-run-workflow" %}

* Explain how to define a subset of the sequences to use with the {% include ARTIFACT name="interactive" text="anvi'o interactive interface" %} enabled by {% include PROGRAM name="anvi-display-pan" %}.

* Calculate a {% include ARTIFACT name="pan-graph" text="anvi'o pangenome graph" %} based on these sequences using {% include PROGRAM name="anvi-pan-graph" %}

* Step-by-step explanation of graph motives in the {% include ARTIFACT name="interactive" text="anvi'o interactive interface" %} using {% include PROGRAM name="anvi-display-pan-graph" %}

### Dependencies

If your system is properly setup, this {% include PROGRAM name="anvi-self-test" %} command should run without any errors:

``` bash
conda activate anvio-dev
anvi-self-test --suite pangenomics
```

## Get the data

For this workflow we will use four genomes of *Candidatus Lucifugimonas marina* published by Lim et al. in 2023. Two complete genomes (GCA_029593895 & GCA_029593915) and two draft genomes (GCA_029532165, GCA_029532145). First we download all four genomes from the NCBI FTP server.

``` bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/593/915/GCF_029593915.1_ASM2959391v1/GCF_029593915.1_ASM2959391v1_genomic.fna.gz -O GCF_029593915.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/593/895/GCF_029593895.1_ASM2959389v1/GCF_029593895.1_ASM2959389v1_genomic.fna.gz -O GCF_029593895.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/532/145/GCF_029532145.1_ASM2953214v1/GCF_029532145.1_ASM2953214v1_genomic.fna.gz -O GCF_029532145.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/532/165/GCF_029532165.1_ASM2953216v1/GCF_029532165.1_ASM2953216v1_genomic.fna.gz -O GCF_029532165.fna.gz
gzip -d *.gz
```

The current folder structure should look like this.

/path/to/\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── GCF_029532145.fna\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── GCF_029532165.fna\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── GCF_029593895.fna\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└── GCF_029593915.fna

## Order and reorient the draft genomes

In the next step we have to order and reorient the draft genome contigs. We use the program {% include PROGRAM name="anvi-reorient-contigs" %} to order and reorient the contigs based on the longest complete genome present in our folder.

``` bash
anvi-reorient-fasta -f ./
```

The fasta files containing draft genomes are now ready to be used on the pangenomics workflow, if we see some blast result files in the folder.

/path/to/\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── GCF_029532145.fna\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── GCF_029532145_blast\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── GCF_029532165.fna\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── GCF_029532165_blast\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── GCF_029593895.fna\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── GCF_029593895_blast\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└── GCF_029593915.fna

## Running the anvi'o pangenomics workflow

From the workflow we will receive multiple databases that we need for further downstream analysis. A single {% include ARTIFACT name="contigs-db" %} for every genome used, containig informations like genecalls, annotations, etc. The {% include ARTIFACT name="genomes-storage-db" %} a special anvi'o database that stores information about genomes, generated from {% include ARTIFACT name="external-genomes" %}, {% include ARTIFACT name="internal-genomes" %}, or both. And lastly the {% include ARTIFACT name="pan-db" %} created from the {% include ARTIFACT name="genomes-storage-db" %} including all features calculated during the pangenomics analysis. For a more detailed description of these special anvi'o databases please read about it in the [anvi'o pangenomics workflow](https://merenlab.org/2016/11/08/pangenomics-v2/).

``` bash
echo -e 'name\tpath' > fasta.txt
for filename in *.fna; do
    echo -e $(basename $filename | cut -d. -f1)'\t'$(realpath $filename) >> fasta.txt
done
```

fasta.txt:
``` txt
name	path
GCF_029532145	/path/to/GCF_029532145.fna
GCF_029532165	/path/to/GCF_029532165.fna
GCF_029593895	/path/to/GCF_029593895.fna
GCF_029593915	/path/to/GCF_029593915.fna
```

```bash
anvi-run-workflow -w pangenomics --get-default-config config.yaml
```

For the sake of simplicity we will keep most values as they are. The only changes we make are on the bottom of the {% include ARTIFACT name="workflow-config" %} by changing the name of the project and add an output path for the {% include ARTIFACT name="external-genomes" %}.

config.yaml:
``` yaml
{...}
    "project_name": "Candidatus_Lucifugimonas_marina",
    "internal_genomes": "",
    "external_genomes": "external-genomes.txt",
    "sequence_source_for_phylogeny": "",
    "output_dirs": {
        "FASTA_DIR": "01_FASTA",
        "CONTIGS_DIR": "02_CONTIGS",
        "PHYLO_DIR": "01_PHYLOGENOMICS",
        "PAN_DIR": "03_PAN",
        "LOGS_DIR": "00_LOGS"
    },
    "max_threads": "",
    "config_version": "3",
    "workflow_name": "pangenomics"
}
```

``` bash
anvi-run-workflow -w pangenomics -c config.yaml
```

[![First pangenome]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_1.png)]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_1.png){:.center-img .width-80}

/path/to/\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── 00_LOGS\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└── ...\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── 01_FASTA\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└── ...\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── 02_CONTIGS\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── GCF_029532145-contigs.db\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── GCF_029532165-contigs.db\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── GCF_029593895-contigs.db\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── GCF_029593915-contigs.db\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└── ...\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── 03_PAN\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── Candidatus_Lucifugimonas_marina-GENOMES.db\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── Candidatus_Lucifugimonas_marina-PAN.db\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└── ...\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── config.yaml\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── external-genomes.txt\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── fasta.txt\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── GCF_029532145.fna\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── GCF_029532165.fna\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── GCF_029593895.fna\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└── GCF_029593915.fna

## Decide on the genomes subset to use

## Create an anvi'o pangenome graph

``` bash
anvi-pan-graph -p 03_PAN/SAR202_Group_1-PAN.db -g 03_PAN/SAR202_Group_1-GENOMES.db -e external-genomes.txt --output-file SAR202_Group_1-JSON.json
anvi-display-pan-graph -p 03_PAN/SAR202_Group_1-PAN.db -g 03_PAN/SAR202_Group_1-GENOMES.db --pan-graph-json SAR202_Group_1-JSON.json
```

[![First pangenome graph]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_2.png)]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_2.png){:.center-img .width-80}

[![Final pangenome graph]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_3.png)]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_3.png){:.center-img .width-80}
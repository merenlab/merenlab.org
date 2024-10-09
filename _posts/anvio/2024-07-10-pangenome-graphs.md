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

* Download and order a four genome dataset containing complete and draft genomes using anvi'o {% include PROGRAM name="anvi-reorient-contigs" %}

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

Anvi'o might complain in case you haven't set up the databases yet and we strongly advice you to do so before running this tutorial. 

```bash
anvi-setup-scg-taxonomy
anvi-setup-ncbi-cogs
anvi-setup-kegg-data
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

In the next step we have to order and reorient the contigs of the draft genome. We use the program {% include PROGRAM name="anvi-reorient-contigs" %} to for this task, based on the longest complete genome present in our folder.

``` bash
anvi-reorient-fasta -f ./ --prioritize-number-of-contigs
```

The fasta files containing draft genomes are now ready to be used on the pangenomics workflow. For every genome excluding the leading genome we see some blast result files in the folder.

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

A {% include ARTIFACT name="fasta-file" %} is required to run the pangenomics workflow. This file has to be tab-seperated with two columns, the first containing the name of the fasta and the second containig the path to the file. We use the power of bash, to speed things up a little. For more informations about anvi'o workflows in general we really recommend reading [Scaling up your analysis with workflows](https://anvio.org/tutorials/scaling-up/).

``` bash
echo -e 'name\tpath' > fasta.txt
for filename in *.fna; do
    echo -e $(basename $filename | cut -d. -f1)'\t'$(realpath $filename) >> fasta.txt
done
```

The resulting fasta.txt file should look like this.

``` txt
name	path
GCF_029532145	/path/to/GCF_029532145.fna
GCF_029532165	/path/to/GCF_029532165.fna
GCF_029593895	/path/to/GCF_029593895.fna
GCF_029593915	/path/to/GCF_029593915.fna
```

Now the first prerequisite to start anvi'o pangenomics workflow is set-up. The next step is to generate and modify the {% include ARTIFACT name="workflow-config" %}. We can create a standard one with this command.

```bash
anvi-run-workflow -w pangenomics --get-default-config config.yaml
```

For the sake of simplicity we will keep most values as they are. The only changes we make are on the bottom of the {% include ARTIFACT name="workflow-config" %} by changing the name of the project and add an output path for the {% include ARTIFACT name="external-genomes" %}, because we need this file later.

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

Now we can start the anvi'o pangenomics workflow. This will probably take a while and you can reward yourself with a coffee while you wait :coffee: 

``` bash
anvi-run-workflow -w pangenomics -c config.yaml
```

The final command of this section will be to calculate the ANI of our dataset with the program, {% include PROGRAM name="anvi-compute-genome-similarity" %}, it uses various tools such as PyANI to compute average nucleotide identity across your genomes, followed by sourmash to compute mash distance across your genomes. We use it now to add these results as an additional layer data to our pangenome. More informations on ANI and the usage in comparative genomics can be found in the tutorial on pangenomics in the section [Computing the average nucleotide identity for genomes](https://merenlab.org/2016/11/08/pangenomics-v2/#computing-the-average-nucleotide-identity-for-genomes-and-other-genome-similarity-metrics-too). Since our genomes are quite similar it should not take longer than a minute.

```bash
anvi-compute-genome-similarity --external-genomes external-genomes.txt \
                               --program pyANI \
                               --output-dir 04_ANI \
                               --num-threads 6 \
                               --pan-db 03_PAN/Candidatus_Lucifugimonas_marina-PAN.db
```

Good job, you are done! Your folder structure should look similar to this.

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
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── 04_ANI\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└── ...\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── config.yaml\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── external-genomes.txt\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── fasta.txt\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── GCF_029532145.fna\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── GCF_029532165.fna\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── GCF_029593895.fna\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└── GCF_029593915.fna

## Decide on the genomes subset to use

We first run {% include ARTIFACT name="interactive" text="anvi'o interactive interface" %} with {% include PROGRAM name="anvi-display-pan" %} to visualize our pangenome.

``` bash
anvi-display-pan -p 03_PAN/Candidatus_Lucifugimonas_marina-PAN.db \
                 -g 03_PAN/Candidatus_Lucifugimonas_marina-GENOMES.db
```

After you click on the draw button you should see a pangenome that looks somewhat similar to this. We strongly recommend to use Google Chrome to offer you the best possible user experience.

[![First pangenome]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_1.png)]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_1.png){:.center-img .width-80}

In the **Main** tab under **Layers** we first order by **gene_cluster_frequencies**. In **Layer Groups** you checkmark the entry **ANI_full_percentage_identity** and in **Display** you change the bottom four entries to Bar and the Min value to at least 0.97.

[![Settings pangenome]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_2.png)]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_2.png){:.center-img .width-80}

After pressing the draw button again, the resulting pangenome should look like this.

[![Final pangenome]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_3.png)]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_3.png){:.center-img .width-80}

In the newly added red squares we see the ANI between the genomes of the pangenome. Aside from the diagonal which contains a similarity of 100% due to comparing the same genomes with each other we see a second very high sqare on the top left. The complete genome GCF_029593915 shares a very high ANI with the draft genome GCF_029532145. Therefore our initial step in creating a pangenome graph is to use those two genomes.

## Create an anvi'o pangenome graph

```bash
mkdir 05_PANGRAPH
```

``` bash
anvi-pan-graph -e external-genomes.txt \
               -p 03_PAN/Candidatus_Lucifugimonas_marina-PAN.db \
               -g 03_PAN/Candidatus_Lucifugimonas_marina-GENOMES.db \
               -o 05_PANGRAPH/Candidatus_Lucifugimonas_marina-JSON.json \
               -G 'GCF_029593915,GCF_029532145'
```

```bash
anvi-display-pan-graph -p 03_PAN/Candidatus_Lucifugimonas_marina-PAN.db \
                       -g 03_PAN/Candidatus_Lucifugimonas_marina-GENOMES.db \
                       -i 05_PANGRAPH/Candidatus_Lucifugimonas_marina-JSON.json
```

[![First pangenome graph]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_4.png)]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_4.png){:.center-img .width-80}

[![Settings first pangenome graph part 1]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_5_1.png)]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_5_1.png){:.width-40} [![Settings first pangenome graph part 2]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_5_2.png)]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_5_2.png){:.width-40}

[![second pangenome graph]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_6.png)]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_6.png){:.center-img .width-80}

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
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── 04_ANI\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└── ...\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── 05_PANGRAPH\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└── Candidatus_Lucifugimonas_marina-JSON.json \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── config.yaml\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── external-genomes.txt\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── fasta.txt\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── GCF_029532145.fna\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── GCF_029532165.fna\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── GCF_029593895.fna\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└── GCF_029593915.fna

## Expanding the pangenome graph with the remaining genomes

``` bash
anvi-pan-graph -e external-genomes.txt \
               -p 03_PAN/Candidatus_Lucifugimonas_marina-PAN.db \
               -g 03_PAN/Candidatus_Lucifugimonas_marina-GENOMES.db \
               -o 05_PANGRAPH/Candidatus_Lucifugimonas_marina-JSON_2.json
```

``` bash
anvi-display-pan-graph -p 03_PAN/Candidatus_Lucifugimonas_marina-PAN.db \
                       -g 03_PAN/Candidatus_Lucifugimonas_marina-GENOMES.db \
                       -i 05_PANGRAPH/Candidatus_Lucifugimonas_marina-JSON_2.json
```

[![Third pangenome graph]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_7.png)]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_7.png){:.center-img .width-80}

[![Settings third pangenome graph part 1]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_8_1.png)]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_8_1.png){:.width-40} [![Settings third pangenome graph part 2]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_8_2.png)]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_8_2.png){:.width-40}

[![Final pangenome graph]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_9.png)]({{site.url}}/images/anvio/2024-07-10-pangenome-graphs/Figure_9.png){:.center-img .width-80}

## Read from the graph - examples and explanaitions

JUST SOME IDEAS:
- What breaks an assembly? - Mobilome, Transposases
- Why is subclustering GCs worth it? - POSCs, Enolases, Paralogs
- Why do we need complex graphs? - Rearrangement (GCF_029593895)
- Why is there so much similar stuff? - Synteny highly conserved
- This circular view is ugly! - Linear view
- What elso can we get from pangenome graphs? - Datatables and values
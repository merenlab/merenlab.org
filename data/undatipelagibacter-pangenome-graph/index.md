---
layout: page
title: A reproducible workflow for Henoch et al, 2026
modified: 2026-05-11
excerpt: "A bioinformatics workflow for our study on the *Undatipelagibacter* pangenome graph"
comments: true
authors: [alex]
---

<div class="extra-info" markdown="1">

<span class="extra-info-header">Summary</span>

The purpose of this reproducible bioinformatics workflow is to give access to ad-hoc analyses and Python code that underpin our findings in the study, "**Synteny-aware microbial pangenomes reveal blueprints of genomic variation**" by Henoch et al:

[PAPER WILL BE HERE ONCE THE BIORXIV PRE-PRINT IS OUT]

Here is a list of links for quick access to the raw and intermediate data used in our manuscript:

* [FASTA files](https://cloud.uol.de/public.php/dav/files/jSFTXG3cSQMBjYX) for 29 *Undatipelagibacter* genomes.
* [Anvi'o contigs databases](https://cloud.uol.de/public.php/dav/files/7bRpYznDNBedSRk), as annotated [digital microbe](https://www.nature.com/articles/s41597-024-03778-z) files, for *Undatipelagibacter*.
* The Undatipelagibacter pangenome: [genome-storage-db](https://cloud.uol.de/public.php/dav/files/TN2bxBCbAS5DRDJ) and [pan-db](https://cloud.uol.de/public.php/dav/files/ctRp8xRWwaPSnp5).
* The Undatipelagibacter pangenome-graph: [pan-graph-db](https://cloud.uol.de/public.php/dav/files/8eZZYqNrAdXF4TA). You can visualize this pangenome graph interactively at [https://merenlab.org/upg-interactive/](/upg-interactive/).
</div>

{:.notice}
If you have any questions, notice an issue, and/or are unable to find an important piece of information here, please feel free to leave a comment down below, send an e-mail to [us](/people/), or get in touch with us through Discord:

{% include _join-anvio-discord.html %}

## Study description and Introduction

Our study implements a computational workflow to study gene-centric pangenome graphs interactively and applies it to 29 circular isolates from _Undatipelagibacter_ to demonstrate that,


* Synteny-aware pangenome graphs **preserve chromosomal context that conventional pangenomes discard**, revealing an **intricate landscape of variable regions** alternating with backbone stretches across genomes,
* Variable regions are **functionally specialized and coherent units** that are distinct from the genomic backbone, rather than random assortments of genes,
* Genomic variability is distributed as a **structured continuum** rather than a binary partition of hypervariable islands and a static backbone, as evidenced by graph-derived metrics such as the "Composite Variability Score" we have implemented,
* Variable regions differ from one another in **scale, structural topology, functional identity, and evolutionary character**, clustering into distinguishable organizational patterns that likely reflect distinct evolutionary mechanisms,
* Fine-grained, **position-specific sequence divergence gradients** can be detected within otherwise conserved operons (exemplified by a volcano-shaped amino acid identity pattern centered on the Skp), and
* **Deeply conserved variable regions** with consistent functional signatures and genomic contexts are shared across genera within the _Pelagibacterales_, pointing to ancestral hotspots of variation maintained over deep evolutionary timescales.

Our reproducible bioinfromatics workflow picks up from another document related to this study, which [explains how to generate pangenome graphs](https://merenlab.org/tutorials/undatipelagibacter-pangenome-graph/) using the **same set of genomes**. Therefore the output of the reproducible tutorial becomes the input of our reproducible bioinformatics workflow.

The sections below produce the analyses behind each main figure of our paper.

* The top panel of __Figure 1__, which visualizes the pangenome graph, is produced in the section [visualizing the Undatipelagibacter pangenome graph](##visualizing-the-undatipelagibacter-pangenome-graph) section,
* The Sankey diagram at the bottom of __Figure 1__ is produced across the [combining the pangenome graph tables](#combining-the-pangenome-graph-tables-into-one) and [summary statistics](#calculating-summary-statistics-on-the-combined-table) sections,
* __Figure 2__, which visualizes the functional distributions of variable regions, and its associated supplementary panel are produced in the [functional distributions plots and VR/BR comparisons](#functional-distributions-plots-and-vrbr-comparisons) section,
* __Figure 4__, which visualizes the graph-derived metrics and Composite Variability Score, is produced in the [metrics](#metrics-of-the-pangenome-graph) section, and
* __Figure 5__, which visualizes position-wise sequence similarity patterns, is produced in the [position-wise comparisons](#position-wise-sequence-comparisons) section.


## Setting up the stage

Reproduce our study requires a few simple steps to set things up, which will not take more than a few minutes.

This reproducible workflow assumes that you have access to a conda enviornment for the development version of anvi'o (`anvio-dev`), which you can install via [https://anvio.org/install/](https://anvio.org/install/#development-version).

In addition to `anvio-dev`, the reproducible workflow requires a *second conda environment*, since some of the tools used below (such as `holoviews`) are not native to the anvio environment. To keep the two environments separate (so that the anvi'o stack stays intact and isolated from the analysis stack that is only used for downstream plotting and statistics), please run the following commands to generate a second conda enviornment. Running these commands will not take more than a minute on a laptop computer:

```bash
# make sure you are not in the anvi'o environment:
conda deactivate

# create a new environment for the reproducible workflow
conda create -y -n henoch_et_al_2026 -c conda-forge -c bioconda \
        python=3.10 \
        holoviews matplotlib seaborn pandas numpy scipy biopython scikit-learn tqdm
```

At this stage, when you run `conda env list` on your terminal, you should see an output that includes at least the following two items:

```bash
anvio-dev                /[some path to]/miniconda3/envs/anvio-dev
henoch_et_al_2026        /[some path to]/miniconda3/envs/henoch_et_al_2026
```

If that is the case, you are ready to clone the repository of our ad-hoc scripts. For this, you can simply run the following command, which will generate a new directory in your home folder:

```
# go to your home directory:
cd ~

# clone the repostiory:
git clone https://github.com/merenlab/Henoch_et_al_2026_pangenome_graphs.git

# change your work directory to the scripts directory:
cd Henoch_et_al_2026_pangenome_graphs
```

At this stage, your working directory structure should look like this:

```
.
├── README.md
├── 00_SCRIPTS
│   ├── create_all_combined.py
│   ├── functional_distribution_clustering.py
│   ├── metrics_clustering.py
│   ├── similarity_per_position.py
│   └── summary_statistics.py
├── 01_DATA
│   └── 00_README.txt
└── 02_RESULTS
    └── 00_README.txt
```

If that is the case, we are good.

The final step of setting up the stage is to download the files that represent the *Undatipelagibacter* pangenome graph and associated files from our [reproducible tutorial](https://merenlab.org/tutorials/undatipelagibacter-pangenome-graph/). For this, please simply copy-paste these commands into your terminal (while still in the same directory):

```bash
curl -L https://cloud.uol.de/public.php/dav/files/TN2bxBCbAS5DRDJ -o 01_DATA/UNDATIPELAGIBACTER-GENOMES.db
curl -L https://cloud.uol.de/public.php/dav/files/ctRp8xRWwaPSnp5 -o 01_DATA/UNDATIPELAGIBACTER-PAN.db
curl -L https://cloud.uol.de/public.php/dav/files/8eZZYqNrAdXF4TA -o 01_DATA/UNDATIPELAGIBACTER-PAN-GRAPH.db
curl -L https://cloud.uol.de/public.php/dav/files/8snz92oDqJeARDK -o 01_DATA/UNDATIPELAGIBACTER-PAN-GRAPH-SUMMARY.tar.gz
tar -xzf 01_DATA/UNDATIPELAGIBACTER-PAN-GRAPH-SUMMARY.tar.gz -C 01_DATA/ && rm 01_DATA/UNDATIPELAGIBACTER-PAN-GRAPH-SUMMARY.tar.gz

# migrate all DB files in case anvi'o had any database updates in between
# in between
anvi-migrate 01_DATA/*db --migrate-quickly
```

At this stage, your working directory structure should look like this:

```
.
├── README.md
├── 00_SCRIPTS
│   ├── create_all_combined.py
│   ├── functional_distribution_clustering.py
│   ├── metrics_clustering.py
│   ├── similarity_per_position.py
│   └── summary_statistics.py
├── 01_DATA
│   ├── 00_README.txt
│   ├── UNDATIPELAGIBACTER-GENOMES.db
│   ├── UNDATIPELAGIBACTER-PAN-GRAPH.db
│   ├── UNDATIPELAGIBACTER-PAN.db
│   └── UNDATIPELAGIBACTER-PAN-GRAPH-SUMMARY
│       ├── GENESxSYNGCs.txt
│       ├── GENOMES_DIST.newick
│       ├── GENOMES_DIST_MAT.txt
│       ├── REGIONS.txt
│       ├── SYNGCs.txt
│       └── misc_data_items
│           └── default.txt
└── 02_RESULTS
    └── 00_README.txt
```

If that is the case, you are ready.

## Visualizing the _Undatipelagibacter_ pangenome graph

Our reproducicble tutorial [here](https://merenlab.org/tutorials/undatipelagibacter-pangenome-graph/) already details of how `UNDATIPELAGIBACTER-PAN.db` (an anvi'o {% include ARTIFACT name='pan-db' %} artifact) and `UNDATIPELAGIBACTER-PAN-GRAPH.db` (an anvi'o {% include ARTIFACT name='pan-graph-db' %} artifact) are generated from FASTA files. What constitutes the left-upper and right-upper figures of the paper's __Figure 1__ are also the recommended starting point for exploring the data:



You can also visualize the pangenome and pangenome graph by activating `anvio-dev`:

```
conda deactivate
conda activate anvio-dev
```

And running the folowing command to visalize the {% include ARTIFACT name='pan-db' %} in your `01_DATA` directory,

```
anvi-display-pan -g 01_DATA/UNDATIPELAGIBACTER-GENOMES.db \
                 -p 01_DATA/UNDATIPELAGIBACTER-PAN.db
```

which will give you an interactive display for that shows you the 'pangenome' part of the Figure 1,

{% include IMAGE path="images/undatipelagibacter_pangenome.png" width="70" caption="The Undatipelagibacter pangenome" %}

And you can run teh following command to visalize the {% include ARTIFACT name='pan-graph-db' %} in your `01_DATA` directory,

```
anvi-display-pan-graph -g 01_DATA/UNDATIPELAGIBACTER-GENOMES.db \
                       -p 01_DATA/UNDATIPELAGIBACTER-PAN-GRAPH.db
```

which will give you an interactive display for that shows you the 'pangenome graph' part of the Figure 1,

{% include IMAGE path="images/undatipelagibacter_pangenome_graph.png" width="70" caption="The Undatipelagibacter pangenome graph" %}

Using the interactive display you can zoom into any variable region, and inspect individual SynGCs together with their functional annotations, and the genomes that contribute to them.

Backbone SynGCs are colored in blue and variable regions in yellow by default, while the different SynGC types (core, duplication, rearrangement, accessory, singleton, and tRNA/rRNA) carry their own node colors that you can change through the interface. All regions are labeled and enable you to jump directly to specific variable regions of interest (e.g. VR #32, #90, #22, or #180, the four highest-CVS regions in the paper), and the stored states reproduce the exact view used in the figures. From any node, you can pull up the underlying genes, their amino acid sequences, and the functional annotations across all contributing genomes, which is how we drilled into individual VRs throughout the paper.

Similarly the tutorial can be used to generate pangenome graphs from the other _Pelagibacterales_ datasets, that build __Figure 3__.

## Combining the pangenome graph tables into one

For most of the commands below, we will stay in the conda environment `henoch_et_al_2026`. Let's switch to it now:

```bash
conda deactivate
conda activate henoch_et_al_2026
```

For easier downstream analysis we first need to combine two of the pangenome graph tables them together. These tables are not combined from that get go to keep the data as atomary as possible and joining them will create a lot of repetetive information, but this is harmless and makes the analysis a lot easier. The `GENESxSYNGCs.txt` includes information at the gene level and `SYNGCs.txt` at the synteny gene cluster level. By joining them together we get access to all the synteny gene cluster information per gene. The following command generates the joined `all_combined.txt` file.

```bash
python3 00_SCRIPTS/create_all_combined.py \
        -g 01_DATA/UNDATIPELAGIBACTER-PAN-GRAPH-SUMMARY/GENESxSYNGCs.txt \
        -s 01_DATA/UNDATIPELAGIBACTER-PAN-GRAPH-SUMMARY/SYNGCs.txt \
        -d 02_RESULTS/
```

At the same step we also join our definition of simplified COG24 groups to the dataset. In case you want to review our definitions, there is a `cog_functional_groups.txt` that includes the same information as the following table.

|    | __COG24_CATEGORY_ACC__   | __definition__                                                    | __functional_group__                  |
|  0 | A                    | RNA processing and modification                               | Genetic Information Processing    |
|  1 | B                    | Chromatin structure and dynamics                              | Genetic Information Processing    |
|  2 | J                    | Translation, ribosomal structure and biogenesis               | Genetic Information Processing    |
|  3 | K                    | Transcription                                                 | Genetic Information Processing    |
|  4 | L                    | Replication, recombination and repair                         | Genetic Information Processing    |
|  5 | D                    | Cell cycle control, cell division, chromosome partitioning    | Genetic Information Processing    |
|  6 | O                    | Posttranslational modification, protein turnover, chaperones  | Genetic Information Processing    |
|  7 | M                    | Cell wall/membrane/envelope biogenesis                        | Cellular Structure & Organization |
|  8 | N                    | Cell motility                                                 | Cellular Structure & Organization |
|  9 | Z                    | Cytoskeleton                                                  | Cellular Structure & Organization |
| 10 | W                    | Extracellular structures                                      | Cellular Structure & Organization |
| 11 | U                    | Intracellular trafficking, secretion, and vesicular transport | Cellular Structure & Organization |
| 12 | Y                    | Nuclear structure                                             | Cellular Structure & Organization |
| 13 | C                    | Energy production and conversion                              | Metabolism & Energy Production    |
| 14 | E                    | Amino acid transport and metabolism                           | Metabolism & Energy Production    |
| 15 | F                    | Nucleotide transport and metabolism                           | Metabolism & Energy Production    |
| 16 | G                    | Carbohydrate transport and metabolism                         | Metabolism & Energy Production    |
| 17 | H                    | Coenzyme transport and metabolism                             | Metabolism & Energy Production    |
| 18 | I                    | Lipid transport and metabolism                                | Metabolism & Energy Production    |
| 19 | P                    | Inorganic ion transport and metabolism                        | Metabolism & Energy Production    |
| 20 | Q                    | Secondary metabolites biosynthesis, transport and catabolism  | Metabolism & Energy Production    |
| 21 | T                    | Signal transduction mechanisms                                | Cellular Communication & Defense  |
| 22 | V                    | Defense mechanisms                                            | Cellular Communication & Defense  |
| 23 | X                    | Mobilome: prophages, transposons                              | Cellular Communication & Defense  |
| 24 | R                    | General function prediction only                              | Poorly Characterized / Unknown    |
| 25 | S                    | Function unknown                                              | Poorly Characterized / Unknown    |
| 26 | None                 | Function unknown                                              | Poorly Characterized / Unknown    |

The resulting `all_combined.txt` contains the merged information from two of the three summary files together with some extra columns, including the type of the original gene cluster that became a synteny gene cluster and the COG24 group definition, both of which we use later. A snippet of the table looks something like this:

|    | __genome_name__ | __gene_caller_id__ | __source_gene_cluster_id__ | __node_id__ | __region_id__ | __region_type__ | __KEGG_Class_ACC__ | __KEGG_BRITE_ACC__ | __KEGG_Module_ACC__ | __KOfam_ACC__ | __COG24_PATHWAY_ACC__ | __COG24_FUNCTION_ACC__ | __COG24_CATEGORY_ACC__ | __COG24_FUNCTION__ | __node_type__ | __node_x__ | __node_y__ | __genome_count__ | __gene_cluster_type__ | __definition__ | __functional_group__ |
|  0 | HIMB122       |              363 | GC_00000001              | GC_00000001_1 |          65 | backbone      | None             | ko00001          | None              | K03704      | None                | COG1278              | K                    | Cold shock protein, CspA family (CspC) (PDB:1C9O) | duplication |      723 |        0 |             29 | core                | Transcription | Genetic Information Processing |
|  1 | HIMB140       |              355 | GC_00000001              | GC_00000001_1 |          65 | backbone      | None             | ko00001          | None              | K03704      | None                | COG1278              | K                    | Cold shock protein, CspA family (CspC) (PDB:1C9O) | duplication |      723 |        0 |             29 | core                | Transcription | Genetic Information Processing |
|  2 | HIMB1488      |              386 | GC_00000001              | GC_00000001_1 |          65 | backbone      | None             | ko00001          | None              | K03704      | None                | COG1278              | K                    | Cold shock protein, CspA family (CspC) (PDB:1C9O) | duplication |      723 |        0 |             29 | core                | Transcription | Genetic Information Processing |
|  3 | HIMB1491      |              353 | GC_00000001              | GC_00000001_1 |          65 | backbone      | None             | ko00001          | None              | K03704      | None                | COG1278              | K                    | Cold shock protein, CspA family (CspC) (PDB:1C9O) | duplication |      723 |        0 |             29 | core                | Transcription | Genetic Information Processing |
|  4 | HIMB1493      |              332 | GC_00000001              | GC_00000001_1 |          65 | backbone      | None             | ko00001          | None              | K03704      | None                | COG1278              | K                    | Cold shock protein, CspA family (CspC) (PDB:1C9O) | duplication |      723 |        0 |             29 | core                | Transcription | Genetic Information Processing |
| (...) | (...) | (...) | (...) | (...) | (...) | (...) | (...) | (...) | (...) | (...) | (...) | (...) | (...) | (...) | (...) | (...) | (...) | (...) | (...) | (...) | (...) |

### Calculating summary statistics on the combined table

With the `all_combined.txt` file in place we can now generate the summary statistics for the pangenome graph. A dedicated script handles this, and you can run it with the following command once the files are in the right places.

```bash
python3 00_SCRIPTS/summary_statistics.py \
        -c 02_RESULTS/all_combined.txt \
        -d 02_RESULTS/ \
        -pw 9.69 \
        -ph 6.27
```

This will generate four summary tables for you. First the `synteny_gene_clusters_summary.txt` file includes information about the different synteny gene cluster types (`node_type`). We want to make sure to see a nice 100% in the last two columns to make sure that no gene call was left out.

|       | __node_type__     |   __num_syn_cluster__ |   __num_gene_calls__ |   __percent_syn_cluster__ |   __percent_gene_calls__ |
| 0     | accessory     |               855 |             6330 |               22.0816 |              14.3991 |
| 1     | core          |              1186 |            34394 |               30.6302 |              78.2375 |
| 2     | duplication   |               107 |             1065 |               2.76343 |               2.4226 |
| 3     | rearrangement |               350 |              798 |               9.03926 |              1.81525 |
| 4     | singleton     |              1374 |             1374 |               35.4855 |               3.1255 |
| Total |               |              3872 |            43961 |                   100 |                  100 |

Similar to the `synteny_gene_clusters_summary.txt` file we have the second summary table the `gene_clusters_summary.txt` that includes information about the different gene cluster types (`gene_cluster_type`) that created the synteny gene clusters. And again we check for 100% in the last two columns to make sure everything is in place.

|       | __gene_cluster_type__   |   __num_gene_cluster__ |   __num_gene_calls__ |   __percent_gene_cluster__ |   __percent_gene_calls__ |
| 0     | accessory           |               1016 |             7182 |                28.1909 |              16.3372 |
| 1     | core                |               1212 |            35401 |                33.6293 |              80.5282 |
| 2     | singleton           |               1376 |             1378 |                38.1798 |               3.1346 |
| Total |                     |               3604 |            43961 |                    100 |                  100 |

And you probably guessed it, we have a third similar table `regions_summary.txt` that includes the same information but based on the pangenome graphs regions.

|       | __region_type__   |   __num_gene_calls__ |   __num_regions__ |   __num_syn_cluster__ |   __percent_syn_cluster__ |   __percent_gene_calls__ |   __percent_regions__ |
| 0     | backbone      |            35090 |           163 |              1210 |                 31.25 |              79.8208 |           48.9489 |
| 1     | variable      |             8871 |           170 |              2662 |                 68.75 |              20.1792 |           51.0511 |
| Total |               |            43961 |           333 |              3872 |                   100 |                  100 |               100 |

The final summary table `conversion_summary` combines the information of all three other files and shows the conversion from gene cluster to synteny gene cluster to region and the number of gene calls per these conversion.

|       | __gene_cluster_type__   | __node_type__     | __region_type__   |   __num_syn_cluster__ |   __num_gene_cluster__ |   __num_gene_calls__ |   __percent_gene_cluster__ |   __percent_syn_cluster__ |   __percent_gene_calls__ |   __conversion_factor__ |
| 0     | accessory           | accessory     | variable      |               855 |                855 |             6330 |                23.7236 |               22.0816 |              14.3991 |                   1 |
| 1     | accessory           | duplication   | variable      |                52 |                 18 |              170 |               0.499445 |               1.34298 |             0.386706 |             2.88889 |
| 2     | accessory           | rearrangement | variable      |               342 |                143 |              682 |                3.96781 |               8.83264 |              1.55138 |             2.39161 |
| 3     | core                | core          | backbone      |              1183 |               1183 |            34307 |                32.8246 |               30.5527 |              78.0396 |                   1 |
| 4     | core                | core          | variable      |                 3 |                  3 |               87 |              0.0832408 |             0.0774793 |             0.197903 |                   1 |
| 5     | core                | duplication   | Both          |                37 |                 15 |              506 |               0.416204 |              0.955579 |              1.15102 |             2.46667 |
| 6     | core                | duplication   | backbone      |                12 |                  6 |              348 |               0.166482 |              0.309917 |             0.791611 |                   2 |
| 7     | core                | duplication   | variable      |                 2 |                  1 |               37 |              0.0277469 |             0.0516529 |            0.0841655 |                   2 |
| 8     | core                | rearrangement | variable      |                 8 |                  4 |              116 |               0.110988 |              0.206612 |              0.26387 |                   2 |
| 9     | singleton           | duplication   | variable      |                 4 |                  2 |                4 |              0.0554939 |              0.103306 |           0.00909897 |                   2 |
| 10    | singleton           | singleton     | variable      |              1374 |               1374 |             1374 |                38.1243 |               35.4855 |               3.1255 |                   1 |
| Total |                     |               |               |              3872 |               3604 |            43961 |                    100 |                   100 |                  100 |             1.07436 |

The script also generated a figure `conversion_summary_sankey.png` to visualize this exact information as a nice sankey diagram. This figure is the second part of the paper's __Figure 1__.

{% include IMAGE path="images/conversion_summary_sankey.png" width="70" caption="Unedited sankey plot output" %}

## Functional distributions plots and VR/BR comparisons

The __Figure 2__ of the paper shows the functional distributions patterns of the pangenome graph's variable regions. The following script generates these patterns, creates a dendrogram based on the patterns and calculates the Hellinger distance violin plots based on 10,000 subsampling runs.

Under the hood, the script does three things in sequence. (1) for every variable region it counts the genes that fall into each of the five simplified COG24 functional groups defined earlier and normalizes them into a proportion vector, which is what you see as the stacked bar plot in the middle column of __Figure 2__. (2) it clusters these per-region proportion vectors with Ward linkage on Euclidean distances and draws the dendrogram on the left; cutting that dendrogram at the height shown in the figure creates the eight functional clusters that we discussed in the paper. (3) for each VR the script computes the Hellinger distance between its functional proportion vector and the proportion vector of the full backbone, which gives the black dot on the right-hand side of the figure. To assess whether that observed distance is unusual, the script draws 10,000 random samples of backbone genes matched in size to the VR and computes the Hellinger distance between those. The resulting null distribution is plotted as the gray violin behind the dot, and a dot falling outside its own violin indicates a VR whose functional composition differs significantly from what you would expect by simply subsampling the backbone.

```bash
python3 00_SCRIPTS/functional_distribution_clustering.py \
        -c 02_RESULTS/all_combined.txt \
        -d 02_RESULTS/ \
        -pw 6.27 \
        -ph 9.69 \
        -r 01_DATA/UNDATIPELAGIBACTER-PAN-GRAPH-SUMMARY/REGIONS.txt
```

{% include IMAGE path="images/functional_distribution_clustering.png" width="70" caption="Clustering of VRs based on their functional profile" %}

At the same time the script generates the related supplementary figure, showing the difference in distribution patterns between the backbone and variable regions, as well as the different pangenome graph artifacts. This is a useful sanity check that confirms the broad functional divergence between backbone and variable regions.

{% include IMAGE path="images/functional_distribution_by_artifact.png" width="50" %}

## Metrics of the pangenome graph

The following command generates the results we used in the third chapter and especially all figures related to the pangenome graph metrics.

The script reads the per-region Complexity, Expansion, Diversity, Weight and Composite Variability Score values directly from the `REGIONS.txt` summary table. These values are computed inside `anvi-pan-genome-graph` according to the mathematical definitions given in the next subsection. From this table the script produces two complementary outputs. (1) it sorts all variable regions by their CVS and plots the ranked curve shown in the upper-right of __Figure 4__, together with a log curve to see whether the CVS values follow a log like decrease. (2) it clusters variable regions by Complexity and Expansion using Ward linkage and cuts the resulting dendrogram into four groups, which correspond to the four topological categories shown at the bottom of __Figure 4__ ('high complexity / high expansion', 'high complexity / low expansion', 'medium / medium', and 'low / low').

```bash
python3 00_SCRIPTS/metrics_clustering.py \
        -c 02_RESULTS/all_combined.txt \
        -d 02_RESULTS/ \
        -pw 6.27 \
        -ph 6.27 \
        -r 01_DATA/UNDATIPELAGIBACTER-PAN-GRAPH-SUMMARY/REGIONS.txt
```

The script includes, among other calculations, the calculation of the Pearson Correlation Coefficient (r), that we used to test a linear relationship between Complexity and Expansion. The test output will be directly printed in the terminal you used to run the script and will look like this.

```bash
n(VR) = 170

Descriptive stats (VR only):
        complexity_mm_scaled  expansion_mm_scaled
mean                0.090374             0.068984
median              0.045455             0.022727
std                 0.158464             0.132202
min                 0.000000             0.011364
max                 1.000000             1.000000

Spearman rho = +0.527  (p = 1.6e-13)
Pearson  r   = +0.679  (p = 2.61e-24)
Bootstrap 95% CI for rho: [+0.386, +0.644]

--- OLS regression (VR only) ---
========================================================================================
                           coef    std err          t      P>|t|      [0.025      0.975]
----------------------------------------------------------------------------------------
Intercept                0.0178      0.009      2.069      0.040       0.001       0.035
complexity_mm_scaled     0.5664      0.047     11.984      0.000       0.473       0.660
========================================================================================
R-squared = 0.461
```

The script will print the upper part of the paper's __Figure 4__. The lower part was generated from the visualized pangenome graph.

{% include IMAGE path="images/metrics_clustering.png" width="70" %}

### Complexity, expansion, weight, and diversity

This subsection gives the formal definitions of the four metrics we use, along with the supporting notation. You can safely skip ahead if you only want to run the pipeline; the math is here so that anyone who wants to reimplement them can do so directly from the equations.

- Let $\mathbb{G} = \{g_1, g_2, \ldots, g_G\}$ be the set of genomes in the dataset and $G = \|\mathbb{G}\|$ its cardinality.

- Let $\mathbb{H} = \{h_1, h_2, \ldots, h_H\}$ be the set of genomes in which the VR is present and $H = \|\mathbb{H}\|$ its cardinality.

- Let $\mathbb{K} = \{k_1, k_2, \ldots, k_K\}$ be the set of distinct synteny gene clusters in the VR and $K = \|\mathbb{K}\|$ its cardinality.

- Let $\mathbb{P} = \{p_1, p_2, \ldots, p_P\}$ be the set of unique synteny pathways in the VR, ordered such that $\|f(p_1)\| \leq \|f(p_2)\| \leq \cdots \leq \|f(p_P)\|$, where $f(p_i) \subseteq \mathbb{G}$ is the set of genomes in which pathway $p_i$ occurs.

- Let $e_i = \|h(g_i)\|$ be the number of genes contributed by genome $g_i \in \mathbb{G}$.

- Let $n_i = \|t(k_i)\|$ be the number of supporting genomes $t(k_i) \subseteq \mathbb{G}$ containing gene cluster $k_i \in \mathbb{K}$.

**Complexity (C)**: “How many distinct structural realizations (“paths”) occur from the supporting genome?” Answered by estimating the number of events leading to the degree of variation visible in the region. Unique pathways inside the VR are visited in order by the amount of genomes backing it, starting with the ones that are less well represented within the dataset. Every visit of a pathway that includes at least one genome not already seen before with this method, counts as one additional degree of complexity, until the breadth of genomic contribution includes every genome of the dataset. Afterwards we subtract one from the result to set e.g. backbone regions with just a single straight pathway to zero.

$$
\begin{align}
    &X_0 = \emptyset,\quad Y = 0\\
    &\text{for } i = 1, 2, \ldots, P:\\
    &\quad \text{if}\ f(p_i) \not\subseteq X_{i-1}:\\
    &\qquad Y \leftarrow Y + 1,\quad X_i = X_{i-1} \cup f(p_i)\\
    &\quad \text{else}:\\
    &\qquad X_i = X_{i-1}
\end{align}
$$

The result $Y$ is decreased by one to compensate for the fact that a single pathway is always present and then divided by the datasets number of genomes, to weight how high the region can score. Less genomes are less possible pathways.

$$
\begin{align}
    C = \frac{(Y - 1)}{G}
\end{align}
$$

**Expansion (E)**: “How much gene content can be inserted in this variable region?” Answered by calculating the maximum number of newly introduced genes by a single genome in the region.


$$
\begin{align}
    E = \max(e_1, e_2, \ldots, e_G)
\end{align}
$$

**Diversity (D)**: “How heterogeneous is the gene content across supporting genomes?” Answered by describing the overall unevenness of synteny gene cluster prevalence. It first calculates the prevalence proportion mi of every synteny gene cluster in the region and then the reversed population variance, therefore VRs splitting in equal SynGC sizes, in terms of genome contributions, score higher. Reverse population variance is calculated by calculating the maximum possible variance first and subtracting the actual VRs population variance from it.

$$
\begin{align}
    m_i = \frac{n_i}{G}, \quad \bar{m} = \frac{1}{K} \sum_{i=1}^{K} m_i, \quad \bar{v} = \frac{1}{2} \left( \frac{1}{G} + \frac{G}{G} \right)
\end{align}
$$

$$
\begin{align}
    D = \frac{1}{2} \left( \left( \frac{1}{G} - \bar{v} \right)^2 + \left( \frac{G}{G} - \bar{v} \right)^2 \right) - \frac{1}{K} \sum_{i=1}^{K} (m_i - \bar{m})^2
\end{align}
$$

**Weight (W)**: “How high is the variable region's impact?” Answered by describing the potential significance of a given VR within the genomic landscape through the fraction of $H$ and $G$.

$$
\begin{align}
    W = \frac{H}{G}
\end{align}
$$

**Composite Variability Score (CVS)** “How _special_ is the variable region?” Answered by calculating the degree of genomic variation inside a given VR. We use the geometric mean to balance four different terms, requiring higher scores in all metrics to reach a high CVS score.

$$
\begin{align}
    CVS = (C' D' E' W')^{\frac{1}{4}}
\end{align}
$$

For the calculation of the CVS, all terms are normalized according to min-max normalization.

$$
\begin{align}
    Z_{\min} = \min(Z), \quad Z_{\max} = \max(Z)
\end{align}
$$

$$
\begin{align}
    Z' = \frac{(Z - Z_{\min})}{(Z_{\max} - Z_{\min})}
\end{align}
$$

## Position-wise sequence comparisons around the highly divergent Skp gene

The volcano-shaped sequence similarity profile around region #148 in __Figure 5__ is one of the most interesting findings of our paper because it shows that what looks like a single highly variable gene at the pangenome graph level (the _Skp_ gene, which split into twelve distinct SynGCs) is actually the middle of a sequence similarity gradient. To regenerate the data for this we need to drop down to the underlying amino acid sequences and align them column-by-column across all 29 genomes.

We split this analysis into two runs of the same script because of an intermediate step that needs anvi'o programs. The first execution with the `--preprocess` flag scans `all_combined.txt` for the user-supplied region IDs (here `147 148 149 150 151`, which covers the envelope biogenesis operon and a few extra SynGCs on either side) and writes out a list of the conventional gene cluster IDs that those regions correspond to. Anvi'o's gene cluster export program expects GC IDs rather than SynGC or region IDs, so this preprocessing step is what bridges the pangenome graph world back into the conventional pangenome world.

```bash
python3 00_SCRIPTS/similarity_per_position.py \
        -c 02_RESULTS/all_combined.txt \
        -d 02_RESULTS/ \
        -f 147 148 149 150 151 \
        --preprocess
```

Running `anvi-get-sequences-for-gene-clusters` we can then export these gene clusters from the {% include ARTIFACT name='genomes-storage-db' %} and {% include ARTIFACT name='pan-db' %}. The `--split-output-per-gene-cluster` flag produces one FASTA file per gene cluster, each containing the amino acid sequences of every contributing genome aligned within that cluster.

This command requires us to leave the conda environment for `henoch_et_al_2026` and activate `anvio-dev` temporarly:

```bash
# deactivate henoch_et_al_2026 and activate anvio-dev
conda deactivate
conda activate anvio-dev

# get the sequences of interest
anvi-get-sequences-for-gene-clusters -p 01_DATA/UNDATIPELAGIBACTER-PAN.db \
                                     -g 01_DATA/UNDATIPELAGIBACTER-GENOMES.db \
                                     --gene-cluster-ids-file 02_RESULTS/gene_clusters.txt \
                                     --split-output-per-gene-cluster \
                                     -O 02_RESULTS/

# deactivate anvio-dev and activate henoch_et_al_2026
conda deactivate
conda activate henoch_et_al_2026
```

A second execution of the script without the `--preprocess` flag reads these per-cluster FASTA files back in and computes the average pairwise amino acid identity (AAI) at each column of the alignment, using BioPython's pairwise comparison routines. The script then orders the SynGCs along the genomic axis and concatenates their per-position AAI values into the single curve plotted in __Figure 5__. Positions in conserved SynGCs end up at the high end of the curve (~97% on average for the flanking operons), positions in the _Skp_ gene drop to the floor (~40%), and positions in the genes neighboring _Skp_ within the same operon take intermediate values, producing the characteristic volcano pattern.

```bash
python3 00_SCRIPTS/similarity_per_position.py -c 02_RESULTS/all_combined.txt \
                                              -d 02_RESULTS/ \
                                              -pw 6.27 \
                                              -ph 3.23 \
                                              -f 147 148 149 150 151
```

{% include IMAGE path="images/similarity_per_position.png" width="80" caption="The amino acid seuqence identity across graph nodes around Skp" %}


The same procedure can in principle be applied to any other variable region of interest by simply changing the `-f` argument to the desired region IDs, which makes this a generic recipe for inspecting fine-grained sequence variation within and around any region the pangenome graph flags as variable.

The prediction of the protein structures in the upper-part of __Figure 5__ was carried out on an high performance computing cluster with Colabfold and are not part of this reproducible workflow. In case you want to still reproduce these structures, all amino acid sequences for the genes in region #148 are included in `02_RESULTS/position_1401_aa.fa` and instructions on how to run Colabfold can be found [here](https://github.com/sokrypton/colabfold).

An onine, ready-to-use version of Colabfold is also available [here](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb).

## Testing whether Skp divergence mirrors the genome phylogeny

To investigate if the the putative operon that encoded Skp was maintained by recombination, we tested the evolution of Skp with that of the *Undatipelagibacter* genomes. This required a high-resolution phylogenomic analysis of the *Undatipelagibacter* genomes, and a phylogeny of the Skp genes.

### Phylogenomics of *Undatipelagibacter*

To calculate a high-resolution tree for *Undatipelagibacter*, we used the same 165 alphaproteobacterial single-copy core genes used in Freel et al. in the study titled "*[New SAR11 isolate genomes and global marine metagenomes resolve ecologically relevant units within the Pelagibacterales](https://doi.org/10.1038/s41467-025-67043-6)*", following the reproducible workflow [here](https://merenlab.org/data/sar11-phylogenomics/).

Our workflow below download and unpack anvi'o {% include ARTIFACT name="contigs-db" %} files for *Undatipelagibacter*, download the {% include ARTIFACT name="hmm-source" %} for 165 alphaproteobacterial SCGs, annotate *Undatipelagibacter* {% include ARTIFACT name="contigs-db" %} files with them, extract a superalignmet of SCGs, and compute a tree.

---

We downloaded the *Undatipelagibacter* {% include ARTIFACT name="contigs-db" %} files used in our study into our work directory using the following command:


```bash
# make sure you are in the right directory
cd ~/Henoch_et_al_2026_pangenome_graphs

# Download the contigs-db files (which are actualy coming from the
# tutorial section of this reproducible workflow at
# https://merenlab.org/tutorials/undatipelagibacter-pangenome-graph/
curl -l https://cloud.uol.de/public.php/dav/files/7bRpYznDNBedSRk \
     -o 01_DATA/CONTIGS-DBs.tar.gz
```

Next, we unpacked the data file and placed it in a location that follows the structure of our work directory:

```bash
# unpack under 01_DATA
tar -zxvf 01_DATA/CONTIGS-DBs.tar.gz -C 01_DATA

# give it a more descriptive name
mv 01_DATA/CONTIGS-DBs 01_DATA/UNDATIPELAGIBACTER-CONTIGS-DBs

# remove the archive file
rm -rf 01_DATA/UNDATIPELAGIBACTER-CONTIGS-DBs.tar.gz

# migrate the contigs-db files to the latest version
# of anvi'o
anvi-migrate 01_DATA/UNDATIPELAGIBACTER-CONTIGS-DBs/*db \
             --migrate-quickly

# generate an external-genomes file for quick access
anvi-script-gen-genomes-file --input-dir 01_DATA/UNDATIPELAGIBACTER-CONTIGS-DBs/ \
                             -o 01_DATA/UNDATIPELAGIBACTER-CONTIGS-DBs.txt
```

Next, we downloaded the alphaproteobacterial SCGs used to calculate the global tree described in Freel et al. ([2016](https://doi.org/10.1038/s41467-025-67043-6)):

```bash
# Download the hmm-source
curl https://merenlab.org/data/sar11-phylogenomics/files/Alphaproteobacterial_SCGs.tar.gz \
     -o 01_DATA/Alphaproteobacterial_SCGs.tar.gz

# unpack it into the 01_DATA directory
tar -zxvf 01_DATA/Alphaproteobacterial_SCGs.tar.gz -C 01_DATA

# give it a better name
mv 01_DATA/Alphaproteobacterial_SCGs 01_DATA/ALPHAPROTEOBACTERIAL_SCGs

# get rid of the archive
rm -rf 01_DATA/Alphaproteobacterial_SCGs.tar.gz
```

Next, we ran the anvi'o the {% include ARTIFACT name="hmm-source" %} on all *Undatipelagibacter* {% include ARTIFACT name="contigs-db" %} files to identify alphaproteobacterial SCGs in them, which took about 2 seconds per genome:

```bash
for genome in 01_DATA/UNDATIPELAGIBACTER-CONTIGS-DBs/*db
do
    anvi-run-hmms -c $genome \
                  -H 01_DATA/ALPHAPROTEOBACTERIAL_SCGs \
                  --num-threads 8
done
```

Having the `ALPHAPROTEOBACTERIAL_SCGs` as an HMM source (here is making sure using one of the genomes as an example),

```bash
anvi-db-info 01_DATA/UNDATIPELAGIBACTER-CONTIGS-DBs/HIMB122.db

DB Info (no touch)
===============================================
Database Path ................................: 01_DATA/UNDATIPELAGIBACTER-CONTIGS-DBs/HIMB122.db
description ..................................: [Not found, but it's OK]
db_type ......................................: contigs (variant: unknown)
version ......................................: 25


(...)

AVAILABLE HMM SOURCES
===============================================
* 'ALPHAPROTEOBACTERIAL_SCGs' (165 models with 165 hits)
* 'Archaea_76' (76 models with 32 hits)
* 'Bacteria_71' (71 models with 70 hits)
* 'Protista_83' (83 models with 3 hits)
* 'Ribosomal_RNA_16S' (3 models with 1 hit)
* 'Ribosomal_RNA_18S' (1 model with 0 hits)
* 'Transfer_RNAs' (68 models with 30 hits)
```

We next generated a super-alignment for each genome where 165 genes are aligned and concatenated across all genomes (which took about 20 seconds):

```bash
# create a directory in 02_RESULTS to store phylogenomics
# related data
mkdir 02_RESULTS/UNDATIPELAGIBACTER-PHYLOGENOMICS/

anvi-get-sequences-for-hmm-hits --external-genomes 01_DATA/UNDATIPELAGIBACTER-CONTIGS-DBs.txt \
                                -o 02_RESULTS/UNDATIPELAGIBACTER-PHYLOGENOMICS/UNDATIPELAGIBACTER-ALPHASCGs-AA-RAW.fa \
                                --hmm-source ALPHAPROTEOBACTERIAL_SCGs \
                                --return-best-hit \
                                --get-aa-sequences \
                                --concatenate
```

Trimmed columns with excessive gaps:

```
trimal -in  02_RESULTS/UNDATIPELAGIBACTER-PHYLOGENOMICS/UNDATIPELAGIBACTER-ALPHASCGs-AA-RAW.fa \
       -out 02_RESULTS//UNDATIPELAGIBACTER-PHYLOGENOMICS/UNDATIPELAGIBACTER-ALPHASCGs-AA.fa \
       -gt 0.50
```

And ran IQTREE to calculate the final tree for genomes:

```bash
iqtree -s 02_RESULTS/UNDATIPELAGIBACTER-PHYLOGENOMICS//UNDATIPELAGIBACTER-ALPHASCGs-AA.fa \
       -m LG+F+R10 \
       -T 10 \
       -ntmax 25 \
       --alrt 1000 \
       -B 1000
```

Once this is over, we renamed the resulting tree file to a more meaningful name, `UNDATIPELAGIBACTER-ALPHASCGs.newick`,

```
mv 02_RESULTS/UNDATIPELAGIBACTER-PHYLOGENOMICS/UNDATIPELAGIBACTER-ALPHASCGs-AA.fa.treefile \
   02_RESULTS/UNDATIPELAGIBACTER-PHYLOGENOMICS/UNDATIPELAGIBACTER-ALPHASCGs.newick
```

The contents of which looked like this:

```
(HIMB1597:0.0000223763,((((((((((HIMB1770:0.0000224095,HIMB1493:0.0000010069)77.3/97:0.0000222991,
HIMB1518:0.0000215109)100/100:0.0190601071,((((((HIMB1636:0.0000000000,HIMB1526:0.0000000000):0.0000000000,
HIMB1552:0.0000000000):0.0000010069,HIMB1641:0.0000010069)78.8/96:0.0000217393,(HIMB1556:0.0000010069,
HIMB1702:0.0000010069)77.5/95:0.0000218372)100/100:0.0129613967,HIMB1662:0.0118990338)100/100:0.0042673028,
HIMB1723:0.0196013815)100/100:0.0024892878)97.6/98:0.0012401690,HIMB1507:0.0179007591)100/100:0.0018750656,
(HIMB1513:0.0196954887,HIMB1573:0.0208457841)100/88:0.0048378315)39.5/35:0.0007825449,(((HIMB140:0.0207519481,
(HIMB1758:0.0001366421,HIMB1611:0.0002025385)100/100:0.0189204500)92.9/82:0.0017767438,(HIMB1577:0.0000010069,
HIMB1631:0.0000010069)100/100:0.0172580891)71.2/60:0.0011601578,HIMB122:0.0209247152)72.3/32:0.0010367221)100/91:0.0021742263,
(((HIMB1506:0.0180084954,HIMB1488:0.0138845770)100/100:0.0092145862,HIMB1701:0.0288550375)83.6/90:0.0017391483,
HIMB1765:0.0242998974)99.7/99:0.0023063614)100/100:0.0045360204,HIMB1593:0.0162610490)100/100:0.0023453705,
HIMB1491:0.0163187769)26.2/66:0.0021571721,HIMB1685:0.0160545714)100/100:0.0201109324,HIMB1559:0.0000445912);
```

### Recovery and phylogenetics of Skp Genes

To recover Skp gene sequences, we used the anvi'o interactive interface for pangenome graphs.

```bash
anvi-display-pan-graph -p 01_DATA/UNDATIPELAGIBACTER-PAN-GRAPH.db -g 01_DATA/UNDATIPELAGIBACTER-GENOMES.db
```

Once the display showed up, we zoomed into the Skp region:


{% include IMAGE path="images/skp_gene_region.png" width="70" caption="The Skp region" %}

Pressing `Alt` on the keyboard, we selected all Skp genes in a bin:

{% include IMAGE path="images/skp_gene_region_selection.png" width="70" caption="The Skp region selected into a bin" %}

From the bins panel, we clicked on the 'Nodes' column of the Skp bin, and downloaded the FASTA file using the relevant section in the dialog window:

{% include IMAGE path="images/skp_gene_region_download.png" width="70" caption="Downloading the amino acid sequnces for all Skp genes" %}

We moved the downloaded file `UNDATIPELAGIBACTER_Bin_1_GENES_AA.fa` to a better location:

```
mv ~/Downloads/UNDATIPELAGIBACTER_Bin_1_GENES_AA.fa \
   01_DATA/UNDATIPELAGIBACTER_SKP_GENES_AA_UNALIGNED.fa
```

Since the deflines of sequences include gene caller ids, and we only need genome names to be able to compare the phylogeny of the Skp genes to the phylogenomics of the *Undatipelagibacter* genomes, we fixed the deflines using the following `awk` one-liner:

```bash
awk '/^>/{sub(/_[0-9]+$/, "")}1' 01_DATA/UNDATIPELAGIBACTER_SKP_GENES_AA_UNALIGNED.fa > tmp \
    && mv tmp 01_DATA/UNDATIPELAGIBACTER_SKP_GENES_AA_UNALIGNED.fa
```

We created a directory to store Skp phylogeny-related output files:

```bash
mkdir 02_RESULTS/SKP-PHYLOGENETICS
```

We then aligned these sequences using `muscle`,

```bash
muscle -in 01_DATA/UNDATIPELAGIBACTER_SKP_GENES_AA_UNALIGNED.fa \
       -out 02_RESULTS/SKP-PHYLOGENETICS/UNDATIPELAGIBACTER_SKP_GENES_AA_RAW.fa
```

Trimmed excessive gaps,

```bash
trimal -in 02_RESULTS/SKP-PHYLOGENETICS/UNDATIPELAGIBACTER_SKP_GENES_AA_RAW.fa \
       -out 02_RESULTS/SKP-PHYLOGENETICS/UNDATIPELAGIBACTER_SKP_GENES_AA.fa \
       -gt 0.50
```

Ran IQTREE the same way before:

```
iqtree -s 02_RESULTS/SKP-PHYLOGENETICS/UNDATIPELAGIBACTER_SKP_GENES_AA.fa \
       -m LG+F+R10 \
       -T 10 \
       -ntmax 25 \
       --alrt 1000 \
       -B 1000
```

And finally renamed the resulting tree file to a more meaningful name, `UNDATIPELAGIBACTER_SKP_GENES_AA.newick`,

```
mv 02_RESULTS/SKP-PHYLOGENETICS/UNDATIPELAGIBACTER_SKP_GENES_AA.fa.treefile \
   02_RESULTS/SKP-PHYLOGENETICS/UNDATIPELAGIBACTER_SKP_GENES_AA.newick
```

The contents of which looked like this:

```
(HIMB1701:1.1326036100,(((((((HIMB1593:0.8440687974,((HIMB1559:0.0000009947,HIMB1597:0.0000009947)99.9/100:0.5075530029,
HIMB1685:0.4322609089)89.5/91:0.1793674488)74/61:0.0604480327,HIMB1491:0.7283119215)72.9/64:0.0809153224,
HIMB1765:0.6819726624)99.6/100:0.5708229496,(HIMB1577:0.0000009947,HIMB1631:0.0000009947)100/100:0.6940641477)
84.7/65:0.1220357814,(HIMB1723:1.1427802668,((((((HIMB1526:0.0000000000,HIMB1636:0.0000000000):0.0000000000,
HIMB1556:0.0000000000):0.0000000000,HIMB1641:0.0000000000):0.0000000000,HIMB1702:0.0000000000):0.0000009947,
HIMB1552:0.0000009947)92.4/99:0.1120543874,HIMB1662:0.0669078561)99.3/100:0.7530830943)88.5/75:0.2934209924)
1.5/31:0.0500581324,HIMB1507:1.1411752748)74.3/51:0.0772569627,((HIMB1573:1.1078217239,((HIMB1493:0.0000000000,
HIMB1770:0.0000000000):0.0000009947,HIMB1518:0.0000009947)100/100:1.0263878720)1.7/26:0.1530953996,(HIMB1488:0.3343931678,
HIMB1506:0.4578569414)99.9/100:0.6944285056)78.4/30:0.1517070031)64.5/27:0.0990607948,((HIMB1611:0.0000009947,
HIMB1758:0.0000009947)100/100:1.1246136160,(HIMB140:0.5375660781,HIMB1513:0.3438026003)96.4/97:0.3479016466)31.6/42:0.0626461251);
```

### Testing the phylogenetic congruence between Skp genes and genomes

Using the data generated in the previous two steps, we implemented a Phython script ([skp_genome_phylogenetic_congruence.py](https://github.com/merenlab/Henoch_et_al_2026_pangenome_graphs/blob/main/00_SCRIPTS/skp_genome_phylogenetic_congruence.py), available to you under the scripts directory) to test whether Skp divergence mirrors genome ancestry, and if yes, to what extent by asking the following questions:

* Are genome-wide and Skp patristic distances congruent?

* Do per-branch Skp substitutions scale with genome-wide divergence?

Both questions are answered by running the script the following way, which produces the terminal output below with statistics,

```bash
python 00_SCRIPTS/skp_genome_phylogenetic_congruence.py

SKP vs GENOME PHYLOGENETIC CONGRUENCE
===============================================
Genome tree ..................................: 02_RESULTS/UNDATIPELAGIBACTER-PHYLOGENOMICS/UNDATIPELAGIBACTER-ALPHASCGs.newick
Skp gene tree ................................: 02_RESULTS/SKP-PHYLOGENETICS/UNDATIPELAGIBACTER_SKP_GENES_AA.newick
Skp alignment ................................: 02_RESULTS/SKP-PHYLOGENETICS/UNDATIPELAGIBACTER_SKP_GENES_AA.fa
Genomes in analysis ..........................: 29

Genome-wide vs Skp patristic distances
===============================================
Pearson r ....................................: 0.804
Spearman r ...................................: 0.419
Mantel p (49,999 permutations) ...............: 0.0000

Per-branch Skp substitutions vs genome branch length
===============================================
Spearman r ...................................: 0.928
Pearson r ....................................: 0.870

FIGURES
===============================================
PDF ..........................................: 02_RESULTS/skp_genome_phylogenetic_congruence.pdf
PNG ..........................................: 02_RESULTS/skp_genome_phylogenetic_congruence.png
```

And generates the following output figure:

{% include IMAGE path="images/skp_genome_phylogenetic_congruence.png" width="70" caption="Exploring the phylogenetic congruence between Skp genes and genomes that encode them in five panels. Panel a shows genome-wide vs Skp pairwise patristic distance, one point per genome pair. Panel b shows the result of a Mantel permutation null distribution with the observed correlation. Panel c shows a tanglegram that contrasts the genome tree and the Skp gene tree. Panel d shows per-branch Skp substitutions vs genome-wide branch length. And Panel e shows the genome tree with each branch painted by its inferred Skp substitutions" %}

Overall, this result shows that Skp divergence tracks genome ancestry to a large degree, and even though Skp is quite variable across the *Undatipelagibacter* (down to 25% amino-acid identity), it does accumulate substitutions on each lineage in proportion to that lineage's genome-wide divergence given the phylogenomic tree for *Undatipelagibacter* genomes (Spearman r = 0.90), and its pairwise divergences correlate significantly with genome-wide distances (Mantel *p* = 1×10^-4).

Both tests support vertical signal in Skp, and the stronger, metric-robust evidence is per-branch: the number of Skp substitutions mapped onto each genome-tree branch scales tightly with the genome-wide divergence within the branch (Spearman r = 0.90, Pearson r = 0.87, highly concordant results that indicate this is not an artifact of a few long branches). The pairwise distance test between the two trees corroborates this with a significant but moderate monotonic correlation (Spearman ρ = 0.42, Mantel p = 1×10^-4; Pearson r = 0.80). The gap between the two distance metrics reflects the clade's structure: genome-wide distances are strongly bimodal (as we have a tight cluster of near-identical genomes plus a band of divergent pairs (see Panel e)), so both trees agree confidently on deep splits while concording only weakly on the fine ordering among the divergent majority. But the vertical signal for Skp appears to be real and strongest at deeper divergences, rather than a uniform tight tracking of every pairwise relationship.


## Investigating whether Skp divergence matches the divergence of its β-barrel clients

TBD

## Closing notes

If you ran into something that did not behave as described, please open an issue on the [anvi'o GitHub repository](https://github.com/merenlab/anvio/issues) or leave a comment below; both are read regularly and bug reports help us improve the tooling for everyone.

If you would like to **build a pangenome graph from your own genomes** rather than only analyzing the _Undatipelagibacter_ one described here, the [pangenome graph tutorial](https://merenlab.org/tutorials/undatipelagibacter-pangenome-graph/) walks through the same steps starting from raw FASTA files.

Thanks for reading, and feel free to reach out to me with questions, suggestions, or just to share what you discovered in your own pangenome graphs.

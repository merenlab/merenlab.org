---
layout: page
title: Single-cell 'Omics with anvi'o
modified: 2022-03-28
excerpt: "How to use anvi'o to study and investigate single-cell genomics"
categories: [anvio]
authors: [florian]
comments: true
redirect_from: /single-cell-genomics/
---

<div class="extra-info" markdown="1">

<span class="extra-info-header">Summary</span>

**The purpose of this workflow** is to learn how to integrate single-cell genomics data into the anvi'o environment. Here is a list of topics that are covered here:

 * Create a {% include ARTIFACT name="contigs-db" %} with one single amplified genome (SAG).
 * Use single-copy core genes to estimate completion and redundancy.
 * Assign taxonomy using [GTDB](https://gtdb.ecogenomic.org/).
 * Annotate open reading frames with function.
 * Make a phylogenetic tree.
 * Compute a pangenomic analysis.
 * Compare SAGs using average nucleotide identity.
 * Make a phylogenomic tree with the pangenome.

 {:.notice}
 If you have any questions about this exercise, or have ideas to make it better, please feel free to get in touch with the anvi'o community through our Slack channel:

 {% include _join-anvio-slack.html %}
 </div>

 ---
 To reproduce this exercise with your own dataset, you should first follow the instructions [here](/2016/06/26/installation-v2/) to install anvi'o.


## Downloading the data pack

First, open your terminal, `cd` to a working directory, and download the data pack we have stored online for you:

``` bash
curl -L https://figshare.com/ndownloader/files/34713856 \
    -o SCG_workshop.tar.gz
```

Then unpack it, and go into the data pack directory:

``` bash
tar -zxvf SCG_workshop.tar.gz

cd SCG-workshop
```

At this point, if you type `ls` in your terminal, this is what you should be seeing:

```
$ ls
MULTIPLE_SAGs    READ_RECRUITMENT SINGLE_SAG
```

Don't forget to activate your anvi'o envrionment:

```bash
conda activate anvio-7.1
```

## Inspecting a single SAG
We start with the basics: we have a FASTA file containing our assembled SAG.
We can do many things to it, like compute statistics, do gene calling, taxonomy assignment, functional annotation, etc.
And we need to interact with this data and share it with collaborators/reviewers/editors/public.

In the world of anvi'o, we put things in databases, which are basically tables of tables. This allows for great flexibility and very easy sharing. There are different types of databases for different types of data, and for different use-cases.

The most crucial database is the {% include ARTIFACT name="contigs-db" %}. In it we store the contigs' sequences and many features related to those sequences like general metrics (length, N50), gene calls, gene functions and taxonomy information.

Let's move to the right directory:
```bash
cd SINGLE_SAG
```

### How to create a contigs db
To create a contigs db, we need a FASTA file of your contigs, which must have [simple deflines to avoid later issues](https://merenlab.org/2016/06/22/anvio-tutorial-v2/#take-a-look-at-your-fasta-file). We can use {% include PROGRAM name="anvi-script-reformat-fasta" %} to simplify the deflines. The flag `--report-file` will create a TAB-delimited file to map between the new and original deflines.

```bash
# create a directory to store the reformatted fasta file
mkdir -p FASTA

# reformat the fasta file
anvi-script-reformat-fasta DATA/AG-910-K02_contigs.fasta \
                           -o FASTA/AG-910-K02_contigs-fixed.fasta \
                           --simplify-names \
                           --report-file FASTA/AG-910-K02-reformat-report.txt
```

Then we can use the command {% include PROGRAM name="anvi-gen-contigs-database" %} to create the contigs db.
When you run this command, anvi'o will identify open reading frames using [Prodigal](https://doi.org/10.1186/1471-2105-11-119).
```bash
# create a directory for the contigs.db
mkdir -p CONTIGS

# gen the contigs.db
anvi-gen-contigs-database -f FASTA/AG-910-K02_contigs-fixed.fasta \
                          -o CONTIGS/AG-910-K02-contigs.db \
                          -n "AG_910_K02"
```


### Annotate single-copy core genes and ribosomal RNAs
Anvi'o can use Hidden Markov Models ([HMMs](https://en.wikipedia.org/wiki/Hidden_Markov_model)) of your favorite genes to find and annotate open reading frames in your {% include ARTIFACT name="contigs-db" %} with the command {% include PROGRAM name="anvi-run-hmms" %}.
Anvi'o comes with six sets of HMMs for ribosomal RNAs (16S, 23S, 5S, 18S, 28S, 12S), and three sets of single-copy core genes (one for each domain of life). These single-copy core genes are used to compute [completion and redundancy estimates](https://merenlab.org/2016/06/09/assessing-completion-and-contamination-of-MAGs/). There is also an HMM model to identify tRNAs in your SAG.


When you run this command, all sets of HMMs will be used to annotate your database:
```bash
anvi-run-hmms -c CONTIGS/AG-910-K02-contigs.db
```

{:.notice}
If you are interested in tRNAs, you can add the flag `--also-scan-trnas`. But if you don't do this, you can always use the program {% include PROGRAM name="anvi-scan-trnas" %} later.

{:.notice}
`Bacteria_71` and  `Archaea_76` are collections that anvi’o developers curated by taking [Mike Lee’s](https://twitter.com/AstrobioMike) bacterial single-copy core gene collection first released in [GToTree](https://academic.oup.com/bioinformatics/article/35/20/4162/5378708), which is an easy-to-use phylogenomics workflow.
`Protista_83` is [a curated collection](http://merenlab.org/delmont-euk-scgs) of [BUSCO](https://busco.ezlab.org/) by [Tom Delmont](https://twitter.com/tomodelmont).


### Display the contigs' stats

Now that we have created a contigs db and run some annotations, we can have a quick look at our SAG using {% include PROGRAM name="anvi-display-contigs-stats" %}:
```bash
anvi-display-contigs-stats CONTIGS/AG-910-K02-contigs.db
```

On the resulting webpage, we can see some basic metrics about the SAG like the total length, number of contigs, number of genes, N50, etc. We also have some information regarding the HMM hits: number of single-copy core genes found and number of ribosomal RNAs.

{:.notice}
The approximate number of genomes is an estimate based on the frequency of each single-copy core gene. It is mostly useful in metagenomics, where a contigs db contains multiple microbial populations. This estimate is based on the [the mode (or most frequently occurring number) of single-copy core genes](https://merenlab.org/2015/12/07/predicting-number-of-genomes/).
In the context of single-cell genomics, where every contigs db should represent a single population, that value should be 1 for either Archaea, Bacteria or Eukarya. If there is more than one microbial population, i.e. some contamination, that value will be over 1.

To export these metrics as a TAB-delimited file, you can use the flag `--output-file`, and if you don't care about the browser interface you can also add the flag `--report-as-text`:
```bash
anvi-display-contigs-stats CONTIGS/AG-910-K02-contigs.db \
                           --output-file AG-910-K02-metrics.txt \
                           --report-as-text
```
Here is how the output looks like:

|**contigs_db**|**AG_910_K02**|
|:--|:--|
|Total Length|1147435|
|Num Contigs|39|
|Num Contigs > 100 kb|3|
|Num Contigs > 50 kb|7|
|Num Contigs > 20 kb|10|
|Num Contigs > 10 kb|18|
|Num Contigs > 5 kb|26|
|Num Contigs > 2.5 kb|34|
|Longest Contig|293404|
|Shortest Contig|2019|
|Num Genes (prodigal)|1213|
|L50|3|
|L75|8|
|L90|17|
|N50|137780|
|N75|39209|
|N90|11820|
|Archaea_76|24|
|Bacteria_71|55|
|Protista_83|2|
|Ribosomal_RNA_12S|0|
|Ribosomal_RNA_16S|1|
|Ribosomal_RNA_18S|0|
|Ribosomal_RNA_23S|1|
|Ribosomal_RNA_28S|0|
|Ribosomal_RNA_5S|0|
|archaea (Archaea_76)|0|
|eukarya (Protista_83)|0|
|bacteria (Bacteria_71)|1|


### Estimate completion and redundancy
Single-copy core genes (SCGs) are particularly useful and provide a proxy for the completeness and redundancy of a given SAG.
Completeness is estimated based on the number of SCGs found in a SAG. For example, if all bacterial SCGs are found, then the genome's completion is 100%. And if 10/100 SCGs are found in two copies, the the redundancy is 10%.

{:.notice}
Why redundancy and not contamination? The presence of multiple copies of SCGs _could_ be indicative of a contamination in your SAG - ie, the presence of more than one population - but microbial genomes have a bad habit of keeping a few SCGs in multiple copies. So until proven otherwise, anvi'o calls it redundancy and not contamination.

As we already have identified the SCGs with {% include PROGRAM name="anvi-run-hmms" %}, we can now use {% include PROGRAM name="anvi-estimate-genome-completeness" %}:
```bash
anvi-estimate-genome-completeness -c CONTIGS/AG-910-K02-contigs.db
```

You can also use the `--output-file` flag to generate a TAB-delimited file with the results:
```bash
anvi-estimate-genome-completeness -c CONTIGS/AG-910-K02-contigs.db \
                                  --output-file AG-910-K02-completion-redundancy.txt
```

|**bin name**|**domain**|**confidence**|**% completion**|**% redundancy**|**num_splits**|**total length**|
|:--|:--|:--|:--|:--|:--|:--|
|AG-910-K02-contigs|BACTERIA|0.8|77.46|0.00|69|1147435|

### Infer taxonomy
Now that we have identified the SCGs in our SAG, anvi'o can use the Genome Taxonomy Database (GTDB) and Diamond (a fast alternative to NCBI's BLAST) to offer an insight into taxonomy.
In a nutshell, anvi'o downloads both single-copy core genes (SCGs) and the taxonomy (defined by GTDB) of the genome from which they come from. It then builds search databases for a subset of these SCGs based on the domains of life, and finally searches these SCGs from every genome in the anvi'o ecosystem and reports a consensus taxonomy.

To be able to use this function in anvi'o, you must set up SCG taxonomy once on your machine by running the command {% include PROGRAM name="anvi-setup-scg-taxonomy" %}.

Then, to infer taxonomy of the SCGs, we can run the following program:
```bash
anvi-run-scg-taxonomy -c CONTIGS/AG-910-K02-contigs.db
```
This command shows you which and how many SCGs were found in your contigs db, as well as how many were taxonomically annotated. But it does not tell you what is the actual taxonomy estimation for the genome as a whole. There are a few different ways to estimate taxonomy in anvi'o, since a contigs db might represent a complex metagenome or a single genome.

In this case, our contigs db only contains a single SAG so we can simply run the following command, using the `--output-file` flag to store the result:
```bash
anvi-estimate-scg-taxonomy -c CONTIGS/AG-910-K02-contigs.db \
                           --output-file AG-910-K02-taxonomy.txt
```

|**bin_name**|**total_scgs**|**supporting_scgs**|**t_domain**|**t_phylum**|**t_class**|**t_order**|**t_family**|**t_genus**|**t_species**|
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
|AG_910_K02|18|12|Bacteria|Proteobacteria|Alphaproteobacteria|Pelagibacterales|Pelagibacteraceae|Pelagibacter|Pelagibacter sp902575835|

The total number of SCGs represents the number of SCGs found in this SAG that were usable for taxonomy estimation, and the number of supporting SCGs represents the number of SCGs matching the consensus taxonomy.

Only 12 out of 18 SCGs match the consensus taxonomy. We can inspect each SCGs taxonomic affiliation by using the `--debug` flag:

```
anvi-estimate-scg-taxonomy -c CONTIGS/AG-910-K02-contigs.db --debug

Contigs DB ...................................: CONTIGS/AG-910-K02-contigs.db
Metagenome mode ..............................: False

* A total of 18 single-copy core genes with taxonomic affiliations were
successfully initialized from the contigs database 🎉 Following shows the
frequency of these SCGs: Ribosomal_S3_C (1), Ribosomal_S6 (1), Ribosomal_S7 (1),
Ribosomal_S8 (1), Ribosomal_S11 (1), Ribosomal_S20p (1), Ribosomal_L1 (1),
Ribosomal_L2 (1), Ribosomal_L3 (1), Ribosomal_L4 (1), Ribosomal_L6 (1),
Ribosomal_L9_C (1), Ribosomal_L16 (1), Ribosomal_L17 (1), Ribosomal_L20 (1),
Ribosomal_L22 (1), ribosomal_L24 (1), Ribosomal_L27A (1), Ribosomal_S2 (0),
Ribosomal_S9 (0), Ribosomal_L13 (0), Ribosomal_L21p (0).

Hits for AG_910_K02
===============================================
╒════════════════╤════════╤══════════╤══════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╕
│ SCG            │ gene   │ pct id   │ taxonomy                                                                                                                         │
╞════════════════╪════════╪══════════╪══════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╡
│ Ribosomal_L27A │ 643    │ 90.9     │ Bacteria / Proteobacteria / Alphaproteobacteria / Pelagibacterales / Pelagibacteraceae / Pelagibacter / Pelagibacter sp902575835 │
├────────────────┼────────┼──────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ Ribosomal_S20p │ 964    │ 100.0    │ Bacteria / Proteobacteria / Alphaproteobacteria / Pelagibacterales / Pelagibacteraceae / Pelagibacter /                          │
├────────────────┼────────┼──────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ Ribosomal_L20  │ 773    │ 100.0    │ Bacteria / Proteobacteria / Alphaproteobacteria / Pelagibacterales / Pelagibacteraceae / Pelagibacter /                          │
├────────────────┼────────┼──────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ Ribosomal_L6   │ 646    │ 98.9     │ Bacteria / Proteobacteria / Alphaproteobacteria / Pelagibacterales / Pelagibacteraceae / Pelagibacter / Pelagibacter sp902575835 │
├────────────────┼────────┼──────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ Ribosomal_S8   │ 647    │ 98.4     │ Bacteria / Proteobacteria / Alphaproteobacteria / Pelagibacterales / Pelagibacteraceae / Pelagibacter / Pelagibacter sp902575835 │
├────────────────┼────────┼──────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ ribosomal_L24  │ 650    │ 97.4     │ Bacteria / Proteobacteria / Alphaproteobacteria / Pelagibacterales / Pelagibacteraceae / Pelagibacter / Pelagibacter sp902575835 │
├────────────────┼────────┼──────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ Ribosomal_L1   │ 670    │ 99.6     │ Bacteria / Proteobacteria / Alphaproteobacteria / Pelagibacterales / Pelagibacteraceae / Pelagibacter / Pelagibacter sp902575835 │
├────────────────┼────────┼──────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ Ribosomal_L16  │ 654    │ 99.2     │ Bacteria / Proteobacteria / Alphaproteobacteria / Pelagibacterales / Pelagibacteraceae / Pelagibacter /                          │
├────────────────┼────────┼──────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ Ribosomal_S3_C │ 655    │ 99.5     │ Bacteria / Proteobacteria / Alphaproteobacteria / Pelagibacterales / Pelagibacteraceae / Pelagibacter /                          │
├────────────────┼────────┼──────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ Ribosomal_L22  │ 656    │ 100.0    │ Bacteria / Proteobacteria / Alphaproteobacteria / Pelagibacterales / Pelagibacteraceae / Pelagibacter / Pelagibacter sp902575835 │
├────────────────┼────────┼──────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ Ribosomal_L2   │ 658    │ 97.8     │ Bacteria / Proteobacteria / Alphaproteobacteria / Pelagibacterales / Pelagibacteraceae / Pelagibacter / Pelagibacter sp902570695 │
├────────────────┼────────┼──────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ Ribosomal_L4   │ 660    │ 98.9     │ Bacteria / Proteobacteria / Alphaproteobacteria / Pelagibacterales / Pelagibacteraceae / Pelagibacter / Pelagibacter sp902575835 │
├────────────────┼────────┼──────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ Ribosomal_L3   │ 661    │ 100.0    │ Bacteria / Proteobacteria / Alphaproteobacteria / Pelagibacterales / Pelagibacteraceae / Pelagibacter / Pelagibacter sp902575835 │
├────────────────┼────────┼──────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ Ribosomal_S7   │ 664    │ 100.0    │ Bacteria / Proteobacteria / Alphaproteobacteria / Pelagibacterales / Pelagibacteraceae / Pelagibacter /                          │
├────────────────┼────────┼──────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ Ribosomal_L9_C │ 444    │ 98.6     │ Bacteria / Proteobacteria / Alphaproteobacteria / Pelagibacterales / Pelagibacteraceae / Pelagibacter / Pelagibacter sp902575835 │
├────────────────┼────────┼──────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ Ribosomal_L17  │ 637    │ 100.0    │ Bacteria / Proteobacteria / Alphaproteobacteria / Pelagibacterales / Pelagibacteraceae / Pelagibacter / Pelagibacter sp902575835 │
├────────────────┼────────┼──────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ Ribosomal_S6   │ 446    │ 100.0    │ Bacteria / Proteobacteria / Alphaproteobacteria / Pelagibacterales / Pelagibacteraceae / Pelagibacter / Pelagibacter sp902575835 │
├────────────────┼────────┼──────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ Ribosomal_S11  │ 639    │ 100.0    │ Bacteria / Proteobacteria / Alphaproteobacteria / Pelagibacterales / Pelagibacteraceae / Pelagibacter / Pelagibacter sp902575835 │
├────────────────┼────────┼──────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ CONSENSUS      │ --     │ --       │ Bacteria / Proteobacteria / Alphaproteobacteria / Pelagibacterales / Pelagibacteraceae / Pelagibacter / Pelagibacter sp902575835 │
╘════════════════╧════════╧══════════╧══════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╛

Estimated taxonomy for "AG_910_K02"
===============================================
╒════════════╤══════════════╤═══════════════════╤══════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╕
│            │   total_scgs │   supporting_scgs │ taxonomy                                                                                                                         │
╞════════════╪══════════════╪═══════════════════╪══════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╡
│ AG_910_K02 │           18 │                12 │ Bacteria / Proteobacteria / Alphaproteobacteria / Pelagibacterales / Pelagibacteraceae / Pelagibacter / Pelagibacter sp902575835 │
╘════════════╧══════════════╧═══════════════════╧══════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════╛
```

Now we can see that the 6 SCGs that were not matching the consensus `Pelagibacter sp902575835` were matching to the genus `Pelagibacter` or to `Pelagibacter sp902570695`.

Currently, anvi'o's strategy for estimating taxonomy is limited to Archaea and Bacteria, and anvi'o will not estimate taxonomy for a SAG with low completion based on the number of SCGs. It is not a perfect way to estimate taxonomy but it is very fast.

### Functional annotations
There are multiple databases out there that we can use to investigate gene functions. Each of these databases fills a specific need and none of them is perfect. A database like [NCBI's clusters of orthologous groups](https://www.ncbi.nlm.nih.gov/research/cog-project/) (COG) will give general annotations, while the [KEGG database](https://www.genome.jp/kegg/) offers a systematic approach with metabolism-related genes, and others like [Pfam](http://pfam.xfam.org/) focus on protein families. Anvi'o can use these three databases to annotate open reading frames, but you can also import any functional annotations not from these sources.

{:.notice}
**You can import your own functional annotations into anvi'o.** Let's say that you are not happy with the functional annotations provided with anvi'o and you really want to look for carbohydrate-active enzymes using the [CAZy](http://www.cazy.org/) database. Well, you could export the gene sequences (DNA or amino acid) using {% include PROGRAM name="anvi-get-sequences-for-gene-calls" %}, annotate your genes with CAZy (or your favorite tool and/or database), and import the result back into the contigs db with {% include PROGRAM name="anvi-import-functions" %}.

For this analysis, we will use COG and KEGG annotations. To be able to run these annotations, we first need to download and set up the COG and KEGG databases locally using {% include PROGRAM name="anvi-setup-ncbi-cogs" %} and {% include PROGRAM name="anvi-setup-kegg-kofams" %}. Like we did for the SCG taxonomy data, we only need to run these commands once, after which the data will be stored on our machines.

If you have these two databases set up, you can annotate your contigs db like this:

```bash
anvi-run-ncbi-cogs -c CONTIGS/AG-910-K02-contigs.db -T 4
anvi-run-kegg-kofams -c CONTIGS/AG-910-K02-contigs.db -T 4
```

If you want to export the functional annotations from the database as TAB-delimited files, you can use {% include PROGRAM name="anvi-export-functions" %}. For this command you will have to specify which annotation source you want to use. If you are unsure what functional annotations were carried out on your contigs db, you can use the flag `--list-annotation-sources`:

```bash
anvi-export-functions -c CONTIGS/AG-910-K02-contigs.db \
                      --list-annotation-sources


FUNCTIONAL ANNOTATION SOURCES FOUND
===============================================
* COG20_PATHWAY
* COG20_FUNCTION
* KEGG_Module
* KEGG_Class
* KOfam
* COG20_CATEGORY
```

Then we can export the COG annotations:
```bash
anvi-export-functions -c CONTIGS/AG-910-K02-contigs.db \
                      --annotation-sources COG20_FUNCTION \
                      --output-file AG-910-K02-COG20_FUNCTION.txt
```
And the output file should be a table that looks like this:

|**gene_callers_id**|**source**|**accession**|**function**|**e_value**|
|:--|:--|:--|:--|:--|
|0|COG20_FUNCTION|COG4091|Predicted homoserine dehydrogenase, contains C-terminal SAF domain (PDB:3UPL)|9.9e-179|
|1|COG20_FUNCTION|COG4392|Branched-chain amino acid transport protein (AzlD2)|1.2e-21|
|2|COG20_FUNCTION|COG1296|Predicted branched-chain amino acid permease (azaleucine resistance) (AzlC)|2.4e-59|
|3|COG20_FUNCTION|COG2105|Predicted gamma-glutamylamine cyclotransferase YtfP, GGCT/AIG2-like family (YtfP) (PDB:1V30) (PUBMED:16754964;20110353)|6.5e-12|
|4|COG20_FUNCTION|COG1058|ADP-ribose pyrophosphatase domain of DNA damage- and competence-inducible protein CinA (CinA) (PDB:4CT9) (PUBMED:25313401)|5.9e-67|
|\-\-|\-\-|\-\-|\-\-|\-\-|
|1209|COG20_FUNCTION|COG0396|Fe-S cluster assembly ATPase SufC (SufC) (PDB:2D3W)|2e-91|
|1210|COG20_FUNCTION|COG0719|Fe-S cluster assembly scaffold protein SufB (SufB) (PDB:1VH4)|2.2e-204|
|1211|COG20_FUNCTION|COG3808|Na+ or H+-translocating membrane pyrophosphatase (OVP1) (PDB:4A01) (PUBMED:11342551)|7.7e-281|

### Estimate KEGG module completion

What makes the KEGG database unique is its ability to contextualize annotations into functional modules and pathways. We can use the program {% include PROGRAM name="anvi-estimate-metabolism" %} to estimate the completion of KEGG modules for our SAG:
```bash
anvi-estimate-metabolism -c CONTIGS/AG-910-K02-contigs.db
```

It will generate a TAB-delimited file where each row represents a KEGG module, its name, the associated KOs, the completeness of the module in the SAG, and more:

```
unique_id	genome_name	kegg_module	module_name	module_class	module_category	module_subcategory	module_definition	module_completeness	module_is_complete	kofam_hits_in_module	gene_caller_ids_in_module	warnings
0	AG_910_K02	M00001	Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate	Pathway modules	Carbohydrate metabolism	Central carbohydrate metabolism	"(K00844,K12407,K00845,K00886,K08074,K00918) (K01810,K06859,K13810,K15916) (K00850,K16370,K21071,K00918) (K01623,K01624,K11645,K16305,K16306) K01803 ((K00134,K00150) K00927,K11389) (K01834,K15633,K15634,K15635) K01689 (K00873,K12406)"	0.3	False	K00134,K00927,K01623	1045,1046,1047	None
1	AG_910_K02	M00002	Glycolysis, core module involving three-carbon compounds	Pathway modules	Carbohydrate metabolism	Central carbohydrate metabolism	"K01803 ((K00134,K00150) K00927,K11389) (K01834,K15633,K15634,K15635) K01689 (K00873,K12406)"	0.3333333333333333	False	K00134,K00927	1045,1046	None
2	AG_910_K02	M00003	Gluconeogenesis, oxaloacetate => fructose-6P	Pathway modules	Carbohydrate metabolism	Central carbohydrate metabolism	"(K01596,K01610) K01689 (K01834,K15633,K15634,K15635) K00927 (K00134,K00150) K01803 ((K01623,K01624,K11645) (K03841,K02446,K11532,K01086,K04041),K01622)"	0.375	False	K00134,K00927,K01623	1045,1046,1047	None
3	AG_910_K02	M00307	Pyruvate oxidation, pyruvate => acetyl-CoA	Pathway modules	Carbohydrate metabolism	Central carbohydrate metabolism	"((K00163,K00161+K00162)+K00627+K00382-K13997),K00169+K00170+K00171+(K00172,K00189),K03737"	1.0	True	K00163,K00163,K00382,K00627	761,230,467,762	None
4	AG_910_K02	M00009	Citrate cycle (TCA cycle, Krebs cycle)	Pathway modules	Carbohydrate metabolism	Central carbohydrate metabolism	"(K01647,K05942) (K01681,K01682) (K00031,K00030) (K00164+K00658+K00382,K00174+K00175-K00177-K00176) (K01902+K01903,K01899+K01900,K18118) (K00234+K00235+K00236+K00237,K00239+K00240+K00241-(K00242,K18859,K18860),K00244+K00245+K00246-K00247) (K01676,K01679,K01677+K01678) (K00026,K00025,K00024,K00116)"	0.875	True	K00024,K00031,K00164,K00239,K00240,K00241,K00242,K00382,K00658,K01679,K01681,K01902,K01903	472,390,469,475,474,477,476,467,468,369,156,470,471	None
5	AG_910_K02	M00010	Citrate cycle, first carbon oxidation, oxaloacetate => 2-oxoglutarate	Pathway modules	Carbohydrate metabolism	Central carbohydrate metabolism	"(K01647,K05942) (K01681,K01682) (K00031,K00030)"	0.6666666666666666	False	K00031,K01681	390,156	None
6	AG_910_K02	M00011	Citrate cycle, second carbon oxidation, 2-oxoglutarate => oxaloacetate	Pathway modules	Carbohydrate metabolism	Central carbohydrate metabolism	"(K00164+K00658+K00382,K00174+K00175-K00177-K00176) (K01902+K01903,K01899+K01900,K18118) (K00234+K00235+K00236+K00237,K00239+K00240+K00241-(K00242,K18859,K18860),K00244+K00245+K00246-K00247) (K01676,K01679,K01677+K01678) (K00026,K00025,K00024,K00116)"	1.0	True	K00024,K00164,K00239,K00240,K00241,K00242,K00382,K00658,K01679,K01902,K01903	472,469,475,474,477,476,467,468,369,470,471	None
7	AG_910_K02	M00004	Pentose phosphate pathway (Pentose phosphate cycle)	Pathway modules	Carbohydrate metabolism	Central carbohydrate metabolism	"(K13937,((K00036,K19243) (K01057,K07404))) K00033 K01783 (K01807,K01808) K00615 K00616 (K01810,K06859,K13810,K15916)"	0.5714285714285714	False	K00615,K00616,K01783,K01808	1044,1053,550,601	None
8	AG_910_K02	M00007	Pentose phosphate pathway, non-oxidative phase, fructose 6P => ribose 5P	Pathway modules	Carbohydrate metabolism	Central carbohydrate metabolism	"K00615 (K00616,K13810) K01783 (K01807,K01808)"	1.0	True	K00615,K00616,K01783,K01808	1044,1053,550,601	None
```

By default, anvi'o labels a module as 'present' if its completeness is > 0.75 (at least 3 out of 4 KOs must be present). We are using SAGs, which are inherently incomplete, so be very critical about this completeness threshold (which you can change using the flag `--module-completion-threshold`).

Now we have gone through all the basics of analyzing a single SAG in anvi'o. However, in real life we often have many SAGs that we want to analyze all at once. Now we will learn how to do this.

## Working with multiple SAGs

Let's step up our game and repeat the above analyses for the 228 SAGs from the AG-910 sample!

You can change your current working directory to `MULTIPLE_SAGs`:
```bash
cd ../MULTIPLE_SAGs
```

And this is what you should see in it:
```bash
$ ls
ADDITIONAL_DATA DATA
```

With 226 contigs dbs in the `DATA` directory.

{:.notice}
All of the contigs dbs are already provided and were annotated in the same way as shown above with AG_910_K02. The reason why we already did all of this is pretty simple: it takes a lot of time and/or computer resource. But now that you know how to generate and annotate one contigs db, you could do this for many dbs at a time - for instance, by using bash loops. But loops would not be the most effective way to do it. If you have access to a more powerful machine like a computing cluster, you can check out [the anvi'o workflows](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/). Basically, you can use {% include PROGRAM name="anvi-run-workflow" %} to automatically run all the steps to generate a set of contigs dbs, a full metagenomics workflow (including assembly, short-read mapping, anvi'o profiles and more), pangenomics, phylogenomics and more. It uses [snakemake](https://snakemake.readthedocs.io/en/stable/), a very powerful workflow management system.

### General metrics

We've previously used {% include PROGRAM name="anvi-display-contigs-stats" %} to export some information about a single SAG like its length, number of contigs, number of single-copy core genes and more. If you look at the help menu of that command, you will see that it can accept as many contigs dbs as you wish.

Let's run it using all of the contigs dbs available in the datapack:
```bash
anvi-display-contigs-stats --output-file stats.txt DATA/*contigs.db
```

{% include IMAGE path="images/display-contigs-stats.png" width=80 %}

At this point, the interactive interface includes a table with nearly 200 columns, which is a bit overwhelming. To make more sense of the numbers, you can click on `Colorify output` to visually appreciate high and low values in each row.

With the `--output-file` flag, you will save that table as a TAB-delimited file and import it in your favorite table-eating software like R, Python and others (we also use Excel sometimes). You will be able to compute the total size distribution of all the SAGs, calculate the average number of contigs, and select SAGs based on your favorite metrics!


### Completeness and redundancy

While {% include PROGRAM name="anvi-display-contigs-stats" %} was very happy to take as many contigs dbs as you wanted directly as parameters in the command, it is not the case with other programs in anvi'o.

Here is what you can find in the help manual of {% include PROGRAM name="anvi-estimate-genome-completeness" %}:
```
MANDATORY INPUT OPTION #2
  Or you can initiate this with an external genomes file.

  -e FILE_PATH, --external-genomes FILE_PATH
                        A two-column TAB-delimited flat text file that lists
                        anvi'o contigs databases. The first item in the header
                        line should read 'name', and the second should read
                        'contigs_db_path'. Each line in the file should
                        describe a single entry, where the first column is the
                        name of the genome (or MAG), and the second column is
                        the anvi'o contigs database generated for this genome.
                        (default: None)
```

You can use an {% include ARTIFACT name="external-genomes" %} file containing the names and paths to your contigs dbs.
While you _can_ manually generate that file, there is a very convenient command in anvi'o that will do it for you: {% include PROGRAM name="anvi-script-gen-genomes-file" %}.

Let's use it:
```bash
anvi-script-gen-genomes-file --input-dir DATA \
                             -o external-genomes.txt
```

And here is how the file should look like:

```bash
$ head -n 10 external-genomes.txt
name	contigs_db_path
AG_910_A01	/path/to/MULTIPLE_SAGs/DATA/AG_910_A01-contigs.db
AG_910_A02	/path/to/MULTIPLE_SAGs/DATA/AG_910_A02-contigs.db
AG_910_A03	/path/to/MULTIPLE_SAGs/DATA/AG_910_A03-contigs.db
AG_910_A04	/path/to/MULTIPLE_SAGs/DATA/AG_910_A04-contigs.db
AG_910_A06	/path/to/MULTIPLE_SAGs/DATA/AG_910_A06-contigs.db
AG_910_A10	/path/to/MULTIPLE_SAGs/DATA/AG_910_A10-contigs.db
AG_910_A11	/path/to/MULTIPLE_SAGs/DATA/AG_910_A11-contigs.db
AG_910_A13	/path/to/MULTIPLE_SAGs/DATA/AG_910_A13-contigs.db
AG_910_A14	/path/to/MULTIPLE_SAGs/DATA/AG_910_A14-contigs.db
```

We can now use {% include PROGRAM name="anvi-estimate-genome-completeness" %} with the {% include ARTIFACT name="external-genomes" %} file instead of a single contigs db:

```bash
anvi-estimate-genome-completeness -e external-genomes.txt \
                                  --output-file completion-redundancy.txt
```

<div class="extra-info" markdown="1">

<span class="extra-info-header">Use bash to get the five most complete SAGs</span>

Let's imagine that you are interested in the top 5 most complete SAGs. You can use the bash command `sort` and `head` to quickly create and store a list of these SAGs' names in a file.

For `sort`, the `n` flag stands for numerical, the `r` flag is for reverse sorting (descending order), and the `k` flag is for specifying which column you want to do the sorting on. The completeness estimate is in the 4th column in `completion-redundancy.txt`, so we use `k4`.
Then, we can use `head` to select only the first 5 lines.

```bash
sort -nrk4 completion-redundancy.txt | head -n 5

AG_910_D13	BACTERIA	1.0	94.37	0.00	73	1117975
AG_910_M21	BACTERIA	1.0	91.55	0.00	48	910570
AG_910_L17	BACTERIA	0.9	90.14	0.00	58	1089091
AG_910_C02	BACTERIA	0.9	90.14	0.00	55	936322
AG_910_B21	BACTERIA	0.9	90.14	2.82	114	1930694
```

If you want to store the list of SAG names in a file, you can use `cut` to select only the first column:
```bash
sort -nrk4 completion-redundancy.txt | head -n 5 | cut -f 1 > five_most_complete_SAGs.txt
```

</div>

### Taxonomy

Our next step is to investigate the taxonomy of all the SAGs. We could be tempted to use {% include PROGRAM name="anvi-estimate-scg-taxonomy" %} with an external genomes file. But as of anvi'o v7.1, this command cannot use an external genomes file. No one is perfect, but you can find and shame an anvi'o developer by using the slack channel, or [anvi'o github](https://github.com/merenlab/anvio) for a feature request like this.

If you really want the taxonomy of all your SAGs (and it is your right), there is a way with the power of bash!
First of all, you want to loop through every SAG and run {% include PROGRAM name="anvi-estimate-scg-taxonomy" %}.
You can use the `while` loop with the external genomes file. Basically, you create two variables that correspond to the columns of your file. The `if` line is there to skip the header of the external genomes file.

```bash
# create a directory to store all the outputs
mkdir -p TAXONOMY

# going through each SAGs at a time
# it will take a few min
while read name path
do
  if [ "$name" == "name" ]; then continue; fi
  anvi-estimate-scg-taxonomy -c $path -o TAXONOMY/$name.txt
done < external-genomes.txt
```

Once all you have all the output, you can merge them with a similar `while` loop.

```bash
# take the header line and create the general output
head -n 1 TAXONOMY/AG_910_A01.txt > taxonomy.txt

# append each SAG's output to the final file
while read name path; do
  if [ "$name" == "name" ]; then continue; fi
  sed 1d TAXONOMY/$name.txt >> taxonomy.txt
done < external-genomes.txt
```

{:.notice}
`sed 1d` is used to remove the first line (table header).

Now you have the final taxonomy output for all your SAGs! Are you interested in *Pelagibacter*? You can search for all the SAGs with a *Pelagibacter* assignment:

```bash
$ grep Pelagibacter taxonomy.txt
AG_910_A06	19	9	Bacteria	Proteobacteria	Alphaproteobacteria	Pelagibacterales	Pelagibacteraceae	Pelagibacter	Pelagibacter sp902567045
AG_910_A11	1	1	Bacteria	Proteobacteria	Alphaproteobacteria	Pelagibacterales	Pelagibacteraceae	Pelagibacter	Pelagibacter sp902612325
AG_910_A14	4	4	Bacteria	Proteobacteria	Alphaproteobacteria	Pelagibacterales	AG-422-B15	AG-422-B15	None
AG_910_A15	6	4	Bacteria	Proteobacteria	Alphaproteobacteria	Pelagibacterales	Pelagibacteraceae	Pelagibacter	Pelagibacter sp902573345
AG_910_A17	17	17	Bacteria	Proteobacteria	Alphaproteobacteria	Pelagibacterales	Pelagibacteraceae	Pelagibacter	None
---
```


### Phylogenetics

Another way to investigate the taxonomic landscape of all these SAGs is to create a phylogenetic tree using a marker gene.
You could use any gene in one of the HMM models that you ran earlier.

For this tutorial, we will use the 16S rRNA gene and make a phylogenetic tree. The command {% include PROGRAM name="anvi-get-sequences-for-hmm-hits" %} will extract the DNA sequence of a given gene (or all sequences for a set of genes if you want to do some phylogenomics). It will also create an alignment when you use the flag `--concatenate-genes`.

```bash
# make a directory to store the results
mkdir -p PHYLOGENETICS_16S_rRNA

# run the command
anvi-get-sequences-for-hmm-hits -e external-genomes.txt \
                                --hmm-source Ribosomal_RNA_16S \
                                --gene-names 16S_rRNA_bac \
                                --concatenate-genes \
                                --return-best-hit \
                                -o PHYLOGENETICS_16S_rRNA/Ribosomal_RNA_16S_bac.fa
```

Well, there are only 61 bacterial 16S rRNA among the 226 SAGs:
```
Sources ......................................: Ribosomal_RNA_16S
Hits .........................................: 67 HMM hits for 1 source(s)
Genes of interest ............................: 16S_rRNA_bac
Filtered hits ................................: 61 hits remain after filtering for 1 gene(s)
```

Keep in mind that no marker gene or single copy core gene was found in all SCGs and therefore the phylogenetic analysis will be limited to a subset of your dataset.

Now that you have an aligned FASTA file of the bacterial 16S rRNA, you can use the command {% include PROGRAM name="anvi-gen-phylogenomic-tree" %} to genereate a Newick tree file. It uses [FastTree](http://www.microbesonline.org/fasttree/) by default, but you can always use your favorite phylogenetic tree software to make a Newick tree.

The command is quite simple:
```bash
anvi-gen-phylogenomic-tree -f PHYLOGENETICS_16S_rRNA/Ribosomal_RNA_16S_bac.fa \
                           -o PHYLOGENETICS_16S_rRNA/Ribosomal_RNA_16S_bac-tree.txt
```

It is now time to use the interactive interface of anvi'o! To do that, we will use the `--manual` mode of `anvi-interactive`. Anvi'o will automatically generate an empty {% include ARTIFACT name="profile-db" %} to store any settings that you change while working in the interface.

```bash
anvi-interactive -p PHYLOGENETICS_16S_rRNA/phylogenetic-profile-Ribosomal_RNA_16S_bac.db \
                 -t PHYLOGENETICS_16S_rRNA/Ribosomal_RNA_16S_bac-tree.txt \
                 --title "Ribosomal_RNA_16S_bac" \
                 --manual
```

You should see this after clicking the `Draw` button:

{% include IMAGE path="images/phylogenetic-tree-raw.png" width=80 %}

{:.notice}
Trees from FastTree are not rooted. You can use an external software to root your tree, or right-click a branch in the interactive interface and select "reroot the tree/dendrogram here".

This tree is not very informative without the taxonomic assignment that you did earlier. Let's exit the interface by closing the browser tab and press `CTRL+C` in the terminal.

To add additional data for each `item` in the interface (here the items are SAGs), you can use the flag `-A`:

```bash
anvi-interactive -p PHYLOGENETICS_16S_rRNA/phylogenetic-profile-Ribosomal_RNA_16S_bac.db \
                 -t PHYLOGENETICS_16S_rRNA/Ribosomal_RNA_16S_bac-tree.txt \
                 --title "Ribosomal_RNA_16S_bac" \
                 --manual \
                 -A taxonomy.txt
```

The interactive interface is quite modular and you can tinker around to make that tree look better.
Homework for you: try to create a similar figure:

{% include IMAGE path="images/phylogenetic-tree-annotated.png" width=80 %}

{:.notice}
If you open the `Mouse` panel on the right side and then hover your mouse over the figure, you can see the taxonomic annotation in the panel.

From this phylogenetic tree (which only includes a subset of our SAGs), and the taxonomy output, we can see that there is a large proportion of Alphaproteobacteria with mostly Pelagibacterales and HIMB59 genomes.

We have just started to understand the microbial diversity of all the SAGs from that one sample. Of course, if you wanted to do a more comprehensive analysis, you could include SAGs from other samples, or even some reference genomes. You could also use other marker genes.

<div class="extra-info" markdown="1">

<span class="extra-info-header">Distribution of HMM hits across SAGs</span>

You can use the command {% include PROGRAM name="anvi-script-gen-hmm-hits-matrix-across-genomes" %} to check the distribution of genes from the HMM collections.

First, you can list the HMMs sources:
```bash
anvi-script-gen-hmm-hits-matrix-across-genomes -e external-genomes.txt \
                                               --output-file hmm_distribution_matrix.txt \
                                               --list-hmm-sources

HMM SOURCES COMMON TO ALL 226 GENOMES
===============================================
* Archaea_76 [type: singlecopy] [num genes: 976]
* Bacteria_71 [type: singlecopy] [num genes: 885]
* Protista_83 [type: singlecopy] [num genes: 1741]
* Ribosomal_RNA_12S [type: Ribosomal_RNA_12S] [num genes: 12]
* Ribosomal_RNA_16S [type: Ribosomal_RNA_16S] [num genes: 40]
* Ribosomal_RNA_18S [type: Ribosomal_RNA_18S] [num genes: 12]
* Ribosomal_RNA_23S [type: Ribosomal_RNA_23S] [num genes: 26]
* Ribosomal_RNA_28S [type: Ribosomal_RNA_28S] [num genes: 12]
* Ribosomal_RNA_5S [type: Ribosomal_RNA_5S] [num genes: 67]
* Transfer_RNAs [type: Transfer_RNAs] [num genes: 547]
```

And generate a matrix for Bacteria_71:
```bash
anvi-script-gen-hmm-hits-matrix-across-genomes -e external-genomes.txt \
                                               -o Bacteria_71_matrix.txt \
                                               --hmm-source Bacteria_71
```

You can open and sort this matrix to see which are the most abundant bacterial SCGs.
Here is how we could check the top 10 most abundant SCGs from Bacteria_71:

```bash
# this command transpose the matrix
anvi-script-transpose-matrix Bacteria_71_matrix.txt \
                             -o Bacteria_71_matrix_transpose.txt

# the awk command sum the frequency of each SCG
awk 'NR>1{for(i=2;i<=NF;i++) t+=$i; print $1"\t"t; t=0}' Bacteria_71_matrix_transpose.txt | sort -nrk2 | head -n 10
```
Here is what I got:
```
Pept_tRNA_hydro	110
Chorismate_synt	109
ATP-synt	109
RBFA	107
PGK	107
Ham1p_like	106
Ribosomal_L21p	104
ATP-synt_A	104
Ribosomal_S15	103
Ribosomal_L27	102
```

Notice that the most abundant SCG is only found in 110/226 SAGs!

</div>

A next step is to compare the SAGs to each other. In other words - we will do some comparative genomics by comparing their gene content (pangenomics), compute genome similarity, and more!

## Pangenomics

In pangenomics, we examine the distribution of genes across genomes. Genes that are present in all genomes in the analysis are known as 'core genes', while genes that appear in only a subset of the genomes are called 'accessory genes'. Core genes often represent the essential functions that every population needs to have in order to survive, while accessory genes can encompass phenotypes specific to sub-group of organisms. You can read more about pangenomics [here](https://merenlab.org/2016/11/08/pangenomics-v2/).

For our pangenome analysis, we will focus on a set of SAGs from the Pelagibacterales and HIMB59 orders. Here is the list of all SAGs from these orders with a minimum completion estimate of 75%:
```
AG_910_D13
AG_910_I09
AG_910_K17
AG_910_M21
AG_910_M23
AG_910_A06
AG_910_O03
AG_910_J22
AG_910_F07
AG_910_G11
AG_910_C18
AG_910_C08
AG_910_B10
AG_910_K02
AG_910_K08
AG_910_K07
AG_910_G21
AG_910_C21
AG_910_C14
```


To add some context in this comparative genomic approach, we can add some reference genomes. In the directory `REFERENCE_GENOMES`, you will find two contigs dbs which were generated with the genomes of [HIMB083](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1096769) and [HIMB59](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=744985).

In the datapack, there is a new external genomes file which only Includes the list of selected SAGs and the two reference genomes. Let's copy it to our working directory:
```bash
cp ADDITIONAL_DATA/external-genomes-pangenomics.txt .
```

<div class="extra-info" markdown="1">

<span class="extra-info-header">Download and include NCBI genomes in anvi'o</span>
We often need to include genomes from public repositories to generate better comparisons. One common resource is NCBI.
[Here](https://merenlab.org/2019/03/14/ncbi-genome-download-magic/) is a comprehensive tutorial on how to obtain NCBI genomes with anvi'o.

In brief, you can install and use the program [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download) to get genomes from GenBank or RefSeq. It is very flexible: you can search for complete or draft genomes, specify genomes using their taxids, or ask for all genomes from a given genera, etc.

Then you can use the program {% include PROGRAM name="anvi-script-process-genbank-metadata" %} which takes the output of `ncbi-genome-download` to generate FASTA files, external gene calls files and functional annotation files. It also creates a `fasta.txt` file, which can directly be used by {% include PROGRAM name="anvi-run-workflow" %} to generate contigs dbs.

</div>

### Genome storage

A {% include ARTIFACT name="genomes-storage" text="genomes-storage database"%} is a special anvi'o database to store genome information. To create a genomes storage with the genomes from the list above, we can use the command {% include PROGRAM name="anvi-gen-genomes-storage" %} along with the new external genomes file:
```bash
# create a directory to store the pagenomics databases
mkdir -p PANGENOMICS

# create the genomes storage db
anvi-gen-genomes-storage -e external-genomes-pangenomics.txt \
                         -o PANGENOMICS/Pelagibacterales-HIMB59-GENOMES.db
```

Success? Good.

### Pangenome analysis

Once your genomes storage is ready, you can use {% include PROGRAM name="anvi-pan-genome" %} to run the actual analysis:
```bash
anvi-pan-genome -g PANGENOMICS/Pelagibacterales-HIMB59-GENOMES.db \
                -n Pelagibacterales-HIMB59 \
                -o PANGENOMICS\
                -T 4
```

In brief, when you use this program, anvi'o will use [DIAMOND](https://www.wsi.uni-tuebingen.de/lehrstuehle/algorithms-in-bioinformatics/software/diamond/) to compute the similarity between each amino acid sequence in every genome. Then, it uses the [MCL](http://micans.org/mcl) algorithm to identify clusters in the amino acid similarity search results. [Here you can find more information about the parameters of anvi-pan-genome](https://merenlab.org/2016/11/08/pangenomics-v2/).

### Display the pangenome

When the analysis is done, you can use {% include PROGRAM name="anvi-display-pan" %} to start an interactive web page containing the results:
```bash
anvi-display-pan -g PANGENOMICS/Pelagibacterales-HIMB59-GENOMES.db \
                 -p PANGENOMICS/Pelagibacterales-HIMB59-PAN.db
```

This is what you should see when you click the `Draw` button:

{% include IMAGE path="images/pangenome-raw.png" width=80 %}

It is pretty messy at the moment. To improve it, you can order the genomes based on the gene clusters' distribution by selecting the `gene cluster frequencies` tree from the Samples panel and Sample Order menu:

{% include IMAGE path="images/pangenome-order-by-freq.png" width=50 %}

And this is what you should see when you draw again:

{% include IMAGE path="images/pangenome-reorganized-freq.png" width=80 %}

It is a bit more organized. You can save the current state in the main panel. Every time you will start this interactive interface, it will load the `default` state. For now, let's close it and kill the server. We can add some meaningful information onto it, like the taxonomy of each of these SAGs! To do that, we can use the command {% include PROGRAM name="anvi-import-misc-data" %}. Right now, we have the file `taxonomy.txt` which has taxonomy information for all 228 of the SAGs, and anvi'o would be upset if we tried to import this table for a pangenome of only 21 genomes. But it is ok, we can use the flag `--just-do-it` and anvi'o will try its best.

```bash
anvi-import-misc-data -p PANGENOMICS/Pelagibacterales-HIMB59-PAN.db \
                      -t layers \
                      --just-do-it \
                      taxonomy.txt
```

You can rerun {% include PROGRAM name="anvi-display-pan" %} and take some time to explore the `Main` and `Layers` panels to modify some aesthetics.

```bash
anvi-display-pan -g PANGENOMICS/Pelagibacterales-HIMB59-GENOMES.db \
                 -p PANGENOMICS/Pelagibacterales-HIMB59-PAN.db
```

For instance, this is the same pangenome, just prettier:

{% include IMAGE path="images/pangenome-refined.png" width=80 %}

{:,notice}
You can export the figure as an SVG with the save button on the bottom of the setting panel and use your favorite SVG-eating software like [Inkscape](https://inkscape.org/).

There are three distinct groups of SAGs and reference genomes based on their gene content. Each group also has consistent taxonomy annotations, with the blue group (in the figure above) belonging to HIMB59, the red containing *Pelagibacter* genomes, and the green group has two GCA-002704185 sp902517865 genomes, which are members of the Pelagibacterales.


### Bin gene clusters and summarize the pangenomic analysis

A gene cluster can be part of the 'core' genome when it is shared by all the genomes, or part of the 'accessory' genome of a subset of genomes. Finally, if a gene cluster is only found in one genome, then we usually refers to it as a 'singleton' gene cluster. If you are interested in the set of genes that are unique to a subgroup, or to a single genome, you can use the interactive interface to 'bin' gene clusters and later summarize these bins to obtain the amino acid sequences of these genes and their functions.

To create a bin, you can simply click on gene clusters and they will be added to the current bin. A more efficient way to add gene clusters to a bin is to click on branches in the central dendrogram.

{:.notice}
You can use the `Bin` tab to select, rename, and create bins.

Let's select the accessory genomes of GCA-002704185 sp902517865 and HIMB59, and also select the singleton gene clusters of AG-910-G11 to further inspect these genes and therefore the functions that make these genomes unique compared to the rest.

After selecting branches on the dendrogram, and renaming the bins, you should have a figure that looks somewhat like this:

{% include IMAGE path="images/pangenome-binning.png" width=100 %}

{:.notice}
You can have a quick look at the content of a bin by clicking on the box with the number of gene clusters.

You can store the collection of bins by clicking on the "Store bin collection" button. Anvi'o lets you make as many collections of bins as you want. You can give the collection a meaningful name, but for now we can just save it as `default`.

A better way to inspect these bins is to use the program {% include PROGRAM name="anvi-summarize" %}. Kill the interactive interface and run the following command:

```bash
anvi-summarize -g PANGENOMICS/Pelagibacterales-HIMB59-GENOMES.db \
               -p PANGENOMICS/Pelagibacterales-HIMB59-PAN.db \
               -C default \
               -o SUMMARY
```
It will create the directory `SUMMARY`, which contains the file `Pelagibacterales-HIMB59_gene_clusters_summary.txt.gz`. To decompress this file, you can run:

```bash
gunzip SUMMARY/Pelagibacterales-HIMB59_gene_clusters_summary.txt.gz
```

Now you have a TAB-delimited file with **all** the genes occurring in the entire pangenome. There is also a column called `bin_name` to help you investigate the genes unique to GCA-002704185 sp902517865, HIMB59, and AG-910-G11.

### Compute average nucleotide identity

A common metric used to compare genomes is the average nucleotide identity (ANI). Anvi'o has a program, {% include PROGRAM name="anvi-compute-genome-similarity" %}, which uses different similarity metrics to compute ANI across your SAGs. It is a convenient program that can accept the external genomes file and optionally a {% include ARTIFACT name="pan-db" %} to display the result directly on top of the pangenome figure. Today, we will run [PyANI](https://pubs.rsc.org/en/content/articlelanding/2016/AY/C5AY02550H) using this program.

Here is how to run the command:
```bash
anvi-compute-genome-similarity -e external-genomes-pangenomics.txt \
                               -p PANGENOMICS/Pelagibacterales-HIMB59-PAN.db \
                               -o ANI \
                               --program pyANI \
                               -T 4
```

Once this is complete, we can visualize the pangenome with {% include PROGRAM name="anvi-display-pan" %}.

```bash
anvi-display-pan -g PANGENOMICS/Pelagibacterales-HIMB59-GENOMES.db \
                 -p PANGENOMICS/Pelagibacterales-HIMB59-PAN.db
```

At first, there is nothing new. But if you go in the `Layers` tab, you will see this:

{% include IMAGE path="images/pangenome-ani-select.png" width=60 %}

Select `ANI_percentage_identity` and click on `Redraw layer data`. Then you should see the ANI on the display.

{:.notice}
You probably want to change the `min` value from the interface to better appreciate the differences in the ANI heatmap across all the genomes.

{% include IMAGE path="images/pangeome-ani-heatmap.png" width=80 %}

You can now also organize the genomes based on the ANI value using the menu `Order by` in the `Layers` tab.
We can learn a few things with this analysis so far. First of all, AG-910-D13 is quite closely related to the reference genome HIMB59. Second of all, AG-910-G11 is still part of the of the same cluster as all the other *Pelagibacter* genomes, but it is clearly more different with a low ANI compare to the other genomes.



### Exploratory analyses

We know that single amplified genomes are not complete. If a gene is absent from a SAG but present in a closely-related genome, we cannot exclude the possibility that it was simply lost during the assembly and therefore cannot prove its absence. However, if a gene or a set of genes is absent from a reference genome but present in SAGs, then we can do science and make some claims.

If you zoom on this part of the pangenome, you will notice a gene cluster with no genes in HIMB083. You can right click on it and select `Inspect gene cluster`.

{% include IMAGE path="images/pangenome-zoom-inspect.png" width=100 %}
{% include IMAGE path="images/pangenome-inspect-page.png" width=80 %}

Here you can see the gene alignment with color-coded amino acids (more information on the coloring [here](https://merenlab.org/2018/02/13/color-coding-aa-alignments/)). You can also click on the gene call number on the left side to have a look at the functional annotation of these genes.

{% include IMAGE path="images/K13693.png" width=80 %}

A glycosyltransferase that is present in our SAGs but not in the reference genome - that's an interesting insight! But is it truly absent from the reference genome? What if we want to search for that function in all of the gene clusters? Well, it is possible to do with the search tab - select "Search functions" and type the COG accession for that function: COG1215. There are five gene clusters containing genes with the same COG annotation, and if you click on `Highlight splits on tree` you will be able to see where these gene clusters are. And none of them are present in the HIMB083 reference. If you search for the corresponding KEGG accession, you will find out that the gene cluster we inspected was the _only one_ with that annotation! This is a strong evidence that this glycosyltransferase is unique to our set of *Pelagibacter* SAGs.

With a quick search on the glucosyl-3-phosphoglycerate synthase gene, we can learn that it should be associated with a glucosyl-3-phosphoglycerate phosphatase for the production of glucosylglycerate ([Costa et al., 2006](https://doi.org/10.1128/JB.188.3.1022-1030.2006)) and that marine cyanobacteria synthesize glucosylglycerate as a secondary compatible solute in nitrogen-poor environments ([Klähn et al., 2009](https://doi.org/10.1111/j.1462-2920.2009.02045.x)).

Let's save our gene cluster of interest by clicking on it, which will add it to a bin.

To further investigate if our SAGs are potentially able to synthesize glucosylglycerate, we need to find the second gene which should be either just before or after the glucosyl-3-phosphoglycerate synthase gene in the genome. For that, we can order the gene clusters by the synteny of one genome. To do this, you can use the `Items order` menu in the main tab and select `Forced synteny <> AG 910 K02` and redraw the figure.

If you zoom on the bin we just created, you can see that the neighbor genes are all absent in HIMB083:
{% include IMAGE path="images/pangenome-synteny.png" width=80 %}

Now, inspect these gene clusters yourself and try to find an hypothetical 'glucosyl-3-phosphoglycerate phosphatase'. Go for it, I'm waiting.

By now you have probably found a gene cluster, just next to our gene of interest, and found out that it has a 'mannosyl-3-phosphoglycerate phosphatase' annotation. These two genes are found side by side, and are only present in the SAGs, and not in the reference genome! Now it is up to you, as a scientist, to step in and decide what you want to do with this result.


## Phylogenomics

Yet another 'omics! Phylogenomics is very similar to phylogenetics in the sense that the goal is to generate a tree and study genomes from an evolutionary point of view. While phylogenetics compares one gene across a set of genomes, phylogenomics includes multiple genes for a larger alignment and possibly a more accurate description of the evolutionary relationships.

Let's make a phylogenomic tree of the 21 genomes we used in the pangenomic analysis.

We can use {% include PROGRAM name="anvi-get-sequences-for-hmm-hits" %} like we did for the phylogenetic analysis, this time to  get multiple genes instead of just one. But in closely-related genomes, some single-copy core genes will be 100% similar across all genomes. If we include them in the alignment, they will not contribute to to the phylogenomic signal. In addition, we don't need to use the small set of bacterial single-copy core genes when we can use the 'core' genome specific to our small set of Pelagibacterales.

In the pangenomic analysis, anvi'o computes two metrics for each gene cluster: functional homogeneity and geometric homogeneity. Because not all gene clusters are created equal, we use these metrics to distinguish very homogenous gene clusters with very high amino acid sequence similarity from less homogenous one with divergent amino acid sequences. The latter category would be interesting to use for a phylogenomic analysis.

{:.notice}
**Geometric homogeneity** indicates the degree of geometric configuration between the genes of a gene cluster. Basically, it takes into account the gap/residue pattern. **Functional homogeneity** uses the aligned residues and attempts to quantify their differences considering the biochemical properties of amino acids, which would affect gene structure/function.

### Filter gene clusters

To select a suitable set of genes, we can use the pangenome interface. If you closed it, you can run:

```bash
anvi-display-pan -g PANGENOMICS/Pelagibacterales-HIMB59-GENOMES.db \
                 -p PANGENOMICS/Pelagibacterales-HIMB59-PAN.db
```

You can filter for specific gene clusters using the `Search` tab and the `Search gene clusters using filters` section:

{% include IMAGE path="images/phylogenomics-filter.png" width=50 %}

Let's try different filters to find the best set of genes. We can start by choosing only the single-copy core gene clusters. For that you want to filter for gene clusters occurring in all genomes (min 21), and with a maximum of 1 gene per genome.

{% include IMAGE path="images/phylogenomics-scg.png" width=50 %}

If you search for them and highlight the results on the tree, you will see that they correspond to the 'SCG cluster' layer, which makes total sense. That's 13 gene clusters, which is not a lot. Unfortunately the core genome is pretty small in this pangenome, which is likely due to the fragmented nature of the SAGs. So we can be a little bit more flexible and search for gene clusters occurring in at least 19 out of 21 genomes by changing the appropriate filter. We can also filter for less homogenous gene clusters to keep only the most divergent ones (max functional homogeneity of 0.8):

{% include IMAGE path="images/phylogenomics-correct-filters.png" width=50 %}

You can now add these gene clusters to a bin (make sure you're using an empty bin or create a new one in the `Bin` tab) by clicking on the button `Append splits to selected bin`.

### Compute and display the phylogenomic tree

Now to make the alignment _and_ the phylogenomic tree, you can simply navigate to the `Layers` tab where you will see a button `Generate a new phylogenetic tree`. When you click on it you can name your tree, select the bin you have just made, choose an aligner, and finally check the box to store the tree permanently.

{% include IMAGE path="images/phylogenomics-make-tree.png" width=50 %}

Anvi'o will automatically display the tree and reorganize the genomes accordingly:

{% include IMAGE path="images/phylogenomics-display.png" width=80 %}


Unsurprisingly, we can still see three clusters corresponding to the three genera, and these clusters fit well with the ANI heatmap.
Notice how AG-910-G11 is now more outside of the *Pelagibacter* cluster than we previously observed?


<div class="extra-info" markdown="1">

<span class="extra-info-header">Make a phylogenomic tree outside of the interactive interface</span>

You don't need to use the interactive interface to extract, align, and compute a phylogenomic tree. Instead, you can use the program {% include PROGRAM name="anvi-get-sequences-for-gene-clusters" %}. It accepts a genome storage and a pan db, and you can use the same filters as the ones in the `Search` tab of the interactive interface. For instance, this command would give you the same set of genes as a multiple sequence alignment:

```bash
anvi-get-sequences-for-gene-clusters -g PANGENOMICS/Pelagibacterales-HIMB59-GENOMES.db \
                                     -p PANGENOMICS/Pelagibacterales-HIMB59-PAN.db \
                                     -o concatenated-proteins-Pelagibacterales-HIMB59.fa \
                                     --concatenate-gene-clusters \
                                     --min-num-genomes-gene-cluster-occurs 19 \
                                     --max-num-genes-from-each-genome 1 \
                                     --max-functional-homogeneity-index 0.8
```

You can then compute a tree with {% include PROGRAM name="anvi-gen-phylogenomic-tree" %}:
```bash
anvi-gen-phylogenomic-tree -f concatenated-proteins-Pelagibacterales-HIMB59.fa \
                           -o phylogenomic-tree-Pelagibacterales-HIMB59.txt
```
And use {% include PROGRAM name="anvi-interactive" %} to display the tree:

```bash
anvi-interactive -p phylogenomic-profile-Pelagibacterales-HIMB59.db \
                 -t phylogenomic-tree-Pelagibacterales-HIMB59.txt \
                 --title "Phylogenomic_tree_Pelagibacterales-HIMB59" \
                 -A taxonomy.txt \
                 --manual
```

You can always import that tree into an existing pangenome with {% include PROGRAM name="anvi-import-misc-data" %}, but you will need to reformat the tree to fit the input format - click [here](https://merenlab.org/2017/12/11/additional-data-tables/#layer-orders-additional-data-table) to learn more.

 </div>


## Read recruitment

It is time to dive into the ecology of our SAGs and investigate where they occur! The set of SAGs we are using is coming from the Bermuda-Atlantic Time-series Study (BATS) Station in the Sargasso Sea and we can use the publicly available metagenomes from this region to recruit reads onto the SAGs. For this tutorial, we will focus on 20 metagenomes from samples collected approximately every month between 2003 and 2004. We will also use the metagenomes from the same sample used to generate the AG-910 SAGs that we have been working on (sample collected in July 2009). [Here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA385855) is a link to the NCBI BioProject for the relevant samples.

The list of SRA accessions that we will use is:
```
SRR5720233 SRR5720238 SRR5720327 SRR5720283 SRR5720235 SRR5720286 SRR5720332
SRR5720276 SRR5720262 SRR5720338 SRR5720322 SRR5720337 SRR5720256 SRR5720257
SRR5720260 SRR5720321 SRR5720251 SRR5720307 SRR5720278 SRR5720342 SRR6507279
```

{:.notice}
The read recruitment analysis can be quite computationally intensive and and can generate relatively large files. The datapack only includes the final files generated with anvi'o and we will focus on the visualization of the results in this tutorial.

<div class="extra-info" markdown="1">

<span class="extra-info-header">Generate the profile database in anvi'o</span>
This is a short summary of the workflow used to generate the {% include ARTIFACT name="profile-db" %} that stores the read recruitment results and can be used by {% include PROGRAM name="anvi-interactive" %} to visualize these results.

To download the metagenomes, you can use `fastq-dump` from the [sra toolkit](https://hpc.nih.gov/apps/sratoolkit.html).

Briefly, the analysis consists of using bowtie to map short reads onto our contigs. Then we use samtools to convert SAM files into BAM files and remove any unmapped reads. Then we create a {% include ARTIFACT name="single-profile-db" %} for every sample. _Single_ profile databases contain the key information about the read mapping from a _single_ sample to your contigs: coverage, abundance, single nucleotide variants, single codon variants, single amino acid variants, and insertions/deletions (indels) per nucleotide position. Finally, we merge all the {% include ARTIFACT name="single-profile-db" text="single profile databases" %} into a _merged_ {% include ARTIFACT name="profile-db" %} (which contains the combined information from all samples that were mapped to the same contigs database).

You can run all these steps manually, as described in this [comprehensive step-by-step exercise on read recruitment](https://merenlab.org/momics/exercises/read-recruitment/). But the {% include ARTIFACT name="profile-db" %} in our datapack was automatically generated using the 'metagenomics' snakemake workflow integrated in the command {% include PROGRAM name="anvi-run-workflow" %}.

This workflow requires a {% include ARTIFACT name="fasta-txt" %} file with the list of SAGs onto which you want to do the read recruitment; a {% include ARTIFACT name="samples-txt" %} file including the path to each metagenome's short reads; and a {% include ARTIFACT name="workflow-config" %} file that describes the steps and parameters considered by {% include PROGRAM name="anvi-run-workflow" %}.

You can create a default {% include ARTIFACT name="workflow-config" %} file using the following command:
```bash
anvi-run-workflow --workflow metagenomics \
                  --get-default-config config.json
```

You can modify the `config.json` and choose all the steps that you want to be run on your samples during the workflow, but also the number of threads for each command, output directory names, and more. To visualize the workflow steps, you can generate a workflow graph by using the flag `--save-workflow-graph` when running {% include PROGRAM name="anvi-run-workflow" %}.

Here is the graph for the workflow used to generate the {% include ARTIFACT name="profile-db" %} in the datapack:

{% include IMAGE path="images/metagenomics-workflow.png" width=90 %}

</div>

### The interactive interface

First, let's change to the right directory:

```bash
cd ../READ_RECRUITMENT
```

To start the interactive interface, you can use this command:
```bash
anvi-interactive -c CONTIGS/AG_910_O03-contigs.db -p PROFILE/PROFILE.db
```

If you click on the `Draw` button, you should see this:

{% include IMAGE path="images/read-recruitment-raw.png" width=80 %}

The central dendrogram represents contigs, or more precisely, 'splits' - as anvi'o splits large contigs into smaller (~20,000bp) sequences for a better visualization experience. Every layer represents a metagenome and the bar plot represents the mean coverage of a 'split' in a given metagenome.

{:.notice}
By default the central dendrogram is computed using differential coverage and the sequence composition of the contigs (kmer composition, default k=4). You can also use only differential coverage, or only sequence composition with the `View` menu in the `Main` tab. Don't forget to re-draw the figure to see the new contig organization.

You can polish the figure a little bit, especially by highlighting the percent of reads mapping to appreciate the seasonality pattern of *Pelagibacter*.

{% include IMAGE path="images/read-recruitment-polished.png" width=80 %}

You probably immediately notice that three splits don't have any coverage in all of the metagenomes! To check if it is contamination, you can select or 'bin' these splits with no coverage and open the `Bin` tab. You will see that these three splits only add up to ~42 kbp, nothing close to a full genome. You can also check for the presence of single-copy core genes present in these splits compared to the rest of the genome; you will see that there are none.

Another interesting observation is that one of the splits is actually part of a large contig containing an rRNA operon, so unless there is a misassembly, this split is definitely part of the genome.


### The inspection page
To further investigate the content of splits, you can right click on one split and select `Inspect split`. It will open another tab on your browser with a more detailed view of the coverage. Using the `Settings` panel, you can choose to display SNVs and indels, overlay GC content, and color genes by their COG category or KEGG class.

{% include IMAGE path="images/read-recruitment-inspect.png" width=80 %}

If you select a gene, you will see its functional annotation (if any) and you will be able to get its DNA and amino acid sequences. You can also directly start various BLAST jobs against NCBI nucleotide collection (nr) or the RefSeq genome database, using either `blastn` for the nucleotide sequence or `blastx` for the translated amino acid sequence.

Try inspecting the two splits with no coverage and have a look at the gene annotations (or lack thereof) to better understand why it is not detected in the environment. You will see that they contain genes associated with the cell wall and membrane!

### Summarizing the read recruitment

Similar to the pangenomic analysis, we can summarize the results of the read recruitment using {% include PROGRAM name="anvi-summarize" %}.
To summarize the whole SAG, you need to tell anvi'o that you want everything in that contigs db, not just a subset of contigs. To do that, we can create a 'bin' which contains all of the contigs. You could do it in the interactive interface and save it in a collection of your choice, or you can use the program {% include PROGRAM name="anvi-script-add-default-collection" %}:

```bash
anvi-script-add-default-collection -c CONTIGS/AG_910_O03-contigs.db -p PROFILE/PROFILE.db
```

And summarize the collection called 'DEFAULT' (or whatever name you gave it):
```bash
anvi-summarize -c CONTIGS/AG_910_O03-contigs.db \
               -p PROFILE/PROFILE.db \
               -C DEFAULT \
               -o SUMMARY
```

Here are all the files generated for this SAG ('EVERYTHING' is the default bin name created by {% include PROGRAM name="anvi-script-add-default-collection" %}; it could be different if you named your bin in the interactive interface):

```bash
$ ls SUMMARY/bin_by_bin/EVERYTHING/
EVERYTHING-Archaea_76-hmm-sequences.txt        EVERYTHING-Ribosomal_RNA_18S-hmm-hits.txt      EVERYTHING-Transfer_RNAs-hmm-sequences.txt     EVERYTHING-percent_completion.txt
EVERYTHING-Bacteria_71-hmm-sequences.txt       EVERYTHING-Ribosomal_RNA_18S-hmm-sequences.txt EVERYTHING-abundance.txt                       EVERYTHING-percent_redundancy.txt
EVERYTHING-GC_content.txt                      EVERYTHING-Ribosomal_RNA_23S-hmm-hits.txt      EVERYTHING-contigs.fa                          EVERYTHING-scg_domain.txt
EVERYTHING-N50.txt                             EVERYTHING-Ribosomal_RNA_23S-hmm-sequences.txt EVERYTHING-detection.txt                       EVERYTHING-scg_domain_confidence.txt
EVERYTHING-Protista_83-hmm-sequences.txt       EVERYTHING-Ribosomal_RNA_28S-hmm-hits.txt      EVERYTHING-gene_calls.txt                      EVERYTHING-scg_taxonomy_details.txt
EVERYTHING-Ribosomal_RNA_12S-hmm-hits.txt      EVERYTHING-Ribosomal_RNA_28S-hmm-sequences.txt EVERYTHING-mean_coverage.txt                   EVERYTHING-std_coverage.txt
EVERYTHING-Ribosomal_RNA_12S-hmm-sequences.txt EVERYTHING-Ribosomal_RNA_5S-hmm-hits.txt       EVERYTHING-mean_coverage_Q2Q3.txt              EVERYTHING-total_length.txt
EVERYTHING-Ribosomal_RNA_16S-hmm-hits.txt      EVERYTHING-Ribosomal_RNA_5S-hmm-sequences.txt  EVERYTHING-num_contigs.txt                     EVERYTHING-variability.txt
EVERYTHING-Ribosomal_RNA_16S-hmm-sequences.txt EVERYTHING-Transfer_RNAs-hmm-hits.txt          EVERYTHING-original_split_names.txt
```

If you were interested in only a part of the genome, like the contigs with no coverage, you can always bin these contigs to create a different collection and then summarize that collection. This would give you the contigs' sequences, gene calls, functional annotations, and more.

{:.notice}
In this example, we have only used one SAG for read recruitment. One of the drawbacks with this is the large amount of non-specific mapping (from different populations in the environment) to conserved genes/operons. You could recruit reads like this to every single SAG, one by one. But a better strategy to reduce non-specific mapping would be to concatenate all SAGs (or a subset of them) into a single FASTA file and map short-reads onto all of them at the same time. In that case, the read mapping would be competitive and you would be able to appreciate the variation of coverage between SAGs in the same sample. In order to reconcile which contigs in the contigs db belong to which SAG, you can import a collection listing all of the contigs and their associated SAG with {% include PROGRAM name="anvi-import-collection" %}.

## Other anvi'o resources

Thanks for following this tutorial which shows you how to integrate single-cell genomics in anvi'o. There are many possible analysis that you can perform and we have just scratched the surface here. In the end, it is up to you to choose what you want to do with your data and what question you want to address.

There are many more tutorials and examples analysis with anvi'o and here is a list of some of them which could be relevant for your science:
- [The infant gut tutorial](https://merenlab.org/tutorials/infant-gut/) is an all in one tutorial using a small set of metagenomes. I covers genome-resolved metagenomics, binning, phylogenomics, pangenomics, metabolism prediction, population genetics and protein structure variation.
- [Studying microbial population genetics with anvi'o](https://merenlab.org/2015/07/20/analyzing-variability/), an insight into SNVs, SCVs (codon) and SAAVs (amino acid) variation analysis.
- [A simple read recruitment exercise](https://merenlab.org/tutorials/read-recruitment/).
- [An exercise on metabolic reconstruction](https://merenlab.org/tutorials/fmt-mag-metabolism/).
- [A tutorial on the anvi'o interactive interface](https://merenlab.org/tutorials/interactive-interface/).

And more at [https://anvio.org/#learn](https://anvio.org/#learn).

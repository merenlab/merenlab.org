---
layout: page
title: Metabolic enrichment of high-fitness MAGs from a longitudinal FMT study
modified: 2021-10-22
authors: [iva]
excerpt: "A mini-tutorial on how to run metabolism estimation and enrichment in anvi'o"
categories: [anvio]
comments: true
redirect_from: /fmt-metabolism/
---

This is a **mini-tutorial for the metabolism suite of programs** in anvi'o. First, we will learn how to estimate metabolism for a single bacterial isolate, starting from its genome sequence and ending with a file of metabolic module completeness scores. Afterwards, we'll be applying this to a larger, real-world dataset of metagenome-assembled genomes from [our recent FMT study](ttps://doi.org/10.1101/2021.03.02.433653), to learn how to estimate metabolism in a more high-throughput manner as well as how to compute enrichment scores for metabolic modules.

{:.notice}
This tutorial is tailored for anvi'o `v7.1` or later. You can learn the version of your installation by running `anvi-interactive -v` in your terminal. If you are using the development branch of anvi'o, look for the Show/Hide boxes.

## Running the metabolism suite of programs - a single-genome example

For this example, we'll be using a bacterial isolate genome from the NCBI - the representative genome for _Akkermansia muciniphila_. I picked this species because it is a pretty cool gut microbe that can degrade mucin ([Derrien 2004](https://doi.org/10.1099/ijs.0.02873-0)). But if you have your own data or an interest in a different species, feel free to use that instead.

Here is how you can download and unpack the _A. muciniphila_ genome:

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/731/575/GCF_009731575.1_ASM973157v1/GCF_009731575.1_ASM973157v1_genomic.fna.gz && \
gunzip GCF_009731575.1_ASM973157v1_genomic.fna.gz
```

The first step in any anvi'o analysis is to get your data into the form that anvi'o likes to play with, which means we need to make a {% include ARTIFACT name="contigs-db" text="contigs database" %}. To do this, we first need to re-format the FASTA file to make sure all of the contig deflines contain nothing more than alphanumeric characters, underscores, and dashes. We can do this with {% include PROGRAM name="anvi-script-reformat-fasta" %}. Then we can pass the reformatted FASTA to {% include PROGRAM name="anvi-gen-contigs-database" %} and let it work its magic:

```bash
# fix the deflines
anvi-script-reformat-fasta GCF_009731575.1_ASM973157v1_genomic.fna \
                          --simplify-names \
                          --output-file A_muciniphila.fa

# convert FASTA into a contigs db
anvi-gen-contigs-database --contigs-fasta A_muciniphila.fa \
                          --project-name A_muciniphila \
                          --output-db-path A_muciniphila-CONTIGS.db \
                          --num-threads 2
```

{:.warning}
The command above sets the `--num-threads` parameter to 2 so that `anvi-gen-contigs-database` works a bit faster, but please take care that your computer has enough CPUs to accommodate this number of threads. This applies anytime you see the `--num-threads` or `-T` parameter in a command.

Now you have a contigs database for this genome, so we can start to work with it to estimate metabolism. There are 3 required steps in the metabolism estimation process (though the first step is only necessary the very first time you are doing this, so in general there are just 2 steps). Those steps are:

1. Getting the required [KEGG](https://www.genome.jp/kegg/) data set up on your computer (which only needs to be done once) with {% include PROGRAM name="anvi-setup-kegg-kofams" %}
2. Adding functional annotations from the [KOfam database](https://doi.org/10.1093/bioinformatics/btz859) to the gene calls in this genome, using {% include PROGRAM name="anvi-run-kegg-kofams" %}
3. Matching those functional annotations to [definitions of metabolic pathways](https://www.genome.jp/kegg/module.html) (aka modules) and computing module completeness scores, which is done by the program {% include PROGRAM name="anvi-estimate-metabolism" %}

You can learn more about how each of these programs work (and different options for running them) by clicking on any of the highlighted links above to go to their respective help pages. Below we will go through a simple example of each step using the _A. muciniphila_ genome that we just downloaded.

### Setting up KEGG data on your computer

First, if you've never worked with KEGG data through anvi'o before (which is likely, considering what you are reading right now), you will need to get this data onto your computer so that the downstream programs can use it. This is as simple as running the following:

```bash
anvi-setup-kegg-kofams
```

That's it. What this does is download 1) KOfam profile hidden Markov models (pHMMs) and 2) KEGG Module definition files onto your computer. It then organizes the pHMMs into one big file that is ready for running [`hmmsearch`](https://doi.org/10.1371/journal.pcbi.1002195) for annotations, and it parses the Module definitions into a {% include ARTIFACT name="modules-db" text="modules database" %}. Once you have these things on your computer, you won't need to run this program again (until you want to update this database with a new version).

By default, this data goes into the anvi'o directory on your computer. If you don't have permission to modify this folder, you will need to pick a different location (that you _do_ have permission to modify) for the data and specify that folder using the `--kegg-data-dir` parameter. If this is your case, please note that the two subsequent steps will also require you to specify that folder location with `--kegg-data-dir`.

### Annotating the genome with KOfam hits

The contigs database includes predicted gene calls (open reading frames, or ORFs), but we don't know what functions these genes encode. The next step is annotating these genes with hits to the KOfam database of functional orthologs. If a gene is similar enough to a protein family in this database, it will be annotated with the family's KEGG Ortholog (KO) number. This number is what allows us to match genes to the metabolic pathways that they belong to.

This step is the most time-consuming in the workflow. If you have enough resources on your computer, you can give it additional threads to speed up the process (below we use just 4 threads, which should work for most laptop models these days).

```bash
anvi-run-kegg-kofams -c A_muciniphila-CONTIGS.db \
                     -T 4
```

When this finishes running (it took about 5 minutes on my computer), the output on your terminal should tell you how many KOfam hits were added to the contigs database. For the _A. muciniphila_ genome, this number should be 1,263:

```
Gene functions ...............................: 1,263 function calls from 1 source (KOfam) for 1,215 unique gene calls have been added to the contigs database.
Gene functions ...............................: 334 function calls from 1 source (KEGG_Module) for 322 unique gene calls have been added to the contigs database.
Gene functions ...............................: 334 function calls from 1 source (KEGG_Class) for 322 unique gene calls have been added to the contigs database.
```

Of these 1,263 annotations, 334 are KOs that belong to one or more metabolic pathways in the KEGG Module database. These are the KOs that will be used to estimate module completeness in the next step.

### Estimating metabolism

Our final step is to calculate the completeness of each pathway in the KEGG Module database. {% include PROGRAM name="anvi-estimate-metabolism" %} will go through the definition of each module and compute the fraction of KOs in this definition that are present in the genome.

```bash
anvi-estimate-metabolism -c A_muciniphila-CONTIGS.db \
                         -O A_muciniphila
```

When you run this command, you should see in your terminal that 58 modules were found to be complete in this genome. Here, 'complete' means that the genome contains at least 75% of the KOs necessary to complete the metabolic pathway. This threshold is mutable - you can change it using the `--module-completion-threshold` parameter. However, while this threshold is useful as a basic filter to narrow down which modules are worth considering, it should not be used to conclusively decide which modules are actually complete or not in this genome, as discussed [here](https://merenlab.org/tutorials/infant-gut/#estimating-metabolism-in-the-enterococcus-genomes).

The program will produce an output file called `A_muciniphila_modules.txt` which describes the completeness of each module. More information about this output (and other available output modes) can be found {% include ARTIFACT name="kegg-metabolism" text="at this link" %}.

When you take a look at the output file, you will see that many of the modules marked as 'complete' (in the `module_is_complete` column) are biosynthesis pathways. The table below shows a few of these.

unique_id | genome_name | kegg_module | module_name | module_class | module_category | module_subcategory | module_definition | module_completeness | module_is_complete | kofam_hits_in_module | gene_caller_ids_in_module | warnings
:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
9 | A_muciniphila | M00005 | PRPP biosynthesis, ribose 5P => PRPP | Pathway modules | Carbohydrate metabolism | Central carbohydrate metabolism | "K00948" | 1.0 | True | K00948 | 1580 | None
16 | A_muciniphila | M00854 | Glycogen biosynthesis, glucose-1P => glycogen/starch | Pathway modules | Carbohydrate metabolism | Other carbohydrate metabolism | "(K00963 (K00693+K00750,K16150,K16153,K13679,K20812)),(K00975 (K00703,K13679,K20812)) (K00700,K16149)" | 1.0 | True | K00700,K00963,K16150 | 1933,1923,2048 | None
51 | A_muciniphila | M00083 | Fatty acid biosynthesis, elongation | Pathway modules | Lipid metabolism | Fatty acid metabolism | "K00665,(K00667 K00668),K11533,((K00647,K09458) K00059 (K02372,K01716,K16363) (K00208,K02371,K10780,K00209))" | 1.0 | True | K00059,K00059,K00059,K00208,K02372,K09458,K16363 | 1055,238,918,1821,2092,1438,2092 | None
57 | A_muciniphila | M00093 | Phosphatidylethanolamine (PE) biosynthesis, PA => PS => PE | Pathway modules | Lipid metabolism | Lipid metabolism | "K00981 (K00998,K17103) K01613" | 1.0 | True | K00981,K01613,K17103 | 1811,1830,1448 | None
58 | A_muciniphila | M00048 | Inosine monophosphate biosynthesis, PRPP + glutamine => IMP | Pathway modules | Nucleotide metabolism | Purine metabolism | "K00764 (K01945,K11787,K11788,K13713) (K00601,K11175,K08289,K11787,K01492) (K01952,(K23269+K23264+K23265),(K23270+K23265)) (K01933,K11787,K11788 (K01587,K11808,K01589 K01588)) (K01923,K01587,K13713) K01756 (K00602,(K01492,K06863 K11176))" | 1.0 | True | K00602,K00764,K01588,K01756,K01923,K01933,K01945,K11175,K23265,K23270 | 1274,506,123,2115,582,505,2041,1389,585,584 | None
59 | A_muciniphila | M00049 | Adenine ribonucleotide biosynthesis, IMP => ADP,ATP | Pathway modules | Nucleotide metabolism | Purine metabolism | "K01939 K01756 (K00939,K18532,K18533,K00944) (K00940,K00873,K12406)" | 1.0 | True | K00873,K00939,K00940,K01756,K01939 | 479,761,1759,2115,2303 | None
60 | A_muciniphila | M00050 | Guanine ribonucleotide biosynthesis IMP => GDP,GTP | Pathway modules | Nucleotide metabolism | Purine metabolism | "K00088 K01951 K00942 (K00940,K18533,K00873,K12406)" | 1.0 | True | K00088,K00873,K00940,K00942,K01951,K01951 | 961,479,1759,1598,960,2129 | None
61 | A_muciniphila | M00051 | Uridine monophosphate biosynthesis, glutamine (+ PRPP) => UMP | Pathway modules | Nucleotide metabolism | Pyrimidine metabolism | "(K11540,(K11541 K01465),((K01954,K01955+K01956) (K00609+K00610,K00608) K01465)) (K00226,K00254,K17828) (K13421,K00762 K01591)" | 0.9166666666666666 | True | K00254,K00609,K00762,K01465,K01591,K01955,K01955,K01956 | 421,121,2288,120,1291,1359,671,1358 | None
62 | A_muciniphila | M00052 | Pyrimidine ribonucleotide biosynthesis, UMP => UDP/UTP,CDP/CTP | Pathway modules | Nucleotide metabolism | Pyrimidine metabolism | "(K13800,K13809,K09903) (K00940,K18533) K01937" | 1.0 | True | K00940,K01937,K09903 | 1759,179,1430 | None
67 | A_muciniphila | M00021 | Cysteine biosynthesis, serine => cysteine | Pathway modules | Amino acid metabolism | Cysteine and methionine metabolism | "(K00640,K23304) (K01738,K13034,K17069)" | 1.0 | True | K00640,K01738,K01738 | 577,1408,2187 | None

(This output was obtained by running the following code: )

```bash
head -n 1 A_muciniphila_modules.txt; awk -F'\t' '$9 > 0.75' A_muciniphila_modules.txt | grep -i 'biosynthesis' | head -n 10
```

Clearly, this is a very talented microbe. It can make a lot of things.

Unfortunately, KEGG does not have a module for mucin degradation, so we won't see evidence of that metabolic capability here. This happens a lot with metabolisms that go beyond the basic ones essential for life, because KEGG is a manually curated resource that hasn't yet gotten to include a lot of the more niche metabolisms out there.

There is a way to get around this limitation, and that is to look at individual KOfam hits for KOs which do not belong to a particular metabolic module, but may be representative of a metabolism of interest. In the case of mucin degradation, the enzymes that break up mucin (by destroying the gylcosidic bonds between the mucin molecules) are called Glycoside hydrolases (GHs). The GH family includes many different types of proteins, including sialidases ([Tailford 2015](https://www.frontiersin.org/articles/10.3389/fgene.2015.00081/full)). There is a KO family for sialidases - [K01186](https://www.genome.jp/dbget-bin/www_bget?ko+K01186) - which means that we can look for genes annotated with this KO as evidence of this microbe's mucin degrading capabilities.

This requires us to obtain a different output type from `anvi-estimate-metabolism`: ["kofam_hits"](https://merenlab.org/software/anvio/help/main/artifacts/kegg-metabolism/#kofam-hits-mode) mode output files have an entry for each gene annotated with a KO in the contigs database, regardless of whether that KO belongs to a metabolic module or not.

This is how you get that file:

```bash
anvi-estimate-metabolism -c A_muciniphila-CONTIGS.db \
                         -O A_muciniphila \
                         --kegg-output-modes kofam_hits
```

<details markdown="1"><summary>Show/hide anvi'o development version </summary>
```bash
anvi-estimate-metabolism -c A_muciniphila-CONTIGS.db \
                         -O A_muciniphila \
                         --output-modes hits
```
</details>

And we can search for the sialidase KO by running the following code:

```bash
head -n 1 A_muciniphila_kofam_hits.txt; \
grep K01186 A_muciniphila_kofam_hits.txt
```

<details markdown="1"><summary>Show/hide anvi'o development version </summary>
```bash
head -n 1 A_muciniphila_hits.txt; \
grep K01186 A_muciniphila_hits.txt
```
</details>

You should see the following output:

unique_id | genome_name | ko | gene_caller_id | contig | modules_with_ko | ko_definition
:---|:---|:---|:---|:---|:---|:---|
1029 | A_muciniphila | K01186 | 678 | c_000000000001 | None | sialidase-1 [EC:3.2.1.18]
1030 | A_muciniphila | K01186 | 2015 | c_000000000001 | None | sialidase-1 [EC:3.2.1.18]

We already know that this organism can degrade mucin, so it is not surprising that there are two copies of the sialidase enzyme encoded in this genome.

## Metabolism estimation and enrichment on a real-world dataset

Now that we know how to work with this suite of programs, let's apply them to a larger set of genomes.

The data we'll be using for this is a real dataset from one of our recent studies, ["Metabolic competency drives microbial colonization and resilience in health and disease”](https://doi.org/10.1101/2021.03.02.433653) by Watson et al. In fact, this post is doing double-duty as a reproducible workflow for one of the analyses in that study. :) The rest of the reproducible workflow can be found [here](https://merenlab.org/data/fmt-gut-colonization/) for anyone who is interested in how we did the other analyses discussed in the paper.

{:.notice}
This part of the tutorial **only** works with anvio-7.1 and not with the current dev version of anvi'o. The reason is KEGG database version compatibility. When you run `anvi-setup-kegg-kofams` with v7.1, anvi'o will download the latest snapshot of KEGG at the time of v7.1 release (v2020-12-23). But if you are using the development version of anvi'o, every user *could* have different version of the KEGG database depending on when you ran `anvi-setup-kegg-kofams`. This makes it more difficult to share annotated contigs.db between collaborators. More information [here](https://anvio.org/help/main/programs/anvi-setup-kegg-kofams/#how-do-i-share-this-data).

### A dataset of high- and low-fitness MAGs

You can download the [datapack](https://figshare.com/ndownloader/files/31120057) for this tutorial by running the following code:

```bash
wget https://figshare.com/ndownloader/files/31120057 -O FMT_MAGS_FOR_METABOLIC_ENRICHMENT.tar.gz
tar -xvf FMT_MAGS_FOR_METABOLIC_ENRICHMENT.tar.gz && cd FMT_MAGS_FOR_METABOLIC_ENRICHMENT/
```

This dataset includes anvi'o {% include ARTIFACT name="contigs-db" text="contigs databases" %} for 40 MAGs of gut microbes, labeled as either "high-fitness" or "low-fitness" according to their [colonization ability](https://merenlab.org/data/fmt-gut-colonization/#defining-colonization-success-and-failure) and prevalence in healthy gut metagenomes (there are 20 MAGs in each group). You can learn the full details of how and why we got them by reading the [study](https://doi.org/10.1101/2021.03.02.433653), but for the purposes of this mini-tutorial, here is what you need to know about these MAGs:

- they were binned from a co-assembly of **longitudinally-sampled gut metagenomes** taken from a **healthy adult who donated stool for fecal microbiota transplantation (FMT)**
- the **high-fitness MAGs** represent microbial **populations that were able to colonize all FMT recipients** who received stool from this donor. They were detected (with sufficient abundance) in recipient gut metagenomes at least 7 days (and up to 1 year) post-FMT
- the **low-fitness MAGs** represent **populations that were NOT able to colonize** FMT recipients who received stool from this donor. They were not detected in recipient gut metagenomes post-FMT
- the **high-fitness MAGs** were the 20 MAGs with the **highest prevalence in gut metagenomes from healthy Canadian adults**, while the **low-fitness MAGs were less prevalent**

The "high-fitness" MAGs were labeled that way because we hypothesized that something about these populations increased their fitness such that they were able to survive the stress of being transplanted into new gut environments and become long-term colonizers (in comparison to the "low-fitness" populations that were unable to survive for long in the recipients). In our study, we sought to learn what distinguishes these two groups from each other - what enables one group to survive while the other does not? What do the "high-fitness" populations have that the "low-fitness" ones don't (or vice versa)?

One way to answer this question is to look at the metabolic potential, or genomically-encoded metabolic capabilities, of these MAGs. Luckily, we have just learned how to do that.

### Estimating metabolism for these MAGs

You can run {% include PROGRAM name="anvi-estimate-metabolism" %} on all 40 MAGs in this dataset at once by utilizing the {% include ARTIFACT name="external-genomes" text="external genomes file" %} provided in the datapack, as shown below:

```bash
anvi-estimate-metabolism -e external-genomes.txt \
                         -O FMT_MAG_metabolism
```

This will give you one [modules mode](https://merenlab.org/software/anvio/help/main/artifacts/kegg-metabolism/#modules-mode) output file called `FMT_MAG_metabolism_modules.txt`, which describes the completeness of each KEGG Module in each MAG.

You could look through this file manually to see what metabolisms are encoded in these genomes, but it will be difficult to tell which pathways best distinguish between our two groups of MAGs. For that task, we need the help of a statistical test.

### Finding enriched metabolic pathways

Anvi'o has a program for computing enrichment of metabolic modules in different groups of genomes, and that program is {% include PROGRAM name="anvi-compute-metabolic-enrichment" %}. It will compute an enrichment score and a list of associated groups for each module that is present in at least one genome (modules are considered 'present' in a genome if they have a high enough completeness score in that genome).

To run this program, you must provide it with the modules mode output file we generated in the last section, as well as a {% include ARTIFACT name="groups-txt" %} file that matches each genome to its group name. The latter file is provided in the datapack, so you don't need to generate it yourself. Here is the code to run the enrichment program:

```bash
anvi-compute-metabolic-enrichment -M FMT_MAG_metabolism_modules.txt \
                                  -G MAG_groups.txt \
                                  -o metabolic-enrichment.txt
```

The result will be a {% include ARTIFACT name="functional-enrichment-txt" %} file describing the name, enrichment score, associated groups, and other information about each metabolic module.

Here are the first 10 lines of this file (scroll to the right to see more columns):

KEGG_MODULE | enrichment_score | unadjusted_p_value | adjusted_q_value | associated_groups | accession | sample_ids | p_LOW_FITNESS | N_LOW_FITNESS | p_HIG_FITNESS | N_HIG_FITNESS
:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
Isoleucine biosynthesis, threonine => 2-oxobutanoate => isoleucine | 22.556406615904528 | 2.0406287012176094e-6 | 1.018150573124383e-5 | HIG_FITNESS | M00570 | KC_MAG_00002,KC_MAG_00007,KC_MAG_00022,KC_MAG_00051,KC_MAG_00055,KC_MAG_00061,KC_MAG_00062,KC_MAG_00080,KC_MAG_00101,KC_MAG_00110,KC_MAG_00120,KC_MAG_00126,KC_MAG_00137,KC_MAG_00143,KC_MAG_00145,KC_MAG_00147,KC_MAG_00151,KC_MAG_00155,KC_MAG_00161,KC_MAG_00176,KC_MAG_00178 | 0.15 | 20 | 0.9 | 20
Valine/isoleucine biosynthesis, pyruvate => valine / 2-oxobutanoate => isoleucine | 22.556406615904528 | 2.0406287012176094e-6 | 1.018150573124383e-5 | HIG_FITNESS | M00019 | KC_MAG_00002,KC_MAG_00007,KC_MAG_00022,KC_MAG_00051,KC_MAG_00055,KC_MAG_00061,KC_MAG_00062,KC_MAG_00080,KC_MAG_00101,KC_MAG_00110,KC_MAG_00120,KC_MAG_00126,KC_MAG_00137,KC_MAG_00143,KC_MAG_00145,KC_MAG_00147,KC_MAG_00151,KC_MAG_00155,KC_MAG_00161,KC_MAG_00176,KC_MAG_00178 | 0.15 | 20 | 0.9 | 20
Adenine ribonucleotide biosynthesis, IMP => ADP,ATP | 21.53846258298333 | 3.4680278204806215e-6 | 1.1535577282062716e-5 | HIG_FITNESS | M00049 | KC_MAG_00002,KC_MAG_00007,KC_MAG_00022,KC_MAG_00051,KC_MAG_00055,KC_MAG_00057,KC_MAG_00061,KC_MAG_00062,KC_MAG_00080,KC_MAG_00093,KC_MAG_00101,KC_MAG_00110,KC_MAG_00116,KC_MAG_00120,KC_MAG_00121,KC_MAG_00126,KC_MAG_00137,KC_MAG_00143,KC_MAG_00145,KC_MAG_00147,KC_MAG_00151,KC_MAG_00155,KC_MAG_00161,KC_MAG_00162,KC_MAG_00176,KC_MAG_00178 | 0.3 | 20 | 1 | 20
Pentose phosphate pathway, non-oxidative phase, fructose 6P => ribose 5P | 20.416675803269 | 6.22846903390889e-6 | 1.5538150847280372e-5 | HIG_FITNESS | M00007 | KC_MAG_00002,KC_MAG_00007,KC_MAG_00017,KC_MAG_00022,KC_MAG_00051,KC_MAG_00055,KC_MAG_00061,KC_MAG_00062,KC_MAG_00080,KC_MAG_00101,KC_MAG_00110,KC_MAG_00116,KC_MAG_00120,KC_MAG_00121,KC_MAG_00126,KC_MAG_00137,KC_MAG_00143,KC_MAG_00145,KC_MAG_00147,KC_MAG_00151,KC_MAG_00155,KC_MAG_00161,KC_MAG_00162,KC_MAG_00178 | 0.25 | 20 | 0.95 | 20
C5 isoprenoid biosynthesis, non-mevalonate pathway | 18.026729736369713 | 2.178249150838522e-5 | 3.6227162410137385e-5 | HIG_FITNESS | M00096 | KC_MAG_00002,KC_MAG_00007,KC_MAG_00017,KC_MAG_00022,KC_MAG_00051,KC_MAG_00055,KC_MAG_00057,KC_MAG_00061,KC_MAG_00062,KC_MAG_00080,KC_MAG_00101,KC_MAG_00110,KC_MAG_00116,KC_MAG_00121,KC_MAG_00126,KC_MAG_00137,KC_MAG_00143,KC_MAG_00145,KC_MAG_00147,KC_MAG_00151,KC_MAG_00155,KC_MAG_00161,KC_MAG_00162,KC_MAG_00176,KC_MAG_00178 | 0.3 | 20 | 0.95 | 20
Inosine monophosphate biosynthesis, PRPP + glutamine => IMP | 18.026729736369713 | 2.178249150838522e-5 | 3.6227162410137385e-5 | HIG_FITNESS | M00048 | KC_MAG_00002,KC_MAG_00007,KC_MAG_00022,KC_MAG_00051,KC_MAG_00055,KC_MAG_00057,KC_MAG_00061,KC_MAG_00062,KC_MAG_00080,KC_MAG_00093,KC_MAG_00101,KC_MAG_00110,KC_MAG_00116,KC_MAG_00120,KC_MAG_00121,KC_MAG_00126,KC_MAG_00137,KC_MAG_00143,KC_MAG_00147,KC_MAG_00151,KC_MAG_00155,KC_MAG_00161,KC_MAG_00162,KC_MAG_00176,KC_MAG_00178 | 0.3 | 20 | 0.95 | 20
F-type ATPase, prokaryotes and chloroplasts | 17.289003389888027 | 3.210393689132397e-5 | 4.5765505959647144e-5 | HIG_FITNESS | M00157 | KC_MAG_00007,KC_MAG_00017,KC_MAG_00022,KC_MAG_00051,KC_MAG_00055,KC_MAG_00057,KC_MAG_00061,KC_MAG_00062,KC_MAG_00080,KC_MAG_00093,KC_MAG_00110,KC_MAG_00116,KC_MAG_00120,KC_MAG_00126,KC_MAG_00137,KC_MAG_00143,KC_MAG_00145,KC_MAG_00147,KC_MAG_00151,KC_MAG_00155,KC_MAG_00161,KC_MAG_00176,KC_MAG_00178 | 0.25 | 20 | 0.9 | 20
Coenzyme A biosynthesis, pantothenate => CoA | 15.824346512060853 | 6.950242168239118e-5 | 8.669378513988786e-5 | HIG_FITNESS | M00120 | KC_MAG_00002,KC_MAG_00007,KC_MAG_00017,KC_MAG_00022,KC_MAG_00051,KC_MAG_00055,KC_MAG_00057,KC_MAG_00061,KC_MAG_00062,KC_MAG_00080,KC_MAG_00093,KC_MAG_00101,KC_MAG_00110,KC_MAG_00116,KC_MAG_00120,KC_MAG_00121,KC_MAG_00126,KC_MAG_00137,KC_MAG_00143,KC_MAG_00145,KC_MAG_00147,KC_MAG_00151,KC_MAG_00155,KC_MAG_00161,KC_MAG_00162,KC_MAG_00178 | 0.35 | 20 | 0.95 | 20
Guanine ribonucleotide biosynthesis IMP => GDP,GTP | 15.172414714793092 | 9.812648253560582e-5 | 1.0728293795381548e-4 | HIG_FITNESS | M00050 | KC_MAG_00002,KC_MAG_00007,KC_MAG_00017,KC_MAG_00022,KC_MAG_00051,KC_MAG_00055,KC_MAG_00057,KC_MAG_00061,KC_MAG_00062,KC_MAG_00080,KC_MAG_00093,KC_MAG_00101,KC_MAG_00110,KC_MAG_00116,KC_MAG_00120,KC_MAG_00121,KC_MAG_00122,KC_MAG_00126,KC_MAG_00137,KC_MAG_00143,KC_MAG_00145,KC_MAG_00147,KC_MAG_00151,KC_MAG_00155,KC_MAG_00157,KC_MAG_00161,KC_MAG_00162,KC_MAG_00176,KC_MAG_00178 | 0.45 | 20 | 1 | 20

The modules are organized so that those with higher enrichment scores (and lower significance values) are at the top. You can filter the output using the `unadjusted_p_value` or `adjusted_q_value` columns to make sure you only keep the modules that are most enriched in one group or another (the `adjusted_q_value` column is arguably the best one to filter with as this significance value is adjusted for multiple hypothesis testing). For example, in our study we considered any metabolic module with a q-value less than 0.05 to be enriched in its associated group, as long as it was also at least 75% complete in at least 50% of the group members.

The modules in the partial output above represent the metabolic pathways with the highest enrichment scores in our MAGs. All of the modules shown passed our filtering criteria - they have q-values less than 0.05 and are present in at least 10/20 of their associated group. In total, our group found 33 different KEGG modules that were enriched in these MAGs, including the above 9. You can find out which modules these are by taking a look at our [supplementary table 7](https://figshare.com/articles/dataset/Supplementary_Tables/14138405/2?file=26827175), sheet (d).

You will see that all 9 of the above modules are enriched in the "high-fitness" group of MAGs - and indeed, all 33 modules that passed our filters were enriched in this group. You might also notice that most of the enriched modules are biosynthesis pathways, particularly of amino acids and essential cofactors. As we discuss in the paper, this indicated to us that "high-fitness" populations were successful in colonizing the FMT recipients because they were metabolically independent - that is, they were able to individually produce the molecules they needed to survive and grow, and thus had a competitive advantage compared to the "low-fitness" populations, which generally did not have these biosynthesis capabilities.

## Conclusion

This metabolism estimation and enrichment analysis allowed us to form a clear hypothesis as to why some microbial populations were able to colonize FMT recipients while others weren't. For a much more robust discussion of the analysis and our conclusions, please read our [study](https://doi.org/10.1101/2021.03.02.433653).

We hope the mini-tutorial above helped you to learn how to run metabolism analyses in anvi'o on your own data. If anything was not clear, please feel free to comment below or reach out to us on Discord with questions. Thanks for reading!

{% include _join-anvio-discord.html %}

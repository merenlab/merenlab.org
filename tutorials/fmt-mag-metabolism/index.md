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

This is a **mini-tutorial for the metabolism suite of programs** in anvi'o. We'll be going through how to estimate metabolism and compute enrichment scores for metabolic modules in two different groups of MAGs.

The data we'll be using for this is a real dataset from one of our recent studies, ["Metabolic competency drives microbial colonization and resilience in health and disease”](https://doi.org/10.1101/2021.03.02.433653) by Watson et al. In fact, this post is doing double-duty as a reproducible workflow for one of the analyses in that study. :) The rest of the reproducible workflow can be found [here](https://merenlab.org/data/fmt-gut-colonization/) for anyone who is interested in how we did the other analyses discussed in the paper.

{:.notice}
This tutorial is tailored for anvi'o `v7.1` or later. You can learn the version of your installation by running `anvi-interactive -v` in your terminal.

## A dataset of high- and low-fitness MAGs

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

One way to answer this question is to look at the metabolic potential, or genomically-encoded metabolic capabilities, of these MAGs.

## Estimating metabolism for these MAGs

The program {% include PROGRAM name="anvi-estimate-metabolism" %} computes the completeness of metabolic pathways in genomes, MAGs, or metagenomes by matching gene annotations from the [KOfam database](https://academic.oup.com/bioinformatics/article/36/7/2251/5631907) against [KEGG definitions of metabolic modules](https://www.genome.jp/kegg/module.html). You can run it on all 40 MAGs in this dataset at once by utilizing the {% include ARTIFACT name="external-genomes" text="external genomes file" %} provided in the datapack, as shown below:

```bash
anvi-estimate-metabolism -e external-genomes.txt -O FMT_MAG_metabolism
```

This will give you one [modules mode](https://merenlab.org/software/anvio/help/main/artifacts/kegg-metabolism/#modules-mode) output file called `FMT_MAG_metabolism_modules.txt`, which describes the completeness of each KEGG Module in each MAG.

You could look through this file manually to see what metabolisms are encoded in these genomes, but it will be difficult to tell which pathways best distinguish between our two groups of MAGs. For that task, we need the help of a statistical test.

## Finding enriched metabolic pathways

Anvi'o has a program for computing enrichment of metabolic modules in different groups of genomes, and that program is {% include PROGRAM name="anvi-compute-metabolic-enrichment" %}. It will compute an enrichment score and a list of associated groups for each module that is present in at least one genome (modules are considered 'present' in a genome if they have a high enough completeness score in that genome).

To run this program, you must provide it with the modules mode output file we generated in the last section, as well as a {% include ARTIFACT name="groups-txt" %} file that matches each genome to its group name. The latter file is provided in the datapack, so you don't need to generate it yourself. Here is the code to run the enrichment program:

```bash
anvi-compute-metabolic-enrichment -M FMT_MAG_metabolism_modules.txt -G MAG_groups.txt -o metabolic-enrichment.txt
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

We hope the mini-tutorial above helped you to learn how to run metabolism analyses in anvi'o on your own data. If anything was not clear, please feel free to comment below or reach out to us on Slack with questions. Thanks for reading!

{% include _join-anvio-slack.html %}

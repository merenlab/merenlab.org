---
layout: post
title: Structural biology with anvi'o
modified: 2018-09-04
excerpt: "Combining metagenomic sequence data with modelled protein structures"
comments: true
authors: [evan]
categories: [anvio]
---

{% capture images %}{{site.url}}/images/anvio/2018-09-04-structural-biology-with-anvio{% endcapture %}

{% include _toc.html %}

{:.notice}
The contents of this post will only work with anvi'o `v5` or later.

{:.notice}
This is a theoretical tutorial describing the details of the structure database, how it's created, and its utility. For a more practical tutorial that demonstrates some of the same concepts, and a thorough walkthrough of the `anvi-display-structure` interactive interface, please visit [http://merenlab.org/tutorials/infant-gut/#chapter-vii-linking-genomic-heterogeneity-to-protein-structures](http://merenlab.org/tutorials/infant-gut/#chapter-vii-linking-genomic-heterogeneity-to-protein-structures).

{% include _join-anvio-slack.html %}

To me, the most fascinating thing about proteins is how specifically they associate with their ligands to carry out their function, and yet how they all kinda look the same. An alpha helix here, a beta sheet there, couple of loops connecting the two, but overall just a big glob of polypeptide chain. The arbitrariness of protein shape is the reason that it's so surprising to me that these proteins don't just promiscuously interact with everything around themselves, aggregating to become toxic biohazards inside the cell. Yet of course we know they don't do that. We know they are finely tuned by natural selection to interact specifically with their ligands and that even a single mutation in the right place can destroy this functionality, leading to a potentially detrimental fitness decrease of the cell and elimination from the population. This is life.

So it's true that proteins are extremely sensitive to certain mutations. But they are surprisingly robust to others, and it is this permissibility of neutral mutations that leads to the majority of amino acid and codon sequence diversity observed in microbial populations. This tutorial will introduce how you can use full-genome shotgun metagenomics to explore insights into the relationships between the sequence heterogeneity in microbial populations, how it manifests at the level of translated protein, and how it is shaped by principles of evolutionary biochemistry.

The meat of this workflow is just two anvi'o commands. The first command predicts protein structures encoded by genes present in a contigs database, and the resulting structures are stored in a structure database. This is accomplished with the program `anvi-gen-structure-database`, cousin to `anvi-gen-contigs-database`, a command you may be more familiar with. Once this is done, variation from a profile database can be visualized onto the structures with the command `anvi-display-structure`. From within the interface one can interactively filter and visualize variants overlaid on the structure. It's really that easy.

The itinerary for this post is to first discuss how to set up your computational environment. Afterwards, the structure database is discussed in depth. Finally applications of the structure database, including interactively visualizing protein variants, is discussed.

# Dependencies you may need to install

This workflow depends heavily on three dependencies, MODELLER, DSSP, and NGL. If you use this workflow, make sure you cite them.

## MODELLER

{:.notice}
**Citation**: [doi:10.1002/0471250953.bi0506s15](https://doi.org/10.1002/0471250953.bi0506s15)

{:.notice}
**Citation**: [doi:10.1146/annurev.biophys.29.1.291](https://doi.org/10.1146/annurev.biophys.29.1.291)

{:.notice}
**Citation**: [doi:10.1006/jmbi.1993.1626](https://doi.org/10.1006/jmbi.1993.1626)

MODELLER is the program anvi'o uses to predict protein structure based on experimentally solved structures in the Protein Data Bank. We'll talk more specifically about how it accomplishes this in the following section, but for now you need to make sure it's installed on your computer. For that, check out these instructions to see if you have it installed ([click me](http://merenlab.org/2016/06/18/installing-third-party-software/#modeller)), and how to install it if you don't. We've tried to make it as simple for you as possible.

## DSSP

{:.notice}
**Citation**: [doi:10.1002/bip.360221211](http://doi.org/10.1002/bip.360221211),

{:.notice}
**Citation**: [doi:10.1093/nar/gkq1105](https://doi.org/10.1093/nar/gkq1105)

DSSP (Dictionary of Secondary Structure Prediction) is the program anvi'o uses to assign secondary structure and other useful biophysical and biochemical characteristics for each residue in the predicted structure of MODELLER. It is not a strict requirement for the workflow, although you will be missing out if you don't install. Check out these instructions to see if you have it installed, and how to install it if you don't ([click me](http://merenlab.org/2016/06/18/installing-third-party-software/#dssp)).

## NGL

{:.notice}
**Citation**: [doi:10.1093/nar/gkv402](http://doi.org/10.1093/nar/gkv402)

{:.notice}
**Citation**: [doi:10.1093/bioinformatics/bty419](https://doi.org/10.1093/bioinformatics/bty419)

NGL (NGL) is an open-source project for visualizing biomolecules. This browser-based solution to visualization means you don't have to install anything, and you can thank them for that. We are continually being impressed with NGL's excellent code and documentation, as well as to their open-source approach to science.

# The structure database

This section provides an overview of how protein structures are predicted, the commands responsible for doing so, the user tunable parameters to refine how your structure database is created, and finally the effects they have on structure accuracy.

Before getting started, there's a disclaimer that has to be made. Us developers at anvi'o are not structural biologists, and we are not authorities on reference based protein homology modelling. We are just taking MODELLER, which is a long-standing golden standard in protein prediction, and making it accessible to microbial ecologists and their data. If you disagree with what were saying, we're probably wrong and would appreciate you correcting us.

## How MODELLER works

How does MODELLER work?  [This paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4186674/) dissects the MODELLER workflow into four major steps which are summarized here as I understand them.

1. Alignment: sequences from a database of proteins with known structure are aligned to your query protein (aka your target protein) sequence, and the database sequences with notably good alignments are stored as candidate templates. 

2. Template selection: from this list of candidate templates, some number are selected to be structural templates that your target structure will be modelled from. The selection is done either with automated criteria or ideally through meticulous, manual intervention. Note that several templates can be used to model your target structure.

3. 3D-alignment: your target sequence is aligned in 3D space to the known structures of your templates. You can think of this as "threading" your 1D sequence into a 3D space.

4. Satisfy spatial restraints: the alignment structure from step 3 is the initial configuration for a simulation which perturbs the atomic positions of your protein iteratively to solve an optimization problem. The optimization minimizes the likelihood of a many-dimensional probability density function that compares C alpha-C alpha distances, main-chain N-O distances, main-chain and side-chain dihedral angles in your protein to the distribution of values found in a large database of reference protein structures. The simulation is evolved until the structure settles into an equilibrium, and the final state is your predicted structure. If you are specifically interested in this step, you should also check out the original MODELLER paper [here](https://www.ncbi.nlm.nih.gov/pubmed/8254673).

If you think this is amazing, the scientific community agrees with you; MODELLER has been cited more times than QIIME, and has remained relevant to the community for over 25 years.

## Why MODELLER works

MODELLER and similar programs make use of protein homology modelling and reference structures to stitch together a predicted structure. This approach is successful for two main reasons:

The first is the surprising fact that structure is much more conserved than sequence. Let me ask you a question. Shown below are the structures of two nearly equal-length proteins I have overlaid on top of each other. Based on their obvious structural similarity, what percentage of amino acids would you guess are identical? To say it another way, if you move along the sequence of one protein, what percentage of residues would match with the other protein's residue at the same position? 

[![3dfr-vs-4dfr]({{images}}/3dfr-vs-4dfr.gif)]( {{images}}/3dfr-vs-4dfr.gif){:.center-img .width-50}

I'll let you think know about it...

These two proteins are [3DFR](https://www.rcsb.org/structure/3dfr) and [4DFR](https://www.rcsb.org/structure/4dfr) and you may be surprised to know that they only share 29 percent sequence identity, with a root-mean square difference between their main chain atoms (RMSD) of 3.927 Ångströms. Visually, you can see that's not much. In the paper from which I found these proteins, they looked at a total of 34 pairs of homologous proteins (which was totally badass in 1986) and they found the following relationship between sequence similarity and RMSD:

[![sequence-vs-structure-divergence]({{images}}/sequence-vs-structure-divergence.png)]( {{images}}/sequence-vs-structure-divergence.png){:.center-img .width-80}

In red I've circled the above comparison. The RMSD is lower than the previously reported value of 3.927 Ångströms because they only performed the calculation on a subset of the atoms. The results quite boldy show that RMSD is conserved even over significant sequence divergence. My interpretation of these results is that the "metastructure" for a protein family corresponds to an island in sequence space, and protein sequences in the family evolve happily within the confines of the island, which is more than spacious enough for functional diversification. Yet if any protein sequence drifts off of the island (i.e its structure becomes too disimilar from the metastructure) it faces deterimental fitness costs brought on by the open ocean (i.e. being unable to fold or remain stable) and is unable to cope. So virtually all sequences remain on the island, because its too difficult to get off, and hey, there's palm trees so things here aren't so bad.

Anyways, that's the first main reason why template-based homology modelling approaches work. The second reason is that there exists a lot of solved structures that can be used as templates to predict the structure of your target sequence. In my brief experience using E. Faecalis and SAR11 genomes, approximately half of the genes have sequence matches above 30%. However, please keep in mind the following law of the universe, which is that the more you want a gene to have a homologous structure, the less often it does.

Based on the results of [this MODELLER paper](https://www.ncbi.nlm.nih.gov/pubmed/18428767), as well as several protein experts I've talked to, 30% sequence similarity with few extended gap regions is around the lower limit of when structures can be modelled through homology and the results can be trusted to within a couple Ångströms. Obviously, drawing parallels between characteristics such as binding pocket geometry should be taken with a grain of salt at this level of similarity, as functions can differ drastically.


## The vanilla command 

So you want to predict some structures for some genes you've identified as interesting, eh? All you need is a list of those gene IDs and the contigs database they belong to. The most vanilla command is

```bash
anvi-gen-structure-database -c path/to/your/CONTIGS.db \
                            --gene-caller-ids 1,2,3 \
                            -o STRUCTURE.db
```

Where `path/to/your/CONTIGS.db` is the file path of your contigs database, `-o STRUCTURE.db` states that the database will be named `STRUCTURE.db`, and `1,2,3` are three genes you are interested in, identified by their gene caller IDs, and separated by commas (no spaces). Alternatively, you could make a file that is a list of genes caller IDs of interest:

```
1
2
3
```

And run the command

```bash
anvi-gen-structure-database -c path/to/your/CONTIGS.db \
                            --genes-of-interest path/to/your/genes/of/interest.txt \
                            -o STRUCTURE.db
```

Whichever the method, anvi'o will invoke MODELLER, and MODELLER will try its best to predict a structure. Everything from the templates used to the quality scores of the protein models will be stored in the structure database, along with the structures themselves.

## Description of all MODELLER parameters

Okay so you know how to make a vanilla structure database, but MODELLER has many parameters that are waiting to be refined. We have tried our best to make as many parameters tunable as possible to create the databases you want and with the level of stringency that you control. (If you are a structural biologist guru and there's something you want to control but that we don't provide access to, please let us know and we'll try our best to add it). Below are a list of all the parameters that you can use to fine-tune database creation.

1. **`--num-models`**. As [already described](#the-structure-database), once MODELLER has 3D-aligned the target to the templates, it carries out an optimization simulation. This process is deterministic, meaning that if you run the simulation 3 times with the same initial structures, the 3 final structures (aka models) will be the same. Yet if you perturb the atomic positions of the initial structure before running the simulations (see `--deviation`), you may stumble upon different, and perhaps even better, models. Therefore, it makes sense to simulate several models since more of the solution space will be searched. The number of models simulated is controlled with this parameter. It should be kept in mind that the largest determinant of a model's accuracy is determined by the protein templates used, so no need to go overboard with excessively large values for `--num-models`. More than 20 seems excessive. The default is 3.

2. **`--deviation`**. As mentioned in `--num-models`, you can explore more of the solution space by perturbing the atomic positions of the initial structure. This value represents the standard deviation of atomic perturbation of the initial structure in Ångströms. The default is 4.

3. **`--modeller-database`**. MODELLER finds candidate template sequences from a database of proteins with known structure. By default, a database containing all of the PDB structures clustered at 95% sequence identity is used. It's filename is `pdb_95.pir` and can be downloaded [here](https://salilab.org/modeller/supplemental.html). It's a very convenient resource and is periodically updated by the Sali Lab. If you don't have it anvi'o will download it for you when the time is right. If you have your own database it must have either the extension `.bin` or `.pir` and should be placed in the directory `anvio/data/misc/MODELLER/db`.

4. **`--scoring-method`**. If you generate 10 models with `--num-models 10`, how should the best model be decided? The metric used could be any of `GA341_score` ([citation](https://salilab.org/pdf/John_NucleicAcidsRes_2003.pdf)), `DOPE_score` ([citation](https://salilab.org/pdf/Shen_ProteinSci_2006.pdf)), or `molpdf`. `GA341` is an absolute measure, where a good model will have a score near `1.0`, whereas anything below `0.6` can be considered bad. `DOPE_score` and `molpdf` scores are relative energy measures, where lower scores are better. `DOPE` has been generally shown to be a better distinguisher between good and bad models than `molpdf`, which was the original assessment method used in the pre-American Civil War era. By default, `DOPE_score` is used. To learn more about assessment methods, see this MODELLER tutorial [here](https://salilab.org/modeller/tutorial/basic.html). 

5. **`--percent-identical-cutoff`**. If a protein in the database has a proper percent identity to the protein of interest that is greater than or equal to the `--percent-identical-cutoff`, then it is considered as a template. Otherwise it is not. Here we define proper percent identity as the percentage of amino acids in the gene of interest that are identical to an entry in the database given the sequence length of the protein of interest. For example, if there is 100% identity between the protein of interest and the template over the length of the alignment, but the alignment length is only half of the protein of interest sequence length, then the proper percent identical is 50%. This helps us avoid the inflation of identity scores due to only partially good matches. The default is 30. Obviously the higher the percentage identity, the more confident you can be in the accuracy of the model, however the less hits you will have. It's very important to point out that selecting template structures based on a single threshold is very simplistic if you are after high accuracy models using only nominally homologous templates. As mentioned, the primary determinant of a model's accuracy is which template structures are used, and so template selection is perhaps the most critical aspect of the workflow, and requires manual intervention if optimal accuracy is required.

6. **`--max-number-templates`**. Generally speaking it is best to use as many templates as possible given that they have high proper percent identity to the protein of interest. Here ([click me](https://salilab.org/modeller/methenz/andras/node4.html)) is an excerpt from the Sali Lab: 'The use of several templates generally increases the model accuracy. One strength of MODELLER is that it can combine information from multiple template structures, in two ways. First, multiple template structures may be aligned with different domains of the target, with little overlap between them, in which case the modeling procedure can construct a homology-based model of the whole target sequence. Second, the template structures may be aligned with the same part of the target, in which case the modeling procedure is likely to automatically build the model on the locally best template [43,44]. In general, it is frequently beneficial to include in the modeling process all the templates that differ substantially from each other, if they share approximately the same overall similarity to the target sequence.' The default is 5, but if only X candidate templates are found to pass the `--percent-identical-cutoff` threshold, then only X will be used.

7. **`--very-fast`**. Use this option if you're impatient and want only the roughest predicted structures. It's fast because the step size of the simulation his very large and a very low quality requirement is established for reaching equilibrium. Not recommended, but we're not going to the be the ones to deny you of freewill.

8. **`--dump-dir`**. Modelling structures requires a lot of moving parts, each which have their own outputs. The output of this program is a structure database containing the pertinent results of this computation, however a lot of stuff doesn't make the cut. By providing a directory for this parameter you will get, in addition to the structure database, a directory containing the raw output for everything produced by MODELLER.

## A quick case study on the importance of key parameters

How much do these parameters matter? In the following we look at a CDP-D-glucose 4,6-dehydratase from a SAR11 genome, and the effect template selection and model number have on structure prediction of this protein.

### How much does the number of models matter?

Hopefully it's clear that increasing the number of models increases the search space and increases the accuracy of your final model, but when is enough, enough? To test this I predicted the SAR11 dehydratase structure with default parameters except with `--num-models 100`. Here are the DOPE scores of those 100 structure predictions, which took about an hour to compute:

[![100-models-hist]({{images}}/100-models-hist.png)]( {{images}}/100-models-hist.png){:.center-img .width-80}

Since model 79 has the best DOPE score of -42376.1132812 (most negative), it was picked as the best model. But that's just a number, and I certainly don't know how to meaningfully compare models using this score besides comparing their magnitudes. So here is a visual comparison between the best and worst models according to DOPE score:

[![num-model-comparison]({{images}}/num-model-comparison.png)]( {{images}}/num-model-comparison.png){:.center-img .width-80}

They have an RMSD of 1.617. I will let you be the judge of how different these are, but in general keep in mind it will depend on what questions you are interested in. Based on their visual similarity the default value is `--num-models 3`.

### How much do templates matter?

It's clear that the higher the similarity your templates have to your target, the better your model will be. But how does the effect of multiple templates influence model results? Let's see.

For this SAR11 dehydratase, MODELLER finds 2 candidate templates above 30 proper percent identity:

```
Template 1 ........................: Protein ID: 1wvg, Chain A (39.4% identical)
Template 2 ........................: Protein ID: 1rkx, Chain A (38.8% identical)
```

Template 1 comes from *Salmonella typhi* and template 2 from *Yersinia pseudotuberculosis*. They are both ~39% similar to SAR11 but 73% similar to one another, with a RMSD of 2.053 Ångströms between their main-chain atoms. I created 3 different models, one using template 1, a second using template 2, and a third using templates 1 and 2. The best models (using `--num-models 15`) from each of these 3 template selections were selected and their structures aligned against one another and shown below.

[![gene-1248-different-template-selections]({{images}}/gene-1248-different-template-selections.png)]( {{images}}/gene-1248-different-template-selections.png){:.center-img .width-80}

Visually it's clear that there are differences between these models that could be considered either very significant, or very insignificant--it depends what you're interested in. If you want to classify which residues are in the hydrophobic core versus the surface of the protein, these are all essentially equivalent models. If on the other hand, you are trying to define the precise geometry of a binding pocket, these differences will be critical to understand. You should also probably stop trying to do this with 39% templates and potentiallly reconsider your life decisions. 

Here is a table quantifying RMSD between each pair (in Ångströms):

Templates used |  1rkx  |  1wvg   |  both
---------------|--------|---------|-------
1rkx           |  -     |  4.429  |  4.304
1wvg           |  -     |  -      |  2.617
both           |  -     |  -      |  -

So which model is best? Looking at the DOPE scores,

Templates used |   DOPE   |
---------------|----------|
1rkx           |  -42248  |
1wvg           |  -42209  |
both           |  -42362  |

We see that the model based on both templates outperforms both single-template models. This reflects a general rule of thumb proposed by the Sali Lab: 

>The use of several templates generally increases the model accuracy. One strength of MODELLER is that it can combine information from multiple template structures ... the template structures may be aligned with the same part of the target, in which case the modeling procedure is likely to automatically build the model on the locally best template. In general, it is frequently beneficial to include in the modeling process all the templates that differ substantially from each other, if they share approximately the same overall similarity to the target sequence.

In our case, both the salmonella and the pseudotuberculosis proteins shared the same overall similarity to the SAR11 sequence, so utilizing both templates led to a better model. What if there was a third template that shared 90% similarity to our protein? Should we use all three templates or should we just use the one with 90% similarity? In this instance, we should use only the one template with 90% similarity, since the other 2 templates do not share the same overall similarity to the target sequence, and would 'contaminate' our high accuracy model based on 90% sequence similarity with information of lower quality. If on the other hand the third model was 45-50% similar to the target sequence, I would probably use all three.

# Using the structure database as an input to anvi-gen-profile-variability

Once you have structures for proteins of interest, what can you do with them? One utility is to annotate sequence variants generated during `anvi-gen-variability-profile` with structural information. 

`anvi-gen-profile-variability` is a robust program to generate variability profiles for metagenomes at the level of nucleotides (SNVs), codons (SCVs), or amino acids (SAAVs). You can get a hands-on feel for this program in this section of the infant gut tutorial ([click me](http://merenlab.org/tutorials/infant-gut/#chapter-iv-microbial-population-genetics)) or get into the nitty-gritty in this more theoretical tutorial dedicated to the subject ([click me](http://merenlab.org/2015/07/20/analyzing-variability/)).

The output of `anvi-gen-profile-variability` is a table where each row corresponds to a sequence variant found in a metagenome. Each column gives descriptive or quantitative information about the variants. If the flag `-s YOUR/STRUCTURE/DATABASE.db` is given to the program, additional columns will be added for structural information if a) the variant is contained within a gene, b) the gene has a predicted structure and c) the variant is either a SAAV (use `--engine AA`) or a SCV (use `--engine CDN`). All the columns that can be added are listed [here](http://merenlab.org/2015/07/20/analyzing-variability/#matrix-output-those-unique-to-structure-database-integration). This is my favourite feature of the structure database.

# Display metagenomic sequence variants directly on predicted structures

{:.notice}
The anvi'o program `anvi-display-structure` described below visualizes variability from metagenomic sequence data onto a reference 3D structures. The structures are predicted either from a MAG or a reference genome predicted during `anvi-gen-structure-database`. This program unfortunately does not predict how such variants will change the 3D structure.

There are a lot of amazing protein visualization software packages and web servers in existence. In comparison to these projects that have innumerable features, which makes them capable of doing anything you can dream of, we wanted to create something with the very specific purpose of interactively visualizing, filtering, and clustering metagenomic sequence variants directly onto the proteins they encode. These aforementioned programs are certainly capable of a job like this--I know because I've done it myself--but it took me forever to do it for just one protein. Moreover, the lack of automation, scalability, and just the straight up time required for a single protein really impeded me from tackling the proteome-level research questions I was interested in, and I realized that anyone interested in similar questions as me could really benefit from a tool designed specifically for this purpose. That tool is what we've tried to create with `anvi-display-structure`.

## Examples

Details describing the interface itself are provided in a hands-on demonstration within the infant gut tutorial ([click me](http://merenlab.org/tutorials/infant-gut/#chapter-vii-linking-genomic-heterogeneity-to-protein-structures)), which I highly recommend you follow along with, either with the provided data or your own. But the tutorial you're currently reading wouldn't be complete if we didn't at least showcase some images and movies from the interface. With this in mind, check out this video demonstration ([click me](https://www.youtube.com/watch?v=kHHk1qUOYoE)) of some of the interface features. As well, feel free to view the screenshots below with their accompanying descriptions.

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2018-09-04-structural-biology-with-anvio/example.png"><img src="{{ site.url }}/images/anvio/2018-09-04-structural-biology-with-anvio/example_structure1.png" style="border: none; width: 100%;" /></a>
</div>

Above is a protein from a reference SAR11 genome called HIMB83. Each of the 4 views correspond to groups of TARA ocean metagenomes organized by the user. Each sphere on the structure corresponds to a SCV that occurred in at least one of the metagenomes of the group. The size of each sphere relates to the solvent accessibility of the residue (aka how exposed to water it is) and the color of each sphere corresponds to the SCV's synonymity. The trend demonstrates that solvent inaccessible residues are under high purifying selection against non-synonymous mutations, but readily permit synonymous mutations as they do not compromise protein stability.

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2018-09-04-structural-biology-with-anvio/example.png"><img src="{{ site.url }}/images/anvio/2018-09-04-structural-biology-with-anvio/example_structure2.png" style="border: none; width: 100%;" /></a>
</div>

Above is a merged view of all TARA ocean metagenomes for the same protein. The size of the spheres refer to the fraction of metagenomes in which the position was a SCV (e.g. the smallest spheres occur in only one metagenome). We see that one region of this protein is significantly more variable than the other. The color refers to the coverage of that site normalized by the average gene coverage for a metagenome. This indicates that the hypervariable region is also recruiting much more reads than the others, as a result of non-specific mapping.

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2018-09-04-structural-biology-with-anvio/example.png"><img src="{{ site.url }}/images/anvio/2018-09-04-structural-biology-with-anvio/example_structure3.png" style="border: none; width: 100%;" /></a>
</div>

These protein views are interactive with many visualization options. For example, above we have zoomed in to a specific region and visualized the side chains of SCVs that are in close proximity to one another.

## Supplying anvi-display-structure with sequence variability

{:.warning}
Please note that the only required input for creating a structure database with `anvi-gen-structure-database` is a contigs database. But to visualize sequence variants on structures in your structure database, you need a merged (many metagenomes) or single (one metagenome) profile database. Furthermore, when the profile database(s) were created with `anvi-profile`, they **must** have been supplied the flag `--profile-SCVs`, as this flag is what makes reporting SCVs and SAAVs possible. Otherwise, anvi'o will complain, and complaining has never solved any problem. If you have profile databases without this flag and want to visualize variants on structures, you have to re-profile with `--profile-SCVs`. We apologize, but the reason this is not done by default is that it comes at a large computational cost.

The remainder of this tutorial describes the two ways in which you can provide anvi'o with sequence variability, which is required to run the interface. The easiest and most straightforward way is to *not*, and then anvi'o will do it for you. Simply provide the corresponding contigs and profile databases alongside your structure database:

```bash
anvi-display-structure -s PATH/TO/STRUCTURE/DB \
                       -p PATH/TO/PROFILE/DB \
                       -c PATH/TO/CONTIGS/DB
```

Anvi'o will calculate both SCVs and SAAVs on-the-fly and the following messages will appear in your terminal:

```
* SAAVs for gene X are loaded
* SCVs for gene X are loaded
```

Afterwards, the interface will open a display for gene X. Notice that SCVs and SAAVs were only calculated for gene X. If you switch genes in the interface, SCVs and SAAVs will be calculated for the new gene.

This is the standard way to visualize variation and is recommended in most cases. But there is an additional way that could be useful for one of two reasons: 1) For extremely large profiles, profiling variation on-the-fly can take an annoyingly long time. With a 26 GB profile database, it takes about a minute to load each gene. 2) If you merely want to show someone variation in a gene across some metagenomes, you don't want to send them a 26 GB profile database. Why not just send them a variation table? Hence, the second way to provide anvi'o with sequence variability is with a variability table generated from `anvi-gen-variability-profile`. With this approach, anvi'o will generate variability from this table instead of reading it from the profile database. The approach should look something like this:

```bash
# generate a variability table
anvi-gen-variability-profile -p PATH/TO/PROFILE/DB \
                             -c PATH/TO/CONTIGS/DB \
                             -s PATH/TO/STRUCTURE/DB \
                             --engine AA \
                             --only-if-structure \
                             -o variability.txt

# run the display using the variability table
anvi-display-structure -s PATH/TO/STRUCTURE/DB \
                       -V variability.txt
```

The first command has some bells and whistles that are worth describing. First, `-s PATH/TO/STRUCTURE/DB` is provided so that structural information is annotated within the variability table, as described [here](FIXME). Second, `--engine AA` indicates that SAAVs will be reported. Use `CDN` if you want to visualize SCVs, and `NT` for SNVs. Thirdly, `--only-if-structure` subsets your genes of interest to only include those with a structure in the structure database. That means if you have 1000 genes in your contigs database but only 3 have structure, variability will only be generated for those 3. Finally, `-o variability.txt` is the name of the output file.

Once `variability.txt` is generated, the second command feeds it into `anvi-display-structure` along with your structure database. A warning will pop up that is reiterated here:

```
WARNING
===============================================
You opted to work with a variability table previously generated from anvi-gen-
variability-profile. As a word of caution, keep in mind that any filters applied
when the table was generated now persist in the following visualizations.
```

# Conclusion

We hope that this tutorial has been useful for you. Something you feel is missing? If you have any suggestions or questions please contact us. If you discover any bugs, please make an issue on our [github](https://github.com/merenlab/anvio).


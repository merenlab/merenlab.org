---
layout: page
title: A reproducible workflow for Füssel et al. 2025
modified: 2025-06-09
excerpt: "A bioinformatics workflow for genomically guided compound prediction in Roseobacter co-culture experiments"
comments: true
authors: [sam]
---

<div class="extra-info" markdown="1">

<span class="extra-info-header">Summary</span>

**The purpose of this page** is to provide access to our bioinformatics workflow that predicted compound identifications of molecular features in our study titled "**Bacterial interactions shape the molecular composition of dissolved organic matter**" by Füssel et al.

In addition to providing transparency in our methods, this workflow can be used as the basis for genomically guided compound prediction in other metabolomics experiments, which will also help refine and validate the approach.

</div>

{:.notice}
If you have any questions, notice an issue, and/or are unable to find an important piece of information here, please feel free to leave a comment down below, send an e-mail to [us](/people/), or get in touch with us through Discord:

{% include _join-anvio-discord.html %}

## Study description

### Background

The carbon cycle in the surface ocean remineralizes the vast majority of fixed carbon on a time scale less than ~1 year, with the residual ~1% entering the longer-lived reservoir of dissolved organic matter (DOM) that is comparable in size to atmospheric carbon. The molecular composition of this DOM pool is similar in the surface and deep ocean and arises from the collective action of the microbiome. The exometabolomes of isolate cultures bear little resemblance to ocean DOM, which prompted us to study co-cultures of marine microbes to elucidate the process of DOM formation.

### Cultures

We studied factorial co-cultures of four coastal Marine _Roseobacter_ Group isolates obtained from the same seawater sample from the North Sea. The strains were grown individually and in co-cultures of two, three, and all four strains, each treatment done in triplicate. The artificial seawater minimal medium used in the experiment contained 1 g/L of glucose, trace elements, vitamins, and a bicarbonate buffer.

The four strains were found to have divergent metabolic capabilities from their genomes and substrate utilization preferences in culture. *Pelagimonas varians* SH4-1 (**SH4**) has a more extensive set of sugar metabolism genes than the other three strains, and grew on a variety of organic acids and monosaccharides, as well as a few polysaccharides. In the glucose minimal medium, SH4 had a negligible lag phase and grew to a higher optical density compared to the other strains. *Phaeobacter* sp. **SH40** and *Sulfitobacter* sp. SH22-1 (**SH22**) grew well on organic acids and relatively poorly on sugars, and *Sulfitobacter* sp. SH24-1b (**SH24**) exhibited limited growth on all tested substrates. After the longer lag phase compared to SH4, the three other strains also grew more slowly to stationary phase.

Growth in co-culture contrasted with growth in monoculture. Growth curves exceeded modeled curves with competitive glucose consumption, especially in co-cultures with SH4. The more strains in a co-culture, the greater the overall growth, suggesting that exometabolites produced by one strain were beneficial substrates to others.

### Untargeted metabolomics

Exometabolomes were measured in duplicate from the supernatant of each biological replicate at the beginning of the experiment and after 255 hours. DOM was extracted via Priority PolLutant SPE cartridges, which preferentially retain hydrophobic organic compounds. The analytes were then measured by FT-ICR-MS in negative ion mode using electrospray ionization. The mass error was <0.1 ppm for all samples following calibration to endogenous peaks. Only masses detected in all replicates of a culture and not present in blanks were retained. Molecular formulas were assigned to spectra by [ICBM-OCEAN](https://pubs.acs.org/doi/10.1021/acs.analchem.9b05659) software.

The fate of molecular formulas was tracked from pure culture to co-cultures. The fastest growing strain, SH4, yielded 2,216 formulas, or 89% of formulas found in pure cultures of the four strains. Across all co-cultures, 2,066 formulas were also found in pure cultures, while 2,508 were not. A majority of formulas novel to the co-cultures were unique to a single co-culture.

### Compound prediction

A formula can represent various isomers, so we used the metabolic networks predicted for each of the strains and consortia to propose molecular identifications for the formulas, as described in this workflow. This approach involves [anvi'o reaction networks](https://anvio.org/help/main/artifacts/reaction-network/), which are constructed from [KEGG Ortholog](https://www.genome.jp/kegg/ko.html) (KO) [annotations of genes](https://anvio.org/help/main/programs/anvi-run-kegg-kofams/) and associated reaction and compound entries from the [ModelSEED Biochemistry Database](https://github.com/ModelSEED/ModelSEEDDatabase). KOs are often annotated with KEGG reactions and EC numbers, indicating potential reactions that may be catalyzed by a gene protein product. Genomic reaction networks of co-cultured strains were merged to produce networks representing the joint metabolic potential.

For each molecular feature in a culture, we matched its neutral formula, formula with one subtracted hydrogen and charge of -1, and formula with two subtracted hydrogens and charge of -2 to the formulas of compounds in the culture reaction network. The network often contains compounds in the charge state that would exist in aqueous solution, so it is necessary to also search for -1 and -2 charged variants of the neutral formula to capture metabolites such as mono- and dicarboxylates.

#### Criteria

A set of criteria is used to evaluate the validity of formula matches and filter possible compounds. Some of the criteria are implemented automatically while others are broader, requiring careful consideration.

##### Multiple compound matches

A formula can match multiple compounds in a reaction network, but the strength of the evidence supporting each match may vary. It can be useful to retain multiple matching compounds when they are closely related metabolites, such as interconverted isomers occurring in the same metabolic pathway. Otherwise, formulas with multiple matches are ignored. The search for versions of each formula with charges of 0, -1, and -2 increases the likelihood of uncertain matches to multiple compounds, causing the formula to be filtered out.

##### Compound consistency across cultures

If a formula is from multiple cultures, it must match the same compound in all of the cultures' reaction networks. If the formula is from cultures A and B, but the matching compound is only in the reaction network of culture A, then the compound match would be ignored. Likewise, if compound X is found in culture A but not culture B, and isomeric compound Y is found in culture B but not culture A, then formula matches to these compounds would be filtered out.

##### Specificity of reaction annotations

To evaluate the validity of a compound match, the basis of the inclusion of the compound in the reaction network must also be understood. Compounds can be included due to their involvement in extremely broad and therefore uncertain categories of reactions. For example, KOs annotated with EC 1.1.1.1 (alcohol dehydrogenase) result in the addition of numerous ModelSEED alcohol dehydrogenase reactions involving various compounds to the reaction network. Formula matches to compounds from permissive reaction annotations are ameliorated by the following considerations.

- ModelSEED reactions included on the basis of higher EC categories, such as 1.1.1.- or 2.3.-.-, are ignored.
- The number of ModelSEED reactions aliased by EC numbers and KEGG reactions is reported. Some EC numbers, like 1.1.1.1, encompass a large number of reactions, while others are specific to a single reaction. KEGG reactions are typically specific to a single ModelSEED reaction. Reaction network compounds involved in ModelSEED reactions that are specific to EC number and KEGG reaction annotations are relatively trustworthy as metabolites that can be cycled by the organisms.
- The number of EC numbers and KEGG reactions that annotate KOs is reported. Some KOs, like K00128 (aldehyde dehydrogenase), encode a variety of reactions which are not necessarily catalyzed by the particular enzyme. Reaction network compounds that derive from KOs known to catalyze specific reactions are relatively trustworthy.

##### Production pathway

Matching compounds must be produced by reactions in a network, not just consumed. Furthermore, reactions are more likely to occur in the organism when they are connected to other reactions encoded by the network rather isolated from the network, particularly when the substrate and product do not arise from or feed into other reactions in the network. One way to assess reaction connectivity is by checking the connectivity of KOs encoding the reaction in KEGG pathway maps.

There can also be uncertainty in gene KO annotations. Sometimes a lower-ranking KO hit rather than the top hit represents the true protein. Erroneous KOs in the network can result in erroneous reactions and compounds. KOs are more likely to be accurate when assigned to multiple genes and when they co-occur with other KOs in KEGG pathway maps.

##### Compound chemistry

Chemical considerations support the existence of a matching compound. Predicted compounds are more likely to exist when they have a propensity to be retained in sample extraction and to be ionized given the mass spectrometric setup. In our study, SPE cartridges are more likely to retain hydrophobic compounds and negative ion mode is more likely to ionize compounds such as carboxylic and phenolic acids that can attain a -1 charge.

##### Known biological isomers

There is the possibility that the true compound represented by a formula is not encoded in the reaction network. It is therefore useful to compare the number of compounds with the formula in the network to the number in a large database of metabolites. We find the number of isomeric compounds in the ModelSEED Biochemistry compounds database. The ModelSEED database includes pesticides and other synthetic compounds, many of which are not represented in the KEGG compound database, one of the databases integrated into the ModelSEED database. Thus we also subset the isomeric ModelSEED compounds to those in the KEGG database. Furthermore, we count the number of these isomeric KEGG compounds that participate in KEGG reactions, as these tend to be more common biological substrates. All else equal, matching compounds with fewer isomers in the reference databases are more likely to actually be in the culture than those with more isomers.

## Reproducing this workflow


### Computational environment

This workflow uses the development version of [anvi'o](https://anvio.org) (`8-dev`), which you can install and activate following [anvi'o installation instructions](https://anvio.org/install/#development-version). Any more recent version of anvi'o should also work successfully. Load the anvi'o conda environment before running the workflow. The ModelSEED database should be installed in the default location for the anvi'o environment by {% include PROGRAM name="anvi-setup-modelseed-database" %}.

The computational demands of reproducing the workflow are minimal, all commands below should run within a few minutes or less on a modest laptop.


### The data pack

Below you will find brief descriptions of individual files used in our downstream analyses, if you would like to follow this workflow, you can download the following data pack that includes the four Roseobacter genomes and the metabolomics table associated with each culture experiment. For this, please open a terminal, create a work directory, and type the following commands (or replace directory names manually):

``` bash
# make sure there is a Downloads directory at your home
mkdir -p ~/Downloads

# change your current directory
cd ~/Downloads

# download the data pack
curl -o roseobacter-metabolomics.tar.gz https://merenlab.org/data/roseobacter-metabolomics/files/roseobacter-metabolomics.tar.gz

# unpack the data pack
tar -zxvf roseobacter-metabolomics.tar.gz

# go into the resulting data directory:
cd roseobacter-metabolomics
```

If you are here, you should be looking at a directory structure like this:

```
.
├── SH4-CONTIGS.db
├── SH40-CONTIGS.db
├── SH24-CONTIGS.db
├── SH22-CONTIGS.db
├── roseobacter-metabolomics-data.tsv
```

### Genomes

The files with the extension `.db` represent the four isolate genomes sequenced with PacBio Hifi long reads. To include them in our computational workflows we used the anvi'o program {% include PROGRAM name="anvi-gen-contigs-database" %} to turn the FASTA files into so-called {% include ARTIFACT name="contigs-db" %} files for downstream analyses. This file format contains much more information than a FASTA file, including gene coordinates, function annotations, and metabolic module membership of individual genes that will be essential to have in this workflow.

You can use the program {% include PROGRAM name="anvi-db-info" %} to learn more about the contents of a given {% include ARTIFACT name="contigs-db" %}:

``` bash
anvi-db-info SH22-CONTIGS.db

DB Info (no touch)
===============================================
Database Path ................................: SH22-CONTIGS.db
description ..................................: [Not found, but it's OK]
db_type ......................................: contigs (variant: unknown)
version ......................................: 24


DB Info (no touch also)
===============================================
project_name .................................: S_marinus_SH22
contigs_db_hash ..............................: hash52f2e51b
split_length .................................: 20000
kmer_size ....................................: 4
num_contigs ..................................: 4
total_length .................................: 4087537
num_splits ...................................: 201
gene_level_taxonomy_source ...................: None
genes_are_called .............................: 1
external_gene_calls ..........................: 0
external_gene_amino_acid_seqs ................: 0
skip_predict_frame ...........................: 0
splits_consider_gene_calls ...................: 1
trna_taxonomy_was_run ........................: 0
trna_taxonomy_database_version ...............: None
creation_date ................................: 1717747570.84562
modules_db_hash ..............................: a2b5bde358bb
scg_taxonomy_was_run .........................: 1
scg_taxonomy_database_version ................: GTDB: v214.1; Anvi'o: v1
gene_function_sources ........................: COG20_FUNCTION,Transfer_RNAs,CAZyme,KOfam,KEGG_Module,COG20_CATEGORY,KEGG_BRITE,KEGG_Class,COG20_PATHWAY
reaction_network_ko_annotations_hash .........: 1e5748cd73acfd2c24692de4d2c488044059aa32
reaction_network_kegg_database_release .......: 5a9644d40061
reaction_network_modelseed_database_sha ......: 194ac8afe48f8a606c0dd07ba3c7af10c02ba2fd

* Please remember that it is never a good idea to change these values. But in some
  cases it may be absolutely necessary to update something here, and a
  programmer may ask you to run this program and do it. But even then, you
  should be extremely careful.


AVAILABLE GENE CALLERS
===============================================
* 'prodigal' (3,851 gene calls)
* 'Transfer_RNAs' (45 gene calls)
* 'Ribosomal_RNA_23S' (3 gene calls)
* 'Ribosomal_RNA_16S' (3 gene calls)


AVAILABLE FUNCTIONAL ANNOTATION SOURCES
===============================================
* CAZyme (129 annotations)
* COG20_CATEGORY (3,208 annotations)
* COG20_FUNCTION (3,208 annotations)
* COG20_PATHWAY (830 annotations)
* KEGG_BRITE (2,396 annotations)
* KEGG_Class (511 annotations)
* KEGG_Module (511 annotations)
* KOfam (2,400 annotations)
* Transfer_RNAs (45 annotations)


AVAILABLE HMM SOURCES
===============================================
* 'Archaea_76' (76 models with 35 hits)
* 'Bacteria_71' (71 models with 72 hits)
* 'Protista_83' (83 models with 3 hits)
* 'Ribosomal_RNA_12S' (1 model with 0 hits)
* 'Ribosomal_RNA_16S' (3 models with 3 hits)
* 'Ribosomal_RNA_18S' (1 model with 0 hits)
* 'Ribosomal_RNA_23S' (2 models with 3 hits)
* 'Ribosomal_RNA_28S' (1 model with 0 hits)
* 'Ribosomal_RNA_5S' (5 models with 0 hits)
* 'Transfer_RNAs' (61 models with 45 hits)
```

Or get a standard FASTA file for a given genome using the program {% include PROGRAM name="anvi-export-contigs" %}:

``` bash
anvi-export-contigs -c SH22-CONTIGS.db -o SH22.fa
```

### Metabolomics table

The other file in this data pack, `roseobacter-metabolomics-data.tsv`, contains the processed spectral data, including monoisotopic molecular formulas and sample abundances. This is the same file that appears in our publication by Füssel et al. as the SI Table 2b.

This is how the first few lines of this table looks like, so you can browse the individual columns that are included:

|**`mz`**|**`diff`**|**`reference`**|**`formula`**|**`formula_isotopefree`**|**`formula_ion`**|**`homseries`**|**`totalc`**|**`HC`**|**`OC`**|**`C`**|**`H`**|**`O`**|**`N`**|**`S`**|**`P`**|**`MDL_3`**|**`ResPow`**|**`m1`**|**`SE`**|**`present_in`**|**`AI`**|**`AImod`**|**`DBE`**|**`Aromatic`**|**`AromaticO_rich`**|**`AromaticO_poor`**|**`Highlyunsaturated`**|**`HighlyunsaturatedO_rich`**|**`HighlyunsaturatedO_poor`**|**`Unsaturated`**|**`UnsaturatedO_rich`**|**`UnsaturatedO_poor`**|**`UnsaturatedwithN`**|**`Saturated`**|**`SaturatedO_rich`**|**`SaturatedO_poor`**|**`mean_signal_to_MDL`**|**`homnetworkmember`**|**`diff_filter`**|**`alternative_formula`**|**`SH4_Start`**|**`SH22_Start`**|**`SH24_Start`**|**`SH40_Start`**|**`SH22_SH4_Start`**|**`SH24_SH4_Start`**|**`SH4_SH40_Start`**|**`SH22_SH24_Start`**|**`SH22_SH40_Start`**|**`SH24_SH40_Start`**|**`SH22_SH24_SH4_Start`**|**`SH22_SH4_SH40_Start`**|**`SH24_SH4_SH40_Start`**|**`SH22_SH24_SH40_Start`**|**`SH22_SH24_SH4_SH40_Start`**|**`SH4_Final`**|**`SH22_Final`**|**`SH24_Final`**|**`SH40_Final`**|**`SH22_SH4_Final`**|**`SH24_SH4_Final`**|**`SH4_SH40_Final`**|**`SH22_SH24_Final`**|**`SH22_SH40_Final`**|**`SH24_SH40_Final`**|**`SH22_SH24_SH4_Final`**|**`SH22_SH4_SH40_Final`**|**`SH24_SH4_SH40_Final`**|**`SH22_SH24_SH40_Final`**|**`SH22_SH24_SH4_SH40_Final`**|
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
|95.0138555132338|0.0281917|96.021129238|C_5 H_4 O_2|C5H4O2|C_5 H_3 O_2|4610|5|0.800|0.400|5|4|2|0|0|0|2234458.755|2256845.978|95.0138533467|0.0475922619|46|0.67|0.75|4|1|0|1|0|0|0|0|0|0|0|0|0|0|2.40|50|FALSE|NA|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00005244|0.00006392|0.00040420|0.00000000|0.00000000|0.00000000|0.00000000|0.00007395|0.00000000|0.00000000|
|95.0502412690511|0.0320849|96.057514619|C_6 H_8 O_1|C6H8O|C_6 H_7 O_1|4610|6|1.333|0.167|6|8|1|0|0|0|2234458.755|2147563.946|95.0502396249|0.0406344107|56|0.40|0.45|3|0|0|0|1|0|1|0|0|0|0|0|0|0|2.36|50|FALSE|NA|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00005490|0.00000000|0.00000000|0.00000000|0.00000000|0.00010295|0.00010424|0.00000000|0.00000000|0.00000000|0.00015376|0.00000000|0.00000000|0.00000000|0.00000000|
|97.0294959619701|0.0698105|98.036779238|C_5 H_6 O_2|C5H6O2|C_5 H_5 O_2|4610|5|1.200|0.400|5|6|2|0|0|0|2235191.822|2138590.306|97.0294935409|0.0620774329|85|0.33|0.50|3|0|0|0|1|0|1|0|0|0|0|0|0|0|2.70|50|FALSE|NA|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00007827|0.00000000|0.00000000|0.00000000|0.00000000|0.00011741|0.00022261|0.00000000|0.00000000|0.00000000|0.00014084|0.00008886|0.00020904|0.00010036|0.00014362|
|100.0404053413720|0.0349439|101.047678242|C_4 H_7 O_2 N_1|C4H7NO2|C_4 H_6 O_2 N_1|3423|4|1.750|0.500|4|7|2|1|0|0|2236291.873|2138425.208|100.0404041292|0.0737641951|24|0.00|0.00|2|0|0|0|0|0|0|1|0|1|1|0|0|0|1.41|27|FALSE|NA|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00005162|0.00011517|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|
|101.0244171971140|0.0022346|102.031693857|C_4 H_6 O_3|C4H6O3|C_4 H_5 O_3|4610|4|1.500|0.750|4|6|3|0|0|0|2236658.677|2100978.308|101.0244159934|0.0508603586|156|0.00|0.20|2|0|0|0|0|0|0|1|1|0|0|0|0|0|1.82|50|FALSE|NA|0.00000000|0.00000000|0.00009419|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00007653|0.00000000|0.00000000|0.00000000|0.00000000|0.00011193|0.00012835|0.00012838|0.00006791|0.00000000|0.00009100|0.00004771|0.00011538|0.00000000|0.00007668|
|101.0396688839130|0.0459102|102.046950000|C_8 H_6|C8H6|C_8 H_5|4610|8|0.750|0.000|8|6|0|0|0|0|2236658.677|1976775.861|101.0396669676|0.0686485050|36|0.75|0.75|6|1|0|1|0|0|0|0|0|0|0|0|0|0|2.69|50|FALSE|NA|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00006977|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00022576|0.00000000|0.00000000|0.00000000|0.00000000|0.00012265|0.00000000|0.00000000|
|102.0560542629160|0.0237912|103.063328242|C_4 H_9 O_2 N_1|C4H9NO2|C_4 H_8 O_2 N_1|3423|4|2.250|0.500|4|9|2|1|0|0|2237025.541|2026054.176|102.0560530102|0.1558467278|17|0.00|0.00|1|0|0|0|0|0|0|0|0|0|0|0|0|0|1.53|27|FALSE|NA|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00010476|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|
|102.9859214249310|0.0325117|103.993201238|C_3 H_4 O_2 S_1|C3H4O2S|C_3 H_3 O_2 S_1|3715|3|1.333|0.667|3|4|2|0|1|0|2237392.465|1957173.000|102.9859199500|0.2308225502|27|0.00|0.00|2|0|0|0|1|1|0|0|0|0|0|0|0|0|1.51|30|FALSE|NA|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00006466|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|
|103.0400590982740|0.0800309|104.047343857|C_4 H_8 O_3|C4H8O3|C_4 H_7 O_3|4610|4|2.000|0.750|4|8|3|0|0|0|2237392.465|1961126.500|103.0400536960|0.0554939082|100|0.00|0.00|1|0|0|0|0|0|0|1|1|0|0|0|0|0|79.08|50|FALSE|NA|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.02121257|0.00000000|0.00039499|0.00000000|0.00000000|0.00185427|0.00250969|0.00037960|0.00006953|0.00011791|0.00027330|0.00017926|0.00143951|0.00048814|0.00040238|
|103.0553130612920|0.1009681|104.062600000|C_8 H_8|C8H8|C_8 H_7|4610|8|1.000|0.000|8|8|0|0|0|0|2237392.465|2085497.409|103.0553111458|0.0903881368|22|0.62|0.62|5|1|0|1|0|0|0|0|0|0|0|0|0|0|1.56|50|FALSE|NA|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00004035|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|
|104.0505372294390|0.3364148|105.057849004|C_7 H_7 N_1|C7H7N|C_7 H_6 N_1|3423|7|1.000|0.000|7|7|0|1|0|0|2237759.450|1930600.560|104.0505355175|0.5149217715|25|0.67|0.67|5|1|0|1|0|0|0|0|0|0|0|0|0|0|1.39|27|FALSE|NA|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00008673|0.00000000|0.00000000|0.00000000|0.00000000|
|105.0345944057540|0.0586372|106.041864619|C_7 H_6 O_1|C7H6O|C_7 H_5 O_1|4610|7|0.857|0.143|7|6|1|0|0|0|2238126.495|2008009.077|105.0345925729|0.3127859784|39|0.67|0.69|5|1|0|1|0|0|0|0|0|0|0|0|0|0|1.40|50|FALSE|NA|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00006203|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00005143|0.00000000|
|107.0502357662910|0.0224047|108.057514619|C_7 H_8 O_1|C7H8O|C_7 H_7 O_1|4610|7|1.143|0.143|7|8|1|0|0|0|2238860.765|1887046.786|107.0502337204|0.0583066930|126|0.50|0.54|4|1|0|1|0|0|0|0|0|0|0|0|0|0|7.52|50|FALSE|NA|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00026394|0.00000000|0.00000000|0.00000000|0.00019360|0.00067382|0.00078355|0.00070539|0.00000000|0.00011044|0.00061324|0.00024031|0.00106589|0.00070622|0.00057724|
|108.0454820492240|0.0471515|109.052763623|C_6 H_7 O_1 N_1|C6H7NO|C_6 H_6 O_1 N_1|3423|6|1.167|0.167|6|7|1|1|0|0|2239227.990|1887069.400|108.0454815707|0.1086665235|60|0.50|0.56|4|1|0|1|0|0|0|0|0|0|0|0|0|0|3.78|27|FALSE|NA|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00053209|0.00030703|0.00000000|0.00000000|0.00006015|0.00000000|0.00027031|0.00038705|0.00000000|0.00008938|0.00000000|0.00000000|0.00000000|
|109.0658885686900|0.0034613|110.073164619|C_7 H_10 O_1|C7H10O|C_7 H_9 O_1|4610|7|1.429|0.143|7|10|1|0|0|0|2239595.276|1907711.752|109.0658874467|0.0511481349|137|0.33|0.38|3|0|0|0|1|0|1|0|0|0|0|0|0|0|2.82|50|FALSE|NA|0.00005259|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00008990|0.00000000|0.00000000|0.00000000|0.00000000|0.00006936|0.00000000|0.00008719|0.00010021|0.00019081|0.00000000|0.00000000|0.00000000|0.00000000|0.00037872|0.00033969|0.00000000|0.00000000|0.00000000|0.00012937|0.00000000|0.00000000|0.00000000|0.00000000|
|111.0087713970950|0.0354592|112.016043857|C_5 H_4 O_3|C5H4O3|C_5 H_3 O_3|4610|5|0.800|0.600|5|4|3|0|0|0|2240330.028|1860703.739|111.0087700372|0.0676582216|119|0.50|0.71|4|1|1|0|0|0|0|0|0|0|0|0|0|0|1.53|50|FALSE|NA|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00003877|0.00000000|0.00008944|0.00000000|0.00000000|0.00006990|0.00009919|0.00007880|0.00000000|0.00000000|0.00000000|0.00005030|0.00011178|0.00006070|0.00006852|
|114.0560473631960|0.0386483|115.063328242|C_5 H_9 O_2 N_1|C5H9NO2|C_5 H_8 O_2 N_1|3423|5|1.800|0.400|5|9|2|1|0|0|2241432.608|1828463.053|114.0560475526|0.1283682584|95|0.00|0.00|2|0|0|0|0|0|0|1|0|1|1|0|0|0|1.86|27|FALSE|NA|0.00000000|0.00000000|0.00000000|0.00008120|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00012011|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00022755|0.00021206|0.00000000|0.00000000|0.00005687|0.00000000|0.00009448|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|
|115.0036843829750|0.0201619|116.010958476|C_4 H_4 O_4|C4H4O4|C_4 H_3 O_4|4610|4|1.000|1.000|4|4|4|0|0|0|2241800.255|1798563.500|115.0036840529|0.0779625901|94|0.00|0.50|3|0|0|0|1|1|0|0|0|0|0|0|0|0|1.66|50|FALSE|NA|0.00005517|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00006655|0.00000000|0.00000000|0.00000000|0.00007120|0.00005576|0.00000000|
|115.0553195669310|0.0344728|116.062600000|C_9 H_8|C9H8|C_9 H_7|4610|9|0.889|0.000|9|8|0|0|0|0|2241800.255|1757533.933|115.0553173348|0.1109887279|45|0.67|0.67|6|1|0|1|0|0|0|0|0|0|0|0|0|0|3.44|50|FALSE|NA|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00010459|0.00000000|0.00000000|0.00000000|0.00000000|0.00006484|0.00010215|0.00031981|0.00000000|0.00000000|0.00000000|0.00000000|0.00021672|0.00011680|0.00000000|
|117.0345855886080|0.0220176|118.041864619|C_8 H_6 O_1|C8H6O|C_8 H_5 O_1|4610|8|0.750|0.125|8|6|1|0|0|0|2242535.730|1693286.109|117.0345837876|0.1468502474|64|0.71|0.73|6|1|0|1|0|0|0|0|0|0|0|0|0|0|7.57|50|FALSE|NA|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00025536|0.00000000|0.00000000|0.00000000|0.00000000|0.00016141|0.00022868|0.00108416|0.00000000|0.00000000|0.00000000|0.00005057|0.00056176|0.00048878|0.00009414|
|117.0557115866400|0.0494567|118.062993857|C_5 H_10 O_3|C5H10O3|C_5 H_9 O_3|4610|5|2.000|0.600|5|10|3|0|0|0|2242535.730|1735420.716|117.0557092044|0.2877162635|74|0.00|0.00|1|0|0|0|0|0|0|1|1|0|0|0|0|0|2.49|50|FALSE|NA|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00013221|0.00000000|0.00000000|0.00000000|0.00000000|0.00011136|0.00014889|0.00010388|0.00000000|0.00000000|0.00009287|0.00006809|0.00015296|0.00010295|0.00010574|
|117.0709629508030|0.0899234|118.078250000|C_9 H_10|C9H10|C_9 H_9|4610|9|1.111|0.000|9|10|0|0|0|0|2242535.730|1643075.394|117.0709600004|0.2995694356|71|0.56|0.56|5|1|0|1|0|0|0|0|0|0|0|0|0|0|18.74|50|FALSE|NA|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00163395|0.00000000|0.00000000|0.00000000|0.00000000|0.00068428|0.00090488|0.00193652|0.00000000|0.00000000|0.00042827|0.00013660|0.00126878|0.00045145|0.00030455|
|118.0298358415610|0.0113410|119.037113623|C_7 H_5 O_1 N_1|C7H5NO|C_7 H_4 O_1 N_1|3423|7|0.714|0.143|7|5|1|1|0|0|2242903.558|1786864.132|118.0298356003|0.0621286083|38|0.80|0.82|6|1|0|1|0|0|0|0|0|0|0|0|0|0|2.83|27|FALSE|NA|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00059625|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|
|118.0331989036020|0.0864244|119.040485623|C_4 H_9 O_1 N_1 S_1|C4H9NOS|C_4 H_8 O_1 N_1 S_1|2328|4|2.250|0.250|4|9|1|1|1|0|2242903.558|1726953.242|118.0331978601|0.1764849790|33|0.00|0.00|1|0|0|0|0|0|0|0|0|0|0|0|0|0|1.75|8|FALSE|NA|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00009807|0.00000000|0.00000000|0.00000000|0.00000000|0.00005414|0.00010486|0.00000000|0.00000000|
|118.0509610909830|0.0448352|119.058242861|C_4 H_9 O_3 N_1|C4H9NO3|C_4 H_8 O_3 N_1|3423|4|2.250|0.750|4|9|3|1|0|0|2242903.558|1807663.056|118.0509604080|0.0559783664|36|0.00|0.00|1|0|0|0|0|0|0|0|0|0|0|0|0|0|3.16|27|FALSE|NA|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00055533|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00012309|0.00014975|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|
|118.0662130626320|0.0798666|119.073499004|C_8 H_9 N_1|C8H9N|C_8 H_8 N_1|3423|8|1.125|0.000|8|9|0|1|0|0|2242903.558|1749825.800|118.0662119608|0.0620450414|55|0.57|0.57|5|1|0|1|0|0|0|0|0|0|0|0|0|0|2.52|27|FALSE|NA|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00002750|0.00000000|0.00053495|0.00000000|0.00000000|0.00000000|0.00006186|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00009519|0.00008829|0.00000000|
|119.0138366416140|0.1346846|120.021129238|C_7 H_4 O_2|C7H4O2|C_7 H_3 O_2|4610|7|0.571|0.286|7|4|2|0|0|0|2243271.447|1682790.200|119.0138336946|0.1518552092|45|0.80|0.83|6|1|0|1|0|0|0|0|0|0|0|0|0|0|2.05|50|FALSE|NA|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00005108|0.00000000|0.00000000|0.00000000|0.00000000|0.00007889|0.00007796|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|
|119.0172143468010|0.0871489|120.024501238|C_4 H_8 O_2 S_1|C4H8O2S|C_4 H_7 O_2 S_1|3715|4|2.000|0.500|4|8|2|0|1|0|2243271.447|1712802.154|119.0172119940|0.0658483958|26|0.00|0.00|1|0|0|0|0|0|0|1|0|1|0|0|0|0|1.62|30|FALSE|NA|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00006240|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00006782|
|119.0349869692270|0.0410272|120.042258476|C_4 H_8 O_4|C4H8O4|C_4 H_7 O_4|4610|4|2.000|1.000|4|8|4|0|0|0|2243271.447|1748584.976|119.0349868693|0.1931883383|42|0.00|0.00|1|0|0|0|0|0|0|1|1|0|0|0|0|0|1.69|50|FALSE|NA|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00012461|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|0.00011364|0.00012481|0.00000000|0.00000000|0.00000000|0.00000000|0.00000000|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|

## Creating the reaction networks

The workflow was originally implemented in a Jupyter notebook that used anvi'o libraries that processed the anvi'o {% include ARTIFACT name="reaction-network" %}

Python cells of the Jupyter notebook are split up into sections of this workflow document with accompanying explanations and analyses of the output. If you wish to follow along interactively in Python, you can enter `python3` in your terminal and run the code blocks sequentially. The following package imports are first required in the Python interactive shell.

```python
import os
import sys
import rdkit
import itertools
import numpy as np
import pandas as pd

from rdkit import Chem
from copy import deepcopy
from typing import Iterable
from collections import defaultdict
```

Imports from anvi'o will fail if the anvi'o package isn't in the Python module search path, a problem that can arise in Jupyter notebooks. If you have used the standard installation instructions on the anvi'o installation page, running this command should solve that issue:

```python
sys.path.append('~/github/anvio')
```

Then you should be able to run these two lines without any errors:

```python
import anvio
import anvio.reactionnetwork as rn
```


### Genomic networks

The four {% include ARTIFACT name="contigs-db" %} files in our data pack for the four strains we have worked with contain gene calls with KO annotations and reaction networks based on the KOs. Genes were annotated with KOs using {% include PROGRAM name="anvi-run-kegg-kofams" %}, and networks were constructed with {% include PROGRAM name="anvi-reaction-network" %}.

List the strains and their database files. Load reaction networks into memory. The dictionary of reaction networks is keyed by a tuple, as co-culture "metagenomic" reaction networks keyed by tuples of strain IDs will be added to the dictionary.

```python
all_strains = ['SH22', 'SH24', 'SH4', 'SH40']
strain_names = {
    'SH22': 'Sulfitobacter sp. SH22-1',
    'SH24': 'Sulfitobacter sp. SH24-1b',
    'SH4': 'Pelagimonas varians SH4-1',
    'SH40': 'Phaeobacter sp. SH40'
}
all_contigs_dbs = [f'{strain}-CONTIGS.db' for strain in all_strains]
```

```python
con = rn.Constructor()
all_networks: dict[tuple[str], rn.GenomicNetwork] = {}
for contigs_db in all_contigs_dbs:
    strain = contigs_db[: contigs_db.index('-CONTIGS.db')]
    all_networks[(strain, )] = con.load_contigs_database_network(contigs_db, quiet=True)
```

### Remove EC categories from networks

Avoid the inclusion of reactions on the basis of higher EC categories, such as 1.1.1.- or 2.3.-.-, that annotate KOs. Higher categories encompass a range of ModelSEED reactions that cannot be confidently attributed to the particular enzyme. Inclusion of these reactions increases the likelihood of false positive formula matches to compounds that are not actually produced by the organism. Networks filtered to removed EC categories are called "refined" networks. These network are used in formula matching.

```python
all_refined_networks: dict[tuple[str], rn.GenomicNetwork] = {}
for strain_combo, unrefined_network in all_networks.items():
    modelseed_reaction_ids_to_retain = []
    for ko in unrefined_network.kos.values():
        modelseed_reaction_ids_to_check = []
        for modelseed_reaction_id, ec_numbers in ko.ec_number_aliases.items():
            for ec_number in ec_numbers:
                if '-' not in ec_number:
                    modelseed_reaction_ids_to_retain.append(modelseed_reaction_id)
                    break
            else:
                modelseed_reaction_ids_to_check.append(modelseed_reaction_id)
        for modelseed_reaction_id in modelseed_reaction_ids_to_check:
            if modelseed_reaction_id in ko.kegg_reaction_aliases:
                modelseed_reaction_ids_to_retain.append(modelseed_reaction_id)
    modelseed_reaction_ids_to_retain = set(modelseed_reaction_ids_to_retain)
    refined_network = unrefined_network.subset_network(reactions_to_subset=modelseed_reaction_ids_to_retain)
    all_refined_networks[strain_combo] = refined_network
```

### Networks based on KEGG reactions

Compare the sizes of reaction networks constructed in two ways: first, using the default of both KEGG reaction and EC number annotations of KOs, and second, using just KEGG reaction annotations. KEGG reactions are more specific than EC numbers, which often map to a larger group of reactions in the ModelSEED database, as explained above in [Specificity of reaction annotations](#specificity-of-reaction-annotations). The "EC+KEGG" network is prone to more false positive formula-compound matches that must be evaluated and fewer false negatives, or missing formula-compound matches, than the "just KEGG" network.

Although it would be useful to design a flag in {% include PROGRAM name="anvi-reaction-network" %} that allows a network to be constructed from KEGG reactions excluding EC numbers, for now we will remove the parts of the "EC+KEGG" networks that are based solely on EC numbers. This is achieved using the function that subsets networks by select items.

```python
all_kegg_networks: dict[tuple[str], rn.GenomicNetwork] = {}
for strain_combo, ec_kegg_network in all_networks.items():
    modelseed_reaction_ids_to_retain = []
    for ko in ec_kegg_network.kos.values():
        for modelseed_reaction_id in ko.kegg_reaction_aliases:
            modelseed_reaction_ids_to_retain.append(modelseed_reaction_id)
    modelseed_reaction_ids_to_retain = set(modelseed_reaction_ids_to_retain)
    kegg_network = ec_kegg_network.subset_network(reactions_to_subset=modelseed_reaction_ids_to_retain)
    all_kegg_networks[strain_combo] = kegg_network
```

### Co-culture "metagenomic" networks

Merge genomic reaction networks to represent co-culture "metagenomic" reaction networks. The network merge function avoids duplicate entries, such as KOs or reactions shared by both networks. Genes with identical anvi'o gene caller IDs (GCIDs) in different genomes would be considered the same in merging, so the identity of the genes must be maintained by adjusting integer GCIDs to be non-overlapping. Since the number of genes in these genomes is less than 10,000, add 10,000 to SH22 genome GCIDs, 20,000 to SH24 GCIDs, 30,000 to SH4 GCIDs, and 40,000 to SH40 GCIDs. Each gene in the network can thereby be traced back to the source genome, with SH22 genes, for example, have GCIDs between 10,000 and 20,000.

```python
def make_gcids_nonoverlapping(networks: dict[tuple[str], rn.GenomicNetwork], increment: int = 10000) -> None:
    i = increment
    for network in networks.values():
        gcids_to_remove = []
        for gcid, gene in network.genes.items():
            assert gcid < increment
            new_gcid = i + gcid
            gene.gcid = new_gcid
            gcids_to_remove.append(gcid)
        for gcid in gcids_to_remove:
            gene = network.genes.pop(gcid)
            network.genes[gene.gcid] = gene
        i += increment
```

```python
def merge_networks(networks: dict[tuple[str], rn.GenomicNetwork]) -> None:
    merged_networks = {}
    for r in range(2, len(networks) + 1):
        for combo in itertools.combinations(networks.items(), r):
            merged_strains = tuple()
            merged_network = None
            for strains, network in combo:
                merged_strains += strains
                if merged_network is None:
                    merged_network = network
                else:
                    merged_network = merged_network.merge_network(network)
            merged_networks[merged_strains] = merged_network
    networks.update(merged_networks)
```

```python
make_gcids_nonoverlapping(all_networks)
merge_networks(all_networks)

make_gcids_nonoverlapping(all_refined_networks)
merge_networks(all_refined_networks)

make_gcids_nonoverlapping(all_kegg_networks)
merge_networks(all_kegg_networks)
```

List the strain combination tuples identifying the co-culture networks.

```python
all_strain_combos = list(all_networks)
```

### Compare networks constructed with different KO annotations

Compare the three types of networks constructed on the basis of varying KO annotations: KEGG reactions and all EC numbers (default networks), KEGG reactions and EC numbers but not higher EC categories ("refined networks"), and just KEGG reactions ("KEGG networks"). How many compounds are removed from the default networks excluding higher EC categories and EC numbers altogether?

```python
header = ['strains', 'EC+KEGG_network_compounds', 'refined_network_compounds', 'KEGG_network_compounds']
rows = []
for strain_combo, ec_kegg_network in all_networks.items():
    refined_network = all_refined_networks[strain_combo]
    kegg_network = all_kegg_networks[strain_combo]
    row = []
    row.append('_'.join(strain_combo))
    row.append(len(ec_kegg_network.metabolites))
    row.append(len(refined_network.metabolites))
    row.append(len(kegg_network.metabolites))
    rows.append(row)
network_compound_counts = pd.DataFrame(rows, columns=header).set_index('strains')
network_compound_counts['refined_compound_fraction'] = network_compound_counts['refined_network_compounds'] / network_compound_counts['EC+KEGG_network_compounds']
network_compound_counts['KEGG_compound_fraction'] = network_compound_counts['KEGG_network_compounds'] / network_compound_counts['EC+KEGG_network_compounds']
print(network_compound_counts.to_string())
```

```python
mean_refined_compound_fraction = network_compound_counts['refined_compound_fraction'].mean()
mean_kegg_compound_fraction = network_compound_counts['KEGG_compound_fraction'].mean()
print(f"An average of {round((1 - mean_refined_compound_fraction) * 100, 1)}% of compounds in the \"EC+KEGG\" network are removed ignoring higher EC categories in the \"refined\" network")
print(f"{round((1 - mean_kegg_compound_fraction) * 100, 1)}% of compounds in the \"EC+KEGG\" network are removed ignoring EC numbers and only considering KEGG reactions in the \"KEGG\" network")
```

On average 40.8% of compounds in the "EC+KEGG" network are removed ignoring higher EC categories in the "refined" network. On average 73.2% of compounds in the "EC+KEGG" network are removed ignoring EC numbers and only considering KEGG reactions in the "KEGG" network.

## Prepare metabolomics data

Load the metabolomics data table, SI Table 2b from the paper. Each row represents a monoisotopic molecular feature.

```python
roseobacter_metabolomics_df = pd.read_csv('roseobacter-metabolomics-data.tsv', sep='\t', header=0)
```

Confirm that a unique molecular formula was assigned to each feature.

```python
len(roseobacter_metabolomics_df) == roseobacter_metabolomics_df['formula_isotopefree'].nunique()
```

### Molecular feature production and consumption

Most molecular features in the table were not initially present in the culture and produced over the course of the experiment.

```python
culture_feature_change_stats: dict[tuple[str], dict[str, int]] = {}
sum_consumed_fraction = 0
sum_produced_fraction = 0
sum_produced_new_fraction = 0
for strain_combo in all_strain_combos:
    culture_feature_change_stats[strain_combo] = feature_change_stats = {}
    feature_change_stats['changed'] = len(roseobacter_metabolomics_df[roseobacter_metabolomics_df[f"{'_'.join(strain_combo)}_Final"] - roseobacter_metabolomics_df[f"{'_'.join(strain_combo)}_Start"] != 0])
    feature_change_stats['consumed'] = len(roseobacter_metabolomics_df[roseobacter_metabolomics_df[f"{'_'.join(strain_combo)}_Final"] - roseobacter_metabolomics_df[f"{'_'.join(strain_combo)}_Start"] < 0])
    feature_change_stats['produced'] = len(roseobacter_metabolomics_df[roseobacter_metabolomics_df[f"{'_'.join(strain_combo)}_Final"] - roseobacter_metabolomics_df[f"{'_'.join(strain_combo)}_Start"] > 0])
    feature_change_stats['produced_new'] = len(roseobacter_metabolomics_df[(roseobacter_metabolomics_df[f"{'_'.join(strain_combo)}_Final"] > 0) & (roseobacter_metabolomics_df[f"{'_'.join(strain_combo)}_Start"] == 0)])

    sum_consumed_fraction += feature_change_stats['consumed'] / feature_change_stats['changed']
    sum_produced_fraction += feature_change_stats['produced'] / feature_change_stats['changed']
    sum_produced_new_fraction += feature_change_stats['produced_new'] / feature_change_stats['changed']

mean_consumed_fraction = sum_consumed_fraction / len(all_strain_combos)
mean_produced_fraction = sum_produced_fraction / len(all_strain_combos)
mean_produced_new_fraction = sum_produced_new_fraction / len(all_strain_combos)

print(f"On average {round(mean_consumed_fraction * 100, 1)}% of compounds measured per culture were consumed (initially present and decreased in abundance)")
print(f"On average {round(mean_produced_fraction * 100, 1)}% of compounds measured per culture were produced (increased in abundance)")
print(f"On average {round(mean_produced_new_fraction * 100, 1)}% of compounds measured per culture were newly produced (not initially present)")
```

On average:
- 8.6% of features were consumed (initially present and decreased in abundance)
- 91.4% of features were produced (increased in abundance)
- 89.7% of features were newly produced (not initially present)

### Add deprotonated formulas

Add formulas for deprotonated versions of compounds as they may exist in the aqueous solution of cultures and the ModelSEED database used to populate compounds in reaction networks. Allow up to 2 hydrogens, 1 per oxygen, to be removed from each neutral formula. It does not make sense to remove 3 hydrogens in searching for common metabolites, since there are few with a -3 charge -- primarily the tricarboxylic acids citrate, isocitrate, and aconitate in the TCA cycle.

```python
new_rows = []
new_idx = 0
for _, row in feature_table.iterrows():
    new_row = row.drop(['formula_isotopefree_minus_1_H', 'formula_isotopefree_minus_2_H'])
    new_row.name = new_idx
    new_row['search_formula'] = row['formula_isotopefree']
    new_row['search_charge'] = 0
    new_rows.append(new_row)
    new_idx += 1
    break

    if not row['formula_isotopefree_minus_1_H']:
        continue
    new_row = row.drop(['formula_isotopefree_minus_1_H', 'formula_isotopefree_minus_2_H'])
    new_row.name = new_idx
    new_row['search_formula'] = row['formula_isotopefree_minus_1_H']
    new_row['search_charge'] = -1
    new_rows.append(new_row)
    new_idx += 1

    if not row['formula_isotopefree_minus_2_H']:
        continue
    new_row = row.drop(['formula_isotopefree_minus_1_H', 'formula_isotopefree_minus_2_H'])
    new_row.name = new_idx
    new_row['search_formula'] = row['formula_isotopefree_minus_2_H']
    new_row['search_charge'] = -2
    new_rows.append(new_row)
    new_idx += 1

feature_table = pd.DataFrame(new_rows)
last_col_names = ['search_formula', 'search_charge']
first_col_names = feature_table.columns.tolist()[: -2]
feature_table = feature_table[last_col_names + first_col_names]
```

Make a new version of the table with a row per formula protonation state.

```python
new_rows = []
new_idx = 0
for _, row in feature_table.iterrows():
    new_row = row.drop(['formula_isotopefree_minus_1_H', 'formula_isotopefree_minus_2_H'])
    new_row.name = new_idx
    new_row['search_formula'] = row['formula_isotopefree']
    new_row['search_charge'] = 0
    new_rows.append(new_row)
    new_idx += 1
    break

    if not row['formula_isotopefree_minus_1_H']:
        continue
    new_row = row.drop(['formula_isotopefree_minus_1_H', 'formula_isotopefree_minus_2_H'])
    new_row.name = new_idx
    new_row['search_formula'] = row['formula_isotopefree_minus_1_H']
    new_row['search_charge'] = -1
    new_rows.append(new_row)
    new_idx += 1

    if not row['formula_isotopefree_minus_2_H']:
        continue
    new_row = row.drop(['formula_isotopefree_minus_1_H', 'formula_isotopefree_minus_2_H'])
    new_row.name = new_idx
    new_row['search_formula'] = row['formula_isotopefree_minus_2_H']
    new_row['search_charge'] = -2
    new_rows.append(new_row)
    new_idx += 1

feature_table = pd.DataFrame(new_rows)
last_col_names = ['search_formula', 'search_charge']
first_col_names = feature_table.columns.tolist()[: -2]
feature_table = feature_table[last_col_names + first_col_names]
```

### Find database isomers

Find compounds in the ModelSEED Biochemistry database with molecular formulas, including deprotonated formulas. To help evaluate the number of possible biomolecular isomers that could exist as part of controlling false positive compound matches (see the section, [Known biological isomers](#known-biological-isomers)), subset isomeric compounds in the KEGG compound database, and those that participate in KEGG reactions.

```python
# Keys are (<formula>, <charge>), values are {<source of isomers>: [(<ModelSEED compound ID>, <ModelSEED compound name>)]}.
compound_isomers: dict[tuple[str, int], dict[str, list[tuple[str, str]]]] = {}
# Load the ModelSEED database from the default anvi'o installation location.
modelseed_db = rn.ModelSEEDDatabase()
compounds_table = modelseed_db.compounds_table

# Subset compounds with KEGG aliases.
compounds_with_kegg_alias_table = compounds_table[compounds_table['KEGG'].notna()]

# Subset compounds that participate in KEGG reactions.
kegg_reactions_table = modelseed_db.kegg_reactions_table
kegg_reaction_compound_ids = []
for compound_ids in kegg_reactions_table['compound_ids']:
    if not isinstance(compound_ids, str):
        continue
    compound_ids: str
    if compound_ids.strip() == '':
        continue
    for compound_id in compound_ids.split(';'):
        kegg_reaction_compound_ids.append(compound_id)
kegg_reaction_compound_ids = sorted(set(kegg_reaction_compound_ids))
select_rows = []
for row in compounds_with_kegg_alias_table.itertuples():
    if row.Index in kegg_reaction_compound_ids:
        select_rows.append(row)
compounds_with_kegg_reaction_table = pd.DataFrame(select_rows).set_index('Index')

for feature_row in feature_table.itertuples():
    formula = feature_row.search_formula
    charge = feature_row.search_charge
    compound_isomers[(formula, charge)] = isomers = {'modelseed_isomers': [], 'kegg_isomers': [], 'kegg_isomers_with_reaction': []}
    for compound_row in compounds_table[(compounds_table['formula'] == formula) & (compounds_table['charge'] == charge)].itertuples():
        isomers['modelseed_isomers'].append((compound_row.Index, compound_row.name))
    for compound_row in compounds_with_kegg_alias_table[(compounds_with_kegg_alias_table['formula'] == formula) & (compounds_with_kegg_alias_table['charge'] == charge)].itertuples():
        isomers['kegg_isomers'].append((compound_row.Index, compound_row.name))
    for compound_row in compounds_with_kegg_reaction_table[(compounds_with_kegg_reaction_table['formula'] == formula) & (compounds_with_kegg_reaction_table['charge'] == charge)].itertuples():
        isomers['kegg_isomers_with_reaction'].append((row.Index, row.name))
```

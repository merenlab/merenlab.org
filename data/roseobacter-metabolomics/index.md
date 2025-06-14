---
layout: page
title: A reproducible workflow for Füssel et al. 2025
modified: 2025-06-09
excerpt: "A bioinformatics workflow for genomically constrained compound prediction in Roseobacter co-culture experiments"
comments: true
authors: [sam]
---

<div class="extra-info" markdown="1">

<span class="extra-info-header">Summary</span>

**The purpose of this page** is to provide access to our bioinformatics workflow that predicted compound identifications of molecular features in our study titled "**Bacterial interactions shape the molecular composition of dissolved organic matter**" by Füssel et al.

In addition to providing transparency in our methods, this workflow can be used as the basis for genomically constrained compound prediction in other metabolomics experiments, which will also help refine and validate the approach.

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

## Downloading data for this workflow

Files needed for this workflow can be downloaded from the following sources.

### Genomes

The four isolate genomes were sequenced with PacBio Hifi long reads.

<!-- TODO: Submit SRAs. List SRAs. anvi-run-workflow command to download SRAs. samples.txt file listing fastq files. -->

### Metabolomics table

Processed spectral data, including monoisotopic molecular formulas and sample abundances, are contained in Füssel et al. SI Table 2b.

<!-- TODO: Upload tsv version of table, perhaps to Zenodo, and provide download link. -->

## Computational environment

In this workflow we used anvi'o version `8-dev`, which is the development version following the `v8` stable release. Any more recent version of anvi'o should also work successfully. Load the anvi'o conda environment before running the workflow. The ModelSEED database should be installed in the default location for the anvi'o environment by {% include PROGRAM name="anvi-setup-modelseed-database" %}.

Choose a working directory and store it as variable `$WD` in your terminal.

```bash
cd /where/you/want/to/work
WD=$PWD
```

The rest of this document will make use of variable `$WD`. The downloaded contigs databases and metabolomics table should be stored in this directory.

```
.
├── SH4-CONTIGS.db
├── SH40-CONTIGS.db
├── SH24-CONTIGS.db
├── SH22-CONTIGS.db
├── si_table_2b.tsv
```

The computational demands of the workflow are minimal, all commands running within a few minutes on a laptop.

### Python code

The workflow was originally implemented in a Jupyter notebook, which can be downloaded.

<!-- TODO: Upload Jupyter notebook to Zenodo and provide download link. -->

A Python script with the same commands can also be downloaded and run with `python3 compound_matching.py`.

<!-- TODO: Upload Python script to Zenodo and provide download link. -->

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

Imports from anvi'o will fail if the anvi'o package isn't in the Python module search path, a problem that can arise in Jupyter notebooks. Here is how it should be added, substituting the path to your anvi'o installation.

```python
sys.path.append('/path/to/top/directory/of/anvio')
```

```python
import anvio
import anvio.reactionnetwork as rn
```

## Create reaction networks

### Genomic networks

The downloaded contigs databases for the four strains contain gene calls with KO annotations and reaction networks based on the KOs. Genes were annotated with KOs using {% include PROGRAM name="anvi-run-kegg-kofams" %}, and networks were constructed with {% include PROGRAM name="anvi-reaction-network" %}.

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

Load the metabolomics data table, SI Table 2b from the paper. Each row represents a monoisotopic molecular features.

```python
si_table_2b = pd.read_csv('si_table_2b.tsv', sep='\t', header=0)
```

Confirm that a unique molecular formula was assigned to each feature.

```python
len(si_table_2b) == si_table_2b['formula_isotopefree'].nunique()
```

### Add deprotonated formulas

Add formulas for deprotonated versions of compounds as they may exist in the aqueous solution of cultures and the ModelSEED database used to populate compounds in reaction networks. Allow up to 2 hydrogens, 1 per oxygen, to be removed from each neutral formula. It does not make sense to remove 3 hydrogens in searching for common metabolites, since there are few with a -3 charge -- primarily the tricarboxylic acids citrate, isocitrate, and aconitate in the TCA cycle.

```python
min_protons_subtracted = 1
max_protons_subtracted = 2
formula_data = si_table_2b[['formula', 'formula_isotopefree', 'O', 'H']]

deprot_rows = []
for _, row in formula_data.iterrows():
    formula_isotopefree = row.formula_isotopefree

    atom_count = {}
    for atomic_entry in row.formula.split():
        atom, count = atomic_entry.split('_')
        count = int(count)
        atom_count[atom] = count

    deprot_row = []
    for num_protons_subtracted in range(min_protons_subtracted, max_protons_subtracted + 1):
        if num_protons_subtracted > row.O:
            deprot_row.append('')
            continue

        new_atom_count = atom_count.copy()
        new_atom_count['H'] = atom_count['H'] - num_protons_subtracted

        new_formula_isotopefree = ''
        for atom, count in new_atom_count.items():
            new_formula_isotopefree += f'{atom}{count}' if count > 1 else atom
        deprot_row.append(new_formula_isotopefree)
    deprot_rows.append(deprot_row)

header = [f'formula_isotopefree_minus_{num_protons_subtracted}_H' for num_protons_subtracted in range(min_protons_subtracted, max_protons_subtracted + 1)]
deprot_table = pd.DataFrame(deprot_rows, columns=header)

insert_after_col = 'formula_isotopefree'
before = si_table_2b.loc[:, :insert_after_col]
after = si_table_2b.loc[:, si_table_2b.columns > insert_after_col]
feature_table = pd.concat([before, deprot_table, after], axis=1)
```

### Add database isomers



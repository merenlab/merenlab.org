---
layout: page
title: Chapter IV - Reproducing Kiefl et al, 2022
modified: 2021-10-21
excerpt: "A complete reproducible workflow of the manuscript 'Structure-informed microbial population genetics elucidate selective pressures that shape protein evolution' by Kiefl et al"
comments: true
authors: [evan]
redirect_from:
  - 
---

{% capture images %}{{site.url}}/data/anvio-structure/images{% endcapture %}
{% capture command_style %}background: #D7484822; border: 4px solid #D74848;{% endcapture %}
{% capture analysis_style %}background: #E6DBE4{% endcapture %}

{:.warning}
This document is **UNDER CONSTRUCTION**. It is not in a state where you can yet reproduce our work. We anticipate this workflow will be finalized by late March, and will remove this message when it is complete.

## Quick Navigation

- [Chapter I: The prologue]({{ site.url }}/data/anvio-structure/chapter-I)
- [Chapter II: Configure your system]({{ site.url }}/data/anvio-structure/chapter-II)
- [Chapter III: Build the data]({{ site.url }}/data/anvio-structure/chapter-III)
- [Chapter IV: Analyze the data]({{ site.url }}/data/anvio-structure/chapter-IV) ← _you are here_
- [Chapter V: Reproduce every number]({{ site.url }}/data/anvio-structure/chapter-V)

## Analysis X: Setting the stage

From here on out, you won't find any more 'steps' to complete, but rather, 'analyses' you may _choose_ to complete. There is no longer an order to the workflow. Instead, run around and carry out any analyses you deem interesting. You have all the data prepared and ready to tackle any of the following analyses, in whichever order you like.

Throughout the lifespan of this paper, my strategy for carrying out analyses was to create a global R environment that would be loaded for any and all figure generation. FIXME

## Analysis X: Read recruitment summary (21 genomes)

In this analysis, we create Table S_RR, which provides summary-level recruitment information about each of the 21 SAR11 genomes that were used in the read recruitment experiment, including HIMB83.

As a reminder, we used Bowtie2 to recruit reads from each metagenome/metatranscriptome to each of the 21 SAR11 genomes in a competitive manner. The complete mapping information is stored in a series of {% include ARTIFACT name="bam-file" text="bam-files" %} present in `04_MAPPING/`, and the pertinent information from all samples and all genomes has been summarized into the {% include ARTIFACT name="profile-db" %} in `06_MERGED/`. To retrieve recruitment statistics for each genome, we can use the program {% include PROGRAM name="anvi-summarize" %}.

<div class="extra-info" style="{{ command_style  }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
anvi-summarize -p 06_MERGED/SAR11_clade/PROFILE.db \
               -c 03_CONTIGS/SAR11_clade-contigs.db \
               -C GENOMES \
               --init-gene-coverages \
               -o 07_SUMMARY_ALL
```
‣ **Time:** 18 min  
‣ **Storage:** 260 Mb  
</div> 

By default, {% include PROGRAM name="anvi-summarize" %} calculates coverage statistics averaged over {% include ARTIFACT name="bin" text="bins" %} (in our case each bin is a genome). But with the flag `--init-gene-coverages`, {% include PROGRAM name="anvi-summarize" %} takes the extra time to also report per-gene coverage statistics.

When finished, {% include PROGRAM name="anvi-summarize" %} produces a {% include ARTIFACT name="summary" %}, that offers some extensive reporting with tab-delimited files that you can open in Excel or Python/R. Here is the directory structure:

Quite simply, Table S_RR is a copy-paste job of a selection of these files. To create the Excel table, run the script `ZZ_SCRIPTS/gen_table_rr.py`.

<details markdown="1"><summary>Show/Hide Script</summary>
```python
#! /usr/bin/env python

import pandas as pd
from pathlib import Path

tables_dir = Path('WW_TABLES')
tables_dir.mkdir(exist_ok=True)

cov = pd.read_csv("07_SUMMARY_ALL/bins_across_samples/mean_coverage.txt", sep='\t')
q2q3 = pd.read_csv("07_SUMMARY_ALL/bins_across_samples/mean_coverage_Q2Q3.txt", sep='\t')
det = pd.read_csv("07_SUMMARY_ALL/bins_across_samples/detection.txt", sep='\t')
per = pd.read_csv("07_SUMMARY_ALL/bins_across_samples/bins_percent_recruitment.txt", sep='\t')
himb083_genes = pd.read_csv("07_SUMMARY_ALL/bin_by_bin/HIMB083/HIMB083-mean_coverage.txt", sep='\t')

with pd.ExcelWriter(tables_dir/'RR.xlsx') as writer:
    cov.to_excel(writer, sheet_name='Coverage')
    q2q3.to_excel(writer, sheet_name='Coverage Q2Q3 (inner quartiles)')
    det.to_excel(writer, sheet_name='Detection')
    per.to_excel(writer, sheet_name='% recruitment')
    himb083_genes.to_excel(writer, sheet_name='HIMB083 gene coverages')
```
</details> 

<div class="extra-info" style="{{ command_style  }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
python ZZ_SCRIPTS/gen_table_rr.py
```
‣ **Time:** Minimal  
‣ **Storage:** Minimal  
‣ **Memory:** Minimal  
</div> 

This creates the table from the paper and plops it into a directory `WW_TABLES`, which stores all tables from the paper.

## Analysis X: Comparing sequence similarity regimes

This analysis is a behind-the-scenes of the supplemental information entitled "_Regimes of sequence similarity probed by metagenomics, SAR11 cultured geomes, and protein families_", and provides explicit reproducibility steps to create Figure S_PS. Basically, we need to estimate the percent similarity from read recruitment results, from pangenomic comparisons, and from the Pfams that HIMB83 genes match to. Given the eclectic data sources, this is a rather lengthy process that I'll break up into 3 parts: (1) read recruitment, (2) pangenome, and (3) Pfam. Each of these steps creates a file `18_PERCENT_ID*.txt` that forms the data for each of the 3 histograms in Figure S_PS.

### (1) HIMB83 read recruitment

Here is the relevant description from the supplemental information.

<blockquote>
[Percent similarity (PS)] values for each gene were calculated by considering one metagenome at a time. In each metagenome, the reads that aligned to the gene were captured, trimmed (so there were no reads overhanging the gene), and compared to the aligned segment of HIMB83. The PS was calculated by comparing non-gap positions. This was then averaged to yield a PS value for each gene-metagenome pair. To define a single PS value for each gene, PS values were averaged across metagenomes.
<div class="blockquote-author">
  <b>Kiefl et al., November 2021 draft</b>
</div>
</blockquote>

Since anvi'o does not store per-read information, this data has to be grabbed from the {% include ARTIFACT name="bam-file" text="bam-files" %} themselves. I wrote a program `ZZ_SCRIPTS/analysis_gene_percent_id.py` to do this.

<details markdown="1"><summary>Show/Hide Script</summary>
```python
#! /usr/bin/env python

import numpy as np
import pandas as pd
import argparse
import anvio.dbops as dbops
import anvio.bamops as bamops

ap = argparse.ArgumentParser()

ap.add_argument('--contigs-db', help="Filepath to the contigs db", required=True)
ap.add_argument('--goi', help="Genes of interest (file of gene caller ids)", required=True)
ap.add_argument('--bams', help="File of bampaths. First column is sample name, second column is bam path. Header should be sample_id\tpath", required=True)
ap.add_argument('--output', help="Output text file", required=True)

args = ap.parse_args()

# ---------------------------------------------------

bams = pd.read_csv(args.bams, sep='\t')
bams = dict(zip(bams['sample_id'], bams['path']))

cdb = dbops.ContigsDatabase(args.contigs_db)
gene_calls = cdb.db.get_table_as_dataframe('genes_in_contigs').set_index('gene_callers_id')
gene_calls = gene_calls[gene_calls.index.isin([int(x.strip()) for x in open(args.goi).readlines()])]

d = {sample_id: [] for sample_id in bams.keys()}
i = 0
for sample_id, bam_path in bams.items():
    bam = bamops.BAMFileObject(bam_path)

    j = 0
    for gene_id, row in gene_calls.iterrows():
        print(f"Sample {sample_id} ({i}/{len(bams)}); Gene {gene_id} ({j}/{gene_calls.shape[0]})")

        percent_ids = []
        for read in bam.fetch_and_trim(row['contig'], row['start'], row['stop']):
            read.vectorize()
            percent_ids.append(read.get_percent_id())

        percent_ids = np.array(percent_ids)
        d[sample_id].append(np.mean(percent_ids))

        j += 1
    i += 1

    bam.close()

d['gene_callers_id'] = gene_calls.index.tolist()
pd.DataFrame(d).to_csv(args.output, sep='\t', index=False)
```
</details> 

Given a bunch of genes and a bunch of bam files, this script spits out the average percent identity of reads mapping to each gene in each sample. But before it can be ran, we need a list of {% include ARTIFACT name="bam-file" %} paths, which is created with another script, `ZZ_SCRIPTS/analysis_gen_mgx_bam_paths.py`.

<details markdown="1"><summary>Show/Hide Script</summary>
```python
#! /usr/bin/env python

import pandas as pd
from pathlib import Path

bam_dir = Path('04_MAPPING/SAR11_clade')
soi = [sample.strip() for sample in open('soi').readlines()]

bam_paths = dict(
    sample_id = [],
    path = [],
)

for bam_path in bam_dir.glob("*.bam"):
    sample_id = bam_path.stem
    if sample_id not in soi:
        continue
    bam_paths['sample_id'].append(sample_id)
    bam_paths['path'].append(str(bam_path))

pd.DataFrame(bam_paths).to_csv("mgx_bam_paths", sep='\t', index=False)
```
</details> 

There is nothing special about this script, it simply creates the following file, which points to each of the metagenomic {% include ARTIFACT name="bam-file" text="bam-files" %} corresponding to the samples of interest `soi` (generated from Step X).

|**sample_id**|**path**|
|:--|:--|
|IOS_64_05M|04_MAPPING/SAR11_clade/IOS_64_05M.bam|
|ION_39_25M|04_MAPPING/SAR11_clade/ION_39_25M.bam|
|IOS_64_65M|04_MAPPING/SAR11_clade/IOS_64_65M.bam|
|IOS_65_30M|04_MAPPING/SAR11_clade/IOS_65_30M.bam|
|IOS_65_05M|04_MAPPING/SAR11_clade/IOS_65_05M.bam|
|ION_36_17M|04_MAPPING/SAR11_clade/ION_36_17M.bam|
|ANW_142_05M|04_MAPPING/SAR11_clade/ANW_142_05M.bam|
|RED_32_80M|04_MAPPING/SAR11_clade/RED_32_80M.bam|
|ASW_78_05M|04_MAPPING/SAR11_clade/ASW_78_05M.bam|
|(...)|(...)|

Ok, first, generate the file `mgx_bam_paths`:

<div class="extra-info" style="{{ command_style  }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
python ZZ_SCRIPTS/analysis_gen_mgx_bam_paths.py
```
‣ **Time:** Minimal  
‣ **Storage:** Minimal  
</div> 

Then, run `ZZ_SCRIPTS/analysis_gene_percent_id.py`:

<div class="extra-info" style="{{ command_style  }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
python ZZ_SCRIPTS/analysis_gene_percent_id.py --bams mgx_bam_paths \
                                              --contigs-db 03_CONTIGS/SAR11_clade-contigs.db \
                                              --goi goi \
                                              --output 18_PERCENT_ID.txt
```
‣ **Time:** 150 min  
‣ **Storage:** Minimal  
</div> 

This creates the file `18_PERCENT_ID.txt` which is a table quantifying the average percent identity of reads for each gene in each sample.

### (2) SAR11 pangenome sequence similarity

Next up, the sequence similarity observed between SAR11 orthologs from the pangenome. Here is the relevant section in the supplemental:

<blockquote>
[G]ene clusters were calculated for HIMB83 and 20 additional SAR11 isolates using the anvi’o pangenomic workflow. An MSA was built from the sequences of each gene cluster using muscle, and then each non-HIMB83 sequence was compared to the HIMB83 sequence. The PS was determined by calculating the fraction of matches in non-gap positions. Each HIMB83 gene was attributed a single PS value by averaging PS values in each pairwise comparison, weighted by the number of non-gap positions in the pairwise alignment. Gene clusters containing multiple HIMB83 genes were ignored.
<div class="blockquote-author">
  <b>Kiefl et al., November 2021 draft</b>
</div>
</blockquote>

First, the IDs for all gene clusters that housed a HIMB83 gene of interest were collected and stored in a file called `gcoi` (gene clusters of interest), thanks to the work out `ZZ_SCRIPTS/get_HIMB83_gene_clusters.py`.

<details markdown="1"><summary>Show/Hide Script</summary>
```python
#! /usr/bin/env python

import pandas as pd

goi = [int(x.strip()) for x in open('goi').readlines()]
df = pd.read_csv("07_SUMMARY_PAN/SAR11_gene_clusters_summary.txt", sep='\t')

gc = df.loc[(df['genome_name'] == 'HIMB083') & (df['gene_callers_id'].isin(goi)), 'gene_cluster_id']
gc = gc.value_counts()[gc.value_counts() == 1].index.tolist()

with open('gcoi', 'w') as f:
    f.write("\n".join([str(x) for x in gc]))
```
</details> 

<div class="extra-info" style="{{ command_style  }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
python ZZ_SCRIPTS/get_HIMB83_gene_clusters.py
```
‣ **Time:** Minimal  
‣ **Storage:** Minimal  
</div> 

Then, alignments of the DNA sequences for each gene cluster are created with `ZZ_SCRIPTS/get_HIMB83_gene_cluster_alignments.sh`, which uses MUSCLE to perform the alignments:

<details markdown="1"><summary>Show/Hide Script</summary>
```bash
#! /usr/bin/env bash

rm -rf 18_HIMB83_CORE_COMPARED_TO_PANGENOME
mkdir -p 18_HIMB83_CORE_COMPARED_TO_PANGENOME

cat gcoi | while read gc; do
    anvi-get-sequences-for-gene-clusters -p 07_PANGENOME/PANGENOME/SAR11-PAN.db \
                                         -g 07_PANGENOME/SAR11-GENOMES.db \
                                         --report-DNA \
                                         -o 18_HIMB83_CORE_COMPARED_TO_PANGENOME/$gc.fa \
                                         --gene-cluster-id $gc \
                                         --min-num-genomes-gene-cluster-occurs 5
    muscle -in 18_HIMB83_CORE_COMPARED_TO_PANGENOME/$gc.fa -out 18_HIMB83_CORE_COMPARED_TO_PANGENOME/$gc.fa
done
```
</details> 

Go ahead and run this (it will take some time).

<div class="extra-info" style="{{ command_style  }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
bash ZZ_SCRIPTS/get_HIMB83_gene_cluster_alignments.sh
```
‣ **Time:** 100 min  
‣ **Storage:** 16 Mb  
</div> 

This will populate the directory `18_HIMB83_CORE_COMPARED_TO_PANGENOME/` with a bunch of gene cluster multiple sequence alignments (MSAs), which will be used to calculate percent similarity of HIMB83 genes to all of the other orthologs. The final script for pangenomic comparisons looks at each of the MSAs in `18_HIMB83_CORE_COMPARED_TO_PANGENOME/` and calculates the percent of matches in non-gap regions between the HIMB83 gene and the orthologs. To parse the MSAs, I made use of [ProDy](http://prody.csb.pitt.edu/) in the script `ZZ_SCRIPTS/analysis_get_percent_id_from_msa.py`.

<details markdown="1"><summary>Show/Hide Script</summary>
```python
#! /usr/bin/env python

import pandas as pd

from prody.sequence.msafile import parseMSA
from pathlib import Path

genome_name = 'HIMB083'
msa_paths = Path('18_HIMB83_CORE_COMPARED_TO_PANGENOME/').glob('*.fa')

dd = {
    'gene_callers_id': [],
    'percent_id': [],
}

for path in msa_paths:
    msa = parseMSA(str(path), format='FASTA')

    for i, label in enumerate(msa._labels):
        if genome_name in label:
            ref_id = i
    is_genome_of_interest = [genome_name in defline for defline in msa._labels]
    reference = msa._msa[is_genome_of_interest,:][0]
    ref_label = msa._labels[ref_id]
    ref_gene_id = int(ref_label.split('|')[-1].split(':')[-1])

    labels = [defline for defline in msa._labels if genome_name not in defline]
    sequences = msa._msa[[not x for x in is_genome_of_interest],:]
    sequences = dict(zip(labels, sequences))

    d = {
        'gene_callers_id': [],
        'match': [],
        'mismatch': [],
        'total': [],
    }

    total_nucleotides = 0
    for label, sequence in sequences.items():
        match, mismatch = 0, 0

        for ref, seq in zip(reference, sequence):
            ref, seq = ref.decode('utf-8'), seq.decode('utf-8')

            if ref == '-' or seq == '-':
                # only compare fully aligned nucleotides
                continue

            if ref == seq:
                match += 1
            else:
                mismatch += 1

        d['gene_callers_id'].append(ref_gene_id)
        d['match'].append(match)
        d['mismatch'].append(mismatch)
        d['total'].append(match + mismatch)

        total_nucleotides += match + mismatch

    df = pd.DataFrame(d)
    df['percent_id'] = df['match']/df['total']
    df['weight'] = df['total']/df['total'].sum()

    avg_percent_id = (df['percent_id']*df['weight']).sum() * 100

    dd['gene_callers_id'].append(ref_gene_id)
    dd['percent_id'].append(avg_percent_id)

df = pd.DataFrame(dd)
df.to_csv("18_PERCENT_ID_PANGENOME.txt", sep="\t", index=False)
```
</details> 

Give it a run when you're ready.

<div class="extra-info" style="{{ command_style  }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
python ZZ_SCRIPTS/analysis_get_percent_id_from_msa.py
```
‣ **Time:** Minimal  
‣ **Storage:** Minimal  
</div> 

Finally, we get the file `18_PERCENT_ID_PANGENOME.txt` which quantifies how similar each HIMB83 gene is to the correspondingorthologs in the SAR11 pangeome.

### (3) Protein family similarity

The third and last sequence comparison regime is between HIMB83 genes and genes of the Pfams they match to. Here is the relevant section in the supplemental:

<blockquote>
HIMB83 genes were matched to Pfam protein families via the anvi’o program `anvi-run-pfams`. Hits that passed the GA gathering threshold were retained, and the best hit (lowest e-value) for each HIMB83 gene was defined as the associated Pfam. For each HIMB83 gene, the associated Pfam seed sequence MSA was downloaded using the Python package prody and the HIMB83 protein sequence was added to the MSA using muscle. PS values were calculated from the MSAs in a manner identical to that outlined in (b). It is important to note that this comparison used protein sequences, whereas (a) and (b) both used nucleotide sequences.
<div class="blockquote-author">
  <b>Kiefl et al., November 2021 draft</b>
</div>
</blockquote>

Because Pfams span much larger evolutionary scales, and simply for convenience, this comparison was done with amino acid sequences rather than nucleotides. Thanks to the ProDy API, I was able to package this process into a single script, `ZZ_SCRIPTS/analysis_get_percent_id_from_pfam_msa.py`.

<details markdown="1"><summary>Show/Hide Script</summary>
```python
#! /usr/bin/env python

from pathlib import Path
from prody.database.pfam import fetchPfamMSA
from prody.sequence.msafile import parseMSA

import numpy as np
import pandas as pd
import argparse
import subprocess
import anvio.utils as utils
import anvio.dbops as dbops

ap = argparse.ArgumentParser()

ap.add_argument('--contigs-db', help="Filepath to the contigs db", required=True)
ap.add_argument('--genome-name', help="e.g. HIMB083", required=True)

args = ap.parse_args()

# --------------------------------------------------------------------------------------

cdb = dbops.ContigsDatabase(args.contigs_db)
functions = cdb.db.get_table_as_dataframe('gene_functions')
cdb.disconnect()

def get_best_pfam(gene_id):
    try:
        return functions.\
            query(f'gene_callers_id == {gene_id} & source == "Pfam"').\
            sort_values(by='e_value').\
            iloc[0, :]\
            ['accession']
    except IndexError:
        return None

# --------------------------------------------------------------------------------------

dd = {
    'gene_callers_id': [],
    'percent_id': [],
}

goi = [int(x.strip()) for x in open("goi").readlines()]
contigs_db = dbops.ContigsSuperclass(args)
for i, gene_id in enumerate(goi):
    print(f"Gene {gene_id} ({i}/{len(goi)})")
    pfam = get_best_pfam(gene_id)
    if not pfam:
        continue

    Path("18_HIMB83_CORE_COMPARED_TO_PFAM").mkdir(exist_ok=True)
    gene_fasta = Path("18_HIMB83_CORE_COMPARED_TO_PFAM") / f"{gene_id}.fa"
    pfam_msa_path = Path('18_HIMB83_CORE_COMPARED_TO_PFAM') / f"{pfam}_msa_for_{gene_id}.fa"
    default_dump_path = Path(f'{pfam}_seed.fasta')

    # Get the Pfam MSA
    try:
        pfam_versionless = ''.join(pfam.split('.')[:-1])
        url = f"https://pfam.xfam.org/family/{pfam_versionless}/alignment/seed/format?format=fasta&alnType=seed&order=t&case=u&gaps=dashes&download=1"
        utils.download_file(url=url, output_file_path=str(default_dump_path), check_certificate=False)
    except:
        print("FAILED")
        continue

    default_dump_path.replace(pfam_msa_path)

    # Get the gene AA sequence
    contigs_db.get_sequences_for_gene_callers_ids( gene_caller_ids_list=[gene_id], output_file_path=str(gene_fasta), report_aa_sequences=True)

    # massage gene into the MSA
    cmdline = f'muscle -profile -in1 {gene_fasta} -in2 {pfam_msa_path} -out {pfam_msa_path}'
    muscle_stdout = subprocess.call(cmdline, shell=True)

    # load the final MSA
    msa = parseMSA(str(pfam_msa_path), format='FASTA')

    for i, label in enumerate(msa._labels):
        if args.genome_name in label:
            ref_id = i
    is_genome_of_interest = [args.genome_name in defline for defline in msa._labels]
    reference = msa._msa[is_genome_of_interest,:][0]

    labels = [defline for defline in msa._labels if args.genome_name not in defline]
    sequences = msa._msa[[not x for x in is_genome_of_interest],:]
    sequences = dict(zip(labels, sequences))

    d = {
        'gene_callers_id': [],
        'match': [],
        'mismatch': [],
        'total': [],
    }

    total_nucleotides = 0
    for label, sequence in sequences.items():
        match, mismatch = 0, 0

        for ref, seq in zip(reference, sequence):
            ref, seq = ref.decode('utf-8'), seq.decode('utf-8')

            if ref == '-' or seq == '-':
                # only compare fully aligned nucleotides
                continue

            if ref == seq:
                match += 1
            else:
                mismatch += 1

        d['gene_callers_id'].append(gene_id)
        d['match'].append(match)
        d['mismatch'].append(mismatch)
        d['total'].append(match + mismatch)

        total_nucleotides += match + mismatch

    df = pd.DataFrame(d)
    df['percent_id'] = df['match']/df['total']
    df['weight'] = df['total']/df['total'].sum()

    avg_percent_id = (df['percent_id']*df['weight']).sum() * 100

    dd['gene_callers_id'].append(gene_id)
    dd['percent_id'].append(avg_percent_id)

df = pd.DataFrame(dd)
df.to_csv("18_PERCENT_ID_PFAM.txt", sep="\t", index=False)
```
</details> 

This script downloads the seed sequence MSA of hundreds of Pfams, so an internet connection and some patience is required. Run it like so:

<div class="extra-info" style="{{ command_style  }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
python ZZ_SCRIPTS/analysis_get_percent_id_from_pfam_msa.py --contigs-db CONTIGS.db --genome-name HIMB083
```
‣ **Time:** 30 min  
‣ **Storage:** 75 Mb  
‣ **Internet::** Yes  
</div> 

### Putting it all together

At this point, you should have the following 3 files:

```
18_PERCENT_ID.txt
18_PERCENT_ID_PANGENOME.txt
18_PERCENT_ID_PFAM.txt
```

Each file holds a distribution of percent similarity scores that can be visualized with the R script `ZZ_SCRIPTS/figure_s_ps.R`.


<details markdown="1"><summary>Show/Hide Script</summary>
```R
#! /usr/bin/env Rscript

source(file.path("utils.R"))

library(tidyverse)

args <- list()
args$reads <- "../18_PERCENT_ID.txt"
args$pfam <- "../18_PERCENT_ID_PFAM.txt"
args$pangenome <- "../18_PERCENT_ID_PANGENOME.txt"
args$output <- "../YY_PLOTS/FIG_S_PS"

dir.create(args$output, showWarnings=F, recursive=T)

# -----------------------------------------------------------------------------
# Reading in the data
# -----------------------------------------------------------------------------

reads <- read_tsv(args$reads) %>%
    pivot_longer(cols=-gene_callers_id) %>%
    group_by(gene_callers_id) %>%
    summarise(reads = mean(value))
pfam <- read_tsv(args$pfam) %>%
    rename(pfam = percent_id)
pan <- read_tsv(args$pangenome) %>%
    rename(pan = percent_id)
df_wide <- reads %>%
    left_join(pfam) %>%
    left_join(pan)
df <- df_wide %>%
    pivot_longer(cols=-gene_callers_id)

# -----------------------------------------------------------------------------
# Make the plot
# -----------------------------------------------------------------------------

g <- ggplot(data = df) +
    geom_histogram(aes(x=value, group=name, fill=name), bins=100, alpha=1.0) +
    scale_fill_manual(
        labels = c("Pangenome gANI", "Pfam", "Metagenomics rANI"),
        values=c("#AAAAAA", "#444444", "#B66B77")
    ) +
    labs(
        x = "Percent similarity (%)",
        y = "Number of genes",
        fill = "Comparison"
    ) +
    scale_y_continuous(limits = c(0,NA), expand = c(0.005, 0)) +
    theme_classic(base_size=12) +
    theme(
        text=element_text(size=12, family="Helvetica", face="bold"),
        legend.position = c(0.3, 0.7)
    )
w <- 5
display(g, file.path(args$output, "Figure_SPS.png"), h=w/2, w=w, as.png=TRUE)

options(pillar.sigfig = 5)
print(df %>% group_by(name) %>% summarise(mean=mean(value, na.rm=TRUE), median=median(value, na.rm=TRUE)))

```
</details> 

To run this script from the command line, you should do the following.

<div class="extra-info" style="{{ command_style  }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
cd ZZ_SCRIPTS
Rscript figure_s_ps.R
cd ..
```
‣ **Time:** Minimal  
‣ **Storage:** Minimal  

{:.notice}
It is bizarre, yet necessary, to sandwich the `Rscript` command between the two `cd` commands. This is ultimately because I am untalented at R.
</div> 

The output image is `YY_PLOTS/FIG_S_PS/Figure_SPS.png`.

[![s_ps]({{images}}/s_ps.png)]( {{images}}/s_ps.png){:.center-img .width-90}

## Analysis X: Codon usage between SAR11 genomes

<div class="extra-info" style="{{ command_style  }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
bash ZZ_SCRIPTS/codon_freqs_per_genome.sh
```
‣ **Time:** Minimal  
‣ **Storage:** Minimal  
</div> 

The above command generates `codon_freqs_per_genome.txt`. To generate Figure S_GNM_CDN_F, run

<div class="extra-info" style="{{ command_style  }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
cd ZZ_SCRIPTS/
Rscript figure_s_gnm_cdn_f.R
cd ..
```
‣ **Time:** Minimal  
‣ **Storage:** Minimal  
</div> 


## Analysis X: Distributions of environmental parameters

This is how I created Figure S_ENV.


<details markdown="1"><summary>Show/Hide Script</summary>
```R
#! /usr/bin/env Rscript

source("utils.R")

args <- list()
args$output <- "../YY_PLOTS/FIG_S_ENV"
args$meta <- "../TARA_metadata.txt"
args$soi <- "../soi"

library(tidyverse)

dir.create(args$output, showWarnings=F, recursive=T)

# -----------------------------------------------------------------------------
# Reading in the data
# -----------------------------------------------------------------------------

soi <- read_tsv(args$soi, col_names = F)

set.seed(4323)
get_pallette <- function(size) {
    colors <- list()
    for (i in 1:size) {
        colors[[i]] <- paste("#", paste(rep(sample(c("3", "4", "5", "6", "7", "8", "9", "A", "B"), 1), 6), collapse=""), sep="")
    }
    colors
}
cbPalette <- get_pallette(10)

df <- read_tsv(args$meta) %>%
    rename(sample_id=`Sample Id`) %>%
    pivot_longer(cols=c(
        `Depth (m)`,
        `Chlorophyll Sensor s`,
        `Temperature (deg C)`,
        `Salinity (PSU)`,
        `Oxygen (umol/kg)`,
        `Nitrates (umol/L)`,
        `PO4 (umol/L)`,
        `SI (umol/L)`
    ))

# -----------------------------------------------------------------------------
# Plot the thing
# -----------------------------------------------------------------------------

g <- ggplot(data = df, aes(x=value)) +
    geom_histogram(aes(fill=name)) +
    facet_wrap(~ name, ncol=3, scales="free", strip.position="bottom") +
    scale_fill_manual(values=cbPalette) +
    labs(
        x="",
        y="Number of Metagenomes"
    ) +
    theme_classic(base_size=12) +
    theme(
        legend.position = "none",
        text=element_text(size=12, family="Helvetica", face="bold"),
        strip.background = element_blank(),
        strip.placement="outside"
    ) +
    scale_y_continuous(limits = c(0,NA), expand = c(0.005, 0))
display(g, file.path(args$output, "meta.png"), width=6.5, height=5, as.png=TRUE)
```
</details> 

<div class="extra-info" style="{{ command_style  }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
cd ZZ_SCRIPTS/
Rscript figure_s_env.R
cd ..
```
‣ **Time:** Minimal  
‣ **Storage:** Minimal  
</div> 

The output image is `YY_PLOTS/FIG_S_ENV/meta.png`.

[![s_env]({{images}}/s_env.png)]( {{images}}/s_env.png){:.center-img .width-90}

## Analysis X: pN$^{(site)}$ variation across genes and samples

This is how I created Figure S_PN_HIST.

<details markdown="1"><summary>Show/Hide Script</summary>
```R
#! /usr/bin/env Rscript

request_scvs <- TRUE
request_regs <- FALSE
withRestarts(source("load_data.R"), terminate=function() message('load_data.R: data already loaded. Nice.'))

library(tidyverse)
library(latex2exp)

args <- list()
args$output <- "../YY_PLOTS/FIG_S_PN_HIST"

dir.create(args$output, showWarnings=F, recursive=T)

# -----------------------------------------------------------------------------
# Just do it
# -----------------------------------------------------------------------------

col1 <- '#888888'
col2 <- 'red'
coeff <- 12
g <- ggplot() +
    geom_histogram(
        data = scvs %>% group_by(sample_id) %>% summarize(x=sd(pN_popular_consensus, na.rm=T)) %>% mutate(group='group'),
        mapping = aes(x=log10(x), y=..density../coeff, fill=group),
        fill = col2,
        alpha = 0.5,
        bins = 100
    ) +
    geom_histogram(
        data = scvs %>% group_by(unique_pos_identifier) %>% summarize(x=sd(pN_popular_consensus)) %>% mutate(group='group'),
        mapping = aes(x=log10(x), y=..density.., fill=group),
        fill = col1,
        alpha = 0.5,
        bins = 100
    ) +
    theme_classic() +
    theme(
        text=element_text(size=10, family="Helvetica", face="bold"),
        axis.title.y = element_text(color = col1, size=13),
        axis.title.y.right = element_text(color = col2, size=13)
    ) +
    labs(
        y="Density",
        x=TeX("$log_{10}$(standard deviation)", bold=TRUE)
    ) +
    scale_y_continuous(
      name = "Density",
      sec.axis = sec_axis(~.*coeff, name="Density"),
      expand = c(0, 0)
    )
s <- 0.9
display(g, file.path(args$output, "fig.png"), width=s*3.5, height=s*2.4)
```
</details> 

<div class="extra-info" style="{{ command_style  }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
cd ZZ_SCRIPTS/
Rscript figure_s_pn_hist.R
cd ..
```
‣ **Time:** FIXME  
‣ **Storage:** Minimal  
</div> 

The output image is `YY_PLOTS/FIG_S_PN_HIST/fig.png`.

[![pn_hist]({{images}}/pn_hist.png)]( {{images}}/pn_hist.png){:.center-img .width-90}

## Analysis X: Comparing AlphaFold to MODELLER

I started working on this study a few years before the exceedingly recent revolution in structure prediction that has been seeded by AlphaFold and its monumental success seen during the [CASP14](https://predictioncenter.org/casp14/) structure prediction competition. Back then, I developed a program {% include PROGRAM name="anvi-gen-structure-database" %} that predicted protein structures using template-based homology modeling with [MODELLER](https://salilab.org/modeller/).

Template-based homology modeling relies on the existence of a template structure that shares homology with your gene of interest. I wrote a whole thing about it ([click here](https://merenlab.org/2018/09/04/getting-started-with-anvio-structure/)) that you should read if you're interested. In that post I explain everything there is to know about {% include PROGRAM name="anvi-gen-structure-database" %}, how it uses MODELLER to predict protein structures, and how you can interactively explore metagenomic sequence variants in the context of predicted protein structures and ligand-binding sites (we'll briefly get into this specific topic later). But **the gist of template-based homology modelling is that it only works if your protein of interest has a homologous protein with an experimentally-resolved structure**. The higher the percent similarity between your protein and the template protein, the more accurate the structure prediction will be.

And up until October 2021, all of the analyses in this paper used structures derived from template-based homology modelling. Then DeepMind released the source code for [AlphaFold](https://www.nature.com/articles/s41586-021-03819-2), which does _not_ require templates with solved structures, and outperformed all other methods that came before it. It was then that I realized it would be worthwhile to make the switch, and so I did.

With that in mind, the usage of MODELLER in this study is somewhat historical. Before AlphaFold, I actually wrote an entire pipeline of software encapsulated in the program {% include PROGRAM name="anvi-gen-structure-database" %} that predicts protein structures _en masse_, where the user simply provides a {% include ARTIFACT name="contigs-db" %} of interest. Since we have two methods to calculate protein structures, we thought it would be insightful to compare them, which led to the supplemental info entitled, "_Comparing structure predictions between AlphaFold and MODELLER_".

In this document I'll detail (a) how to predict MODELLER structures using {% include PROGRAM name="anvi-gen-structure-database" %} and (b) how I compared MODELLER structures to AlphaFold structures.

### Calculating MODELLER structures

Calculating structures with MODELLER using {% include PROGRAM name="anvi-gen-structure-database" %} is exceedingly easy and in a rather extensive [blog post](https://merenlab.org/2018/09/04/getting-started-with-anvio-structure/), I go into all of the nitty gritty. The net result is that generating structures for genes within a {% include ARTIFACT name="contigs-db" %} has never been easier. It boils down to just one command, which you should feel free to run with as many threads as you can afford to shell out.

<div class="extra-info" style="{{ command_style  }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
anvi-gen-structure-database -c CONTIGS.db -o 09_STRUCTURE_MOD.db -T <NUM_THREADS> --genes-of-interest goi
```
‣ **Time:** ~(18/`<NUM_THREADS>`) hours
‣ **Storage:** 160 Mb  
‣ **Internet:** Maybe  

{:.notice}
Too lazy? Just download the resulting {% include ARTIFACT name="structure-db" %}: `wget FIXME` FIXME this will not work unless they also use the wget CONTIGS.db due to hash mismatch

{:.notice}
If you don't have internet, you'll need an offline database built ahead of time that you can generate with {% include PROGRAM name="anvi-setup-pdb-database" %} (which will itself require internet).

</div>

Running the above command will produce a {% include ARTIFACT name="structure-db" %} of MODELLER-predicted structures.

### Trustworthy MODELLER structures

In Step X, we defined 'trustworthy' AlphaFold structures as those that maintain an average pLDDT score of >80, and we stored the corresponding gene IDs in the file `12_GENES_WITH_GOOD_STRUCTURES`.

Similarly, we take certain measures to try and define 'trustworthy' when using MODELLER structures that are outlined in the Methods section of the paper:

<blockquote>
We discarded any proteins if the best template had a percent similarity of <30%. Unlike more sophisticated homology approaches that make use of multi-domain templates (CITE RaptorX), we used single-domain templates which are convenient and are accurate up to several angstroms, yet can lead to physically inaccurate models when the templates’ domains match to some, but not all of the sequences’ domains. To avoid this, we discarded any templates if the alignment coverage of the protein sequence to the chosen template was <80%. Applying these filters resulted in 408 structures from the 1a.3.V core, which was further refined by requiring that the root mean squared distance (RMSD)  between the predicted structure and the most similar template did not exceed 7.5 angstroms, and that the GA341 model score exceeded 0.95. After applying these constraints, we were left with 346 structures in the 1a.3.V that we assumed to be ‘trustworthy’ structures as predicted by MODELLER
<div class="blockquote-author">
  <b>Kiefl et al., February 2022 draft</b>
</div>
</blockquote>

In the above command, the 30% template homology and 80% alignment coverage have already been applied (though you can change these default cutoffs with the flags `--percent-cutoff` and `--alignment-fraction-cutoff`, respectively), meaning all structures in `09_STRUCTURE_MOD.db` satisfy these constraints. To apply the RMSD and GA341 filters, I wrote the script `ZZ_SCRIPTS/gen_genes_with_good_structures_modeller.py`, which calculates the RMSD between each structure and its 'best' template using ProDy, where best template is defined by the template with the highest alignment coverage $\times$ percent similarity score. It also queries the GA341 score produced by MODELLER, and filters any models with scores <0.95.

<details markdown="1"><summary>Show/Hide gen_genes_with_good_structures_modeller.py</summary>
```python
#! /usr/bin/env python

import argparse
import anvio.structureops as sops
import anvio.filesnpaths as filesnpaths
import anvio.utils as utils

from prody.proteins.pdbfile import parsePDB
from prody.proteins.compare import matchChains
from prody.measure.transform import calcRMSD
from prody.measure.transform import calcTransformation

ap = argparse.ArgumentParser()
ap.add_argument('-s', '--structure')
args = ap.parse_args()

sdb = sops.StructureDatabase(args.structure)
templates = sdb.db.get_table_as_dataframe('templates')
models = sdb.db.get_table_as_dataframe('models')
goi = set([int(x.strip()) for x in open('goi').readlines()])
templates = templates[templates['corresponding_gene_call'].isin(goi)]

# downsize templates to only include the best template, where best template is defined as
# the template with the highest (alignment coverage * percent similarity = proper percent similarity)
templates = templates.groupby('corresponding_gene_call').apply(lambda df: df[df['proper_percent_similarity'] == df['proper_percent_similarity'].max()])
templates = templates.drop_duplicates(subset=['corresponding_gene_call', 'proper_percent_similarity'])

def get_RMSD(template_path, gene_path):
    x = parsePDB(template_path)
    y = parsePDB(gene_path)

    try:
        matches = matchChains(x, y, seqid=30, overlap=30)
    except:
        matches = None

    if matches is None:
        return matches

    calcTransformation(matches[0][0], matches[0][1]).apply(matches[0][0])
    return calcRMSD(matches[0][0], matches[0][1])

for i, row in templates.iterrows():
    template = row['pdb_id'] + row['chain_id']
    gene_id = row['corresponding_gene_call']

    # export structures
    try:
         template_path = utils.download_protein_structure(protein_code=row['pdb_id'], chain=row['chain_id'], output_path=filesnpaths.get_temp_file_path(), raise_if_fail=True)
    except:
        continue
    gene_path = sdb.export_pdb_content(gene_id, filesnpaths.get_temp_file_path(), ok_if_exists=True)

    templates.loc[i, 'RMSD'] = get_RMSD(template_path, gene_path)

models = models[models['picked_as_best'] == 1]
templates.reset_index(drop=True, inplace=True)
df = templates.merge(models, on='corresponding_gene_call', how='left')
df = df[(df.RMSD <= 7.5) & (df.GA341_score >= 0.95)]

with open('12_GENES_WITH_GOOD_STRUCTURES_MODELLER', 'w') as f:
    for gene_id in df['corresponding_gene_call']:
        f.write(f"{gene_id}\n")
```
</details> 

This should take about 5 minutes to complete, and outputs the file `12_GENES_WITH_GOOD_STRUCTURES_MODELLER`.

<div class="extra-info" style="{{ command_style  }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
python ZZ_SCRIPTS/gen_genes_with_good_structures_modeller.py -s 09_STRUCTURE_MOD.db
```
‣ **Time:** 5 min  
‣ **Storage:** Minimal  
‣ **Internet:** Yes  
</div> 

### Method comparison

You should now have these 2 files

```
12_GENES_WITH_GOOD_STRUCTURES
12_GENES_WITH_GOOD_STRUCTURES_MODELLER
```

which indicate which structure predictions are considered trustworthy by each method. To compare the two methods, I created a summary of alignment and structural metrics. Here is the list of metrics considered.

**Alignment metrics**. These are calculated by first aligning the structures determined from each method. Obviously, these metrics are only suitable in the subset of protein sequences where a trustworthy structure was determined via both methods.

1. RMSD - The [root-mean-square deviation](https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions) of alpha carbon (backbone) atoms. The units are angstroms.
2. TM score - The [template modeling score](https://en.wikipedia.org/wiki/Template_modeling_score). This is a popular metric designed to outperform other metrics of global similarity such as RMSD.
2. Contact map MAE - The [mean-absolute-error](https://en.wikipedia.org/wiki/Mean_absolute_error) between residue center-of-mass contact maps. This isn't a widely _used_ metric, just something I made up to de-emphasize RMSD losses introduced by loose linkers between domains. The units are angstroms.

**Structure metrics**. These are calculated directly from the protein structures and are calculated for each structure.

1. mean RSA - The mean relative solvent accessibility.
2. $R_g$ - The [radius of gyration](https://en.wikipedia.org/wiki/Radius_of_gyration). The units are angstroms.
3. End-to-end distance - This is the Euclidean distance from the N-terminus to the C-terminus. The units are angstroms.
3. Fraction of alpha helices - This is the fraction of a gene that corresponded to alpha helices. According to the 8-state secondary structure proposed by DSSP, these correspond to classes G, H, and I.
3. Fraction of beta strands - This is the fraction of a gene that corresponded to beta strands. According to the 8-state secondary structure proposed by DSSP, these correspond to classes E and B.
3. Fraction of loops/unstructured - This is the fraction of a gene that corresponded to loops or otherwise unclassified residues. According to the 8-state secondary structure proposed by DSSP, these correspond to classes S, T, and C.

**AlphaFold metrics**. These are specific to AlphaFold.

1. mean pLDDT - The mean [pLDDT](https://alphafold.ebi.ac.uk/faq)

**MODELLER metrics**. These are specific to MODELLER.

1. template similarity - The mean percent similarity of all templates used for the structure.

For each protein structure, I calculated all of these metrics with the script `ZZ_SCRIPTS/comp_struct_preds.py`, which takes about 15 minutes to run.

<details markdown="1"><summary>Show/Hide Script</summary>
```python
#! /usr/bin/env python

import anvio.utils as utils
import anvio.terminal as terminal
import anvio.structureops as sops
import anvio.filesnpaths as filesnpaths

from prody import confProDy
from prody.measure.measure import calcGyradius, calcDistance
from prody.proteins.pdbfile import parsePDB
from prody.proteins.compare import matchChains
from prody.measure.transform import calcRMSD, calcTransformation

import numpy as np
import pandas as pd
import tmscoring
import matplotlib.pyplot as plt

from pathlib import Path

# -------------------------------------------------------
# Load up the structure databases and auxiliary info
# -------------------------------------------------------

confProDy(verbosity='none')
progress = terminal.Progress()
run = terminal.Run()

db_AF = sops.StructureDatabase('09_STRUCTURE.db')
db_MOD = sops.StructureDatabase('09_STRUCTURE_MOD.db')

plddt = pd.read_csv('09_STRUCTURES_AF/pLDDT_gene.txt', sep='\t').set_index('gene_callers_id')
plddt_residue = pd.read_csv('09_STRUCTURES_AF/pLDDT_residue.txt', sep='\t')

# -------------------------------------------------------
# Get intersection considered trustworthy for both methods
# -------------------------------------------------------

good_AF = set([int(x.strip()) for x in open('12_GENES_WITH_GOOD_STRUCTURES').readlines()])
good_MOD = set([int(x.strip()) for x in open('12_GENES_WITH_GOOD_STRUCTURES_MODELLER').readlines()])
good_intersect = good_AF.intersection(good_MOD)

run.info('Number of trustworthy AlphaFold predictions', len(good_AF))
run.info('Number of trustworthy MODELLER predictions', len(good_MOD))
run.info('Number of shared trustworthy predictions', len(good_intersect))

# -------------------------------------------------------
# Define metrics for assessing ALIGNMENTS
# -------------------------------------------------------

def get_RMSD(path1, path2):
    struct1 = parsePDB(path1)
    struct2 = parsePDB(path2)
    matches = matchChains(struct1, struct2, seqid=30, overlap=30)
    calcTransformation(matches[0][0], matches[0][1]).apply(matches[0][0])
    return calcRMSD(matches[0][0], matches[0][1])

def get_TM_score(path1, path2):
    """Get TM score using tmscoring.

    For whatever reason, there is a bit of an issue with tmscoring. In a rare number of cases, the
    TM score returned depends upon the order of filepaths passed to tmscoring.TMscoring, and by
    trial and error, I have found that taking the max of these values consistently yields the same
    results as the original TM web server https://zhanggroup.org/TM-score/. As such, I return the
    maximum of both values here. See issue: https://github.com/Dapid/tmscoring/issues/6
    """
    aln = tmscoring.TMscoring(path1, path2)
    aln.optimise()
    TM_score = aln.tmscore(**aln.get_current_values())

    aln = tmscoring.TMscoring(path2, path1)
    aln.optimise()
    TM_score_2 = aln.tmscore(**aln.get_current_values())

    return max(TM_score, TM_score_2)

def get_contact_map_mean_absolute_error(gene_id):
    structure_AF = db_AF.get_structure(gene_id)
    contact_map_AF = get_contact_map(structure_AF)
    structure_MOD = db_MOD.get_structure(gene_id)
    contact_map_MOD = get_contact_map(structure_MOD)
    MAE = np.abs(contact_map_MOD - contact_map_AF).sum() / np.product(contact_map_MOD.shape)
    return MAE

def get_contact_map(structure):
    d = {}
    protein_length = len(structure.structure)
    contact_map = np.zeros((protein_length, protein_length))
    for i, residue1 in enumerate(structure.structure):
        if i not in d:
            d[i] = structure.get_residue_center_of_mass(residue1)
        COM1 = d[i]
        for j, residue2 in enumerate(structure.structure):
            if i > j:
                contact_map[i, j] = contact_map[j, i]
            else:
                if j not in d:
                    d[j] = structure.get_residue_center_of_mass(residue2)
                COM2 = d[j]
                contact_map[i, j] = np.sqrt(np.sum((COM1 - COM2)**2))
    return contact_map

# -------------------------------------------------------
# Define metrics for assessing STRUCTURES
# -------------------------------------------------------

def get_Rg(path):
    """Get the radius of gyration"""
    struct = parsePDB(path).select('not hydrogen')
    return calcGyradius(struct)

def get_end_to_end(path):
    """Get the radius of gyration"""
    struct = parsePDB(path).select('not hydrogen')
    nter, cter = struct[0], struct[-1]
    return calcDistance(nter, cter)

def get_mean_RSA(gene_id, db):
    return db.get_residue_info_for_gene(gene_id)['rel_solvent_acc'].mean()

def get_frac_sec_struct(gene_id, db):
    res_info = db.get_residue_info_for_gene(gene_id)
    alpha = res_info['sec_struct'].isin(['G','H','I']).sum()/len(res_info)
    beta = res_info['sec_struct'].isin(['E','B']).sum()/len(res_info)
    loop = res_info['sec_struct'].isin(['S','T','C']).sum()/len(res_info)
    return alpha, beta, loop

def get_length(gene_id, db1, db2):
    try:
        return len(db1.get_structure(gene_id).get_sequence())
    except:
        return len(db2.get_structure(gene_id).get_sequence())

# -------------------------------------------------------
# Run the metrics for each gene
# -------------------------------------------------------

d = dict(
    gene_callers_id = [],
    length = [],
    has_AF = [],
    has_MOD = [],
    RMSD = [],
    TM_score = [],
    contact_MAE = [],
    mean_RSA_AF = [],
    mean_RSA_MOD = [],
    gyradius_AF = [],
    gyradius_MOD = [],
    endtoend_AF = [],
    endtoend_MOD = [],
    frac_alpha_AF = [],
    frac_alpha_MOD = [],
    frac_beta_AF = [],
    frac_beta_MOD = [],
    frac_loop_AF = [],
    frac_loop_MOD = [],
    template_similarity = [],
    mean_pLDDT = [],
    mean_trunc_pLDDT = [],
)

progress.new('Calculating metrics', progress_total_items=len(good_AF.union(good_MOD)))

for gene_id in good_AF.union(good_MOD):
    has_AF = True if gene_id in good_AF else False
    has_MOD = True if gene_id in good_MOD else False

    path_AF = db_AF.export_pdb_content(gene_id, filesnpaths.get_temp_file_path()) if has_AF else np.nan
    # Calculate structure metrics
    mean_RSA_AF = get_mean_RSA(gene_id, db_AF) if has_AF else np.nan
    gyradius_AF = get_Rg(path_AF) if has_AF else np.nan
    endtoend_AF = get_end_to_end(path_AF) if has_AF else np.nan
    frac_alpha_AF, frac_beta_AF, frac_loop_AF = get_frac_sec_struct(gene_id, db_AF) if has_AF else (np.nan, np.nan, np.nan)
    # AF-specific metrics
    mean_pLDDT = plddt.loc[gene_id, 'plddt'] if has_AF else np.nan
    mean_trunc_pLDDT = plddt_residue.loc[plddt_residue['gene_callers_id']==gene_id, 'plddt'][10:-10].mean() if has_AF else np.nan

    path_MOD = db_MOD.export_pdb_content(gene_id, filesnpaths.get_temp_file_path()) if has_MOD else np.nan
    # Calculate structure metrics
    mean_RSA_MOD = get_mean_RSA(gene_id, db_MOD) if has_MOD else np.nan
    gyradius_MOD = get_Rg(path_MOD) if has_MOD else np.nan
    endtoend_MOD = get_end_to_end(path_MOD) if has_MOD else np.nan
    frac_alpha_MOD, frac_beta_MOD, frac_loop_MOD = get_frac_sec_struct(gene_id, db_MOD) if has_MOD else (np.nan, np.nan, np.nan)
    # Calculate MOD-specific metrics
    template_similarity = db_MOD.get_template_info_for_gene(gene_id)['percent_similarity'].mean() if has_MOD else np.nan

    # Calculate alignment metrics
    RMSD = get_RMSD(path_AF, path_MOD) if (has_MOD and has_AF) else np.nan
    TM_score = get_TM_score(path_AF, path_MOD) if (has_MOD and has_AF) else np.nan
    MAE = get_contact_map_mean_absolute_error(gene_id) if (has_MOD and has_AF) else np.nan

    progress.update(f"Gene ID {gene_id}, RMSD {RMSD:.1f}, contact MAE {MAE:.1f}, TM score {TM_score:.2f}")
    progress.increment()

    d['gene_callers_id'].append(gene_id)
    d['length'].append(get_length(gene_id, db_AF, db_MOD))
    d['has_AF'].append(has_AF)
    d['has_MOD'].append(has_MOD)
    d['RMSD'].append(RMSD)
    d['TM_score'].append(TM_score)
    d['contact_MAE'].append(MAE)
    d['mean_RSA_AF'].append(mean_RSA_AF)
    d['mean_RSA_MOD'].append(mean_RSA_MOD)
    d['gyradius_AF'].append(gyradius_AF)
    d['gyradius_MOD'].append(gyradius_MOD)
    d['endtoend_AF'].append(gyradius_AF)
    d['endtoend_MOD'].append(gyradius_MOD)
    d['frac_alpha_AF'].append(frac_alpha_AF)
    d['frac_alpha_MOD'].append(frac_alpha_MOD)
    d['frac_beta_AF'].append(frac_beta_AF)
    d['frac_beta_MOD'].append(frac_beta_MOD)
    d['frac_loop_AF'].append(frac_loop_AF)
    d['frac_loop_MOD'].append(frac_loop_MOD)
    d['template_similarity'].append(template_similarity)
    d['mean_pLDDT'].append(mean_pLDDT)
    d['mean_trunc_pLDDT'].append(mean_trunc_pLDDT)

progress.end()

df = pd.DataFrame(d)
df.to_csv('09_STRUCTURE_comparison.txt', sep='\t', index=False)

tables_dir = Path('WW_TABLES')
tables_dir.mkdir(exist_ok=True)
with pd.ExcelWriter(tables_dir/'STRUCT_COMP.xlsx') as writer:
    df.to_excel(writer, sheet_name='Structure comparison')

```
</details> 

<div class="extra-info" style="{{ command_style  }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
python ZZ_SCRIPTS/comp_struct_preds.py
```
‣ **Time:** 15 min  
‣ **Storage:** Minimal  
</div>

This creates the files `09_STRUCTURE_comparison.txt` and `WW_TABLES/STRUCT_COMP.xlsx`, which are otherwise known as Table S_FIXME. Using this table, I created Figure S_FIXME with the script `ZZ_SCRIPTS/figure_s_comp.R`.

<details markdown="1"><summary>Show/Hide Script</summary>
```R
#! /usr/bin/env Rscript

source(file.path("utils.R"))

library(tidyverse)
library(latex2exp)
library(cowplot)

args <- list()
args$input <- "../09_STRUCTURE_comparison.txt"
args$output <- "../YY_PLOTS/FIG_S_COMP"

dir.create(args$output, showWarnings=F, recursive=T)

df <- read_tsv(args$input) %>% mutate(frac_ss_MOD = frac_alpha_MOD + frac_beta_MOD, frac_ss_AF = frac_alpha_AF + frac_beta_AF)

size <- 1.1
#col1 <- '#D7C49E'
#col2 <- '#343148'
col1 <- '#00203F'
col2 <- '#70ba98'
col3 <- '#FC766A'
col4 <- '#5B84B1'

# -----------------------------------------------------------
# Measures of alignment (row 1)
# -----------------------------------------------------------

g1 <- ggplot(df %>% filter(!is.na(TM_score))) +
    geom_freqpoly(aes(TM_score), color=col1, size=size) +
    my_theme(8) +
    scale_y_continuous(expand=c(0,0)) +
    labs(x = 'TM score')

g2 <- ggplot(df %>% filter(!is.na(RMSD))) +
    geom_freqpoly(aes(RMSD), color=col1, size=size) +
    my_theme(8) +
    scale_y_continuous(expand=c(0,0)) +
    labs(x='RMSD (Å)')

# -----------------------------------------------------------
# Structure metric comparisons within intersection (row 2)
# -----------------------------------------------------------

plot_data <- df %>% filter(has_MOD)

g3 <- ggplot(plot_data) +
    geom_freqpoly(aes(x=frac_alpha_MOD+frac_beta_MOD), col=col1, size=size) +
    geom_freqpoly(aes(x=frac_alpha_AF+frac_beta_AF), col=col2, size=size) +
    my_theme(8) +
    scale_y_continuous(expand=c(0,0)) +
    labs(x="secondary structure fraction (SSF)")

g4 <- ggplot(data=df %>% filter(has_MOD, !is.na(TM_score))) +
    geom_boxplot(
        mapping=aes(x=TM_score<0.8, y=frac_ss_AF/frac_ss_MOD), fill=col1, alpha=0.3, size=0.7*size, outlier.size=0.8) +
    scale_x_discrete(labels=c("TM score > 0.8","TM score < 0.8")) +
    my_theme(8) +
    scale_y_continuous(expand=c(0,0)) +
    labs(x='', y=TeX("$SSF^{AlphaFold} / SSF^{MODELLER}$", bold=T))

# -----------------------------------------------------------
# Comparisons between AlphaFold cohorts (row 3)
# -----------------------------------------------------------

g5 <- ggplot() +
    geom_freqpoly(data=df %>% filter(has_MOD), mapping=aes(x=mean_pLDDT, y=..density..), col=col3, size=size) +
    geom_freqpoly(data=df %>% filter(!has_MOD), mapping=aes(x=mean_pLDDT, y=..density..), col=col4, size=size) +
    my_theme(8) +
    scale_y_continuous(expand=c(0,0)) +
    labs(x="Mean pLDDT")

g6 <- ggplot() +
    geom_freqpoly(data=df %>% filter(has_MOD), mapping=aes(x=length, y=..density..), col=col3, size=size) +
    geom_freqpoly(data=df %>% filter(!has_MOD), mapping=aes(x=length, y=..density..), col=col4, size=size) +
    my_theme(8) +
    scale_y_continuous(expand=c(0,0)) +
    labs(x="Protein sequence length")

# -----------------------------------------------------------
# Squish them together
# -----------------------------------------------------------

g <- plot_grid(g1, g2, g3, g4, g5, g6, nrow=3)
display(g, file.path(args$output, "fig.png"), width=6, height=5, as.png=T)

```
</details> 

<div class="extra-info" style="{{ command_style  }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
cd ZZ_SCRIPTS/
Rscript figure_s_comp.R
cd -
```
‣ **Time:** Minimal  
‣ **Storage:** Minimal  
</div> 

Running this creates Figure S_COMP under the filename `YY_PLOTS/FIG_S_COMP/fig.png`:

[![s_comp]({{images}}/s_comp.png)]( {{images}}/s_comp.png){:.center-img .width-70}

## Analysis X: pN- and pS-weighted RSA and DTL null distributions

In Figure 1, we calculated how per-site pN and pS distribute with respect to the variables RSA and DTL. That data looks like this:

[![fig1ab]({{images}}/fig1ab.png)]({{images}}/fig1ab.png){:.center-img .width-70}

For the purposes of explanation, let's just focus on RSA (left) and pN (red) for the time being. This is a weighted distribution calculated by weighting each residue's RSA by the pN observed at each sample. To gauge how this distribution should look if pN was completely uncorrelated with RSA, we shuffled the data, so that each residue (in each sample) was given the RSA of a randomly chosen residue (in a randomly chosen sample). To avoid biases introduced from a single shuffle, we calculated 10 null distributions and the one displayed in Figure 1 is the average of all 10.

The null distribution for pS-weighted RSA is created just the same way... Yet you should be asking yourself, why is there only one null distribution displayed in Figure 1a? The reason is purely aesthetic. As it turns out, the null distributions of pS-weighted and pN-weighted RSA are nearly identical, and so increase visual clarity I only displayed one, as is mentioned in the figure caption:

<blockquote>
Since the null distribution for pS(site) so closely resembles the null distribution for pN(site), it has been excluded for visual clarity, but can be seen in Figure S_SHUFF_COMP
<div class="blockquote-author">
  <b>Kiefl et al., February 2022 draft</b>
</div>
</blockquote>

In Figure S_SHUFF_COMP I explicitly compare these null distributions to show there's no sleight of hand. Here is the script that creates Figure S_SHUFF_COMP.

<details markdown="1"><summary>Show/Hide Script</summary>
```R
#! /usr/bin/env Rscript

request_scvs <- TRUE
request_regs <- FALSE
withRestarts(source("load_data.R"), terminate=function() message('load_data.R: data already loaded. Nice.'))

library(optparse)
library(tidyverse)
library(cowplot)

args <- list()
args$output <- "../YY_PLOTS/FIG_S_SHUFF_COMP"

# Create directory
dir.create(args$output, showWarnings=F, recursive=T)

# -----------------------------------------------------------------------------
# make the comparison function
# -----------------------------------------------------------------------------

get_plot <- function(variable, xlim=NA, ylim=NA,
                       sqrty=FALSE, replicates=10, shuffle=T) {

    g <- ggplot()
    count_method <- "density"
    bin_count <- 50

    type = 'pN_popular_consensus'
    plot_data <- scvs %>%
        select((!!sym(variable)), (!!sym(type))) %>%
        filter(!is.na((!!sym(type))), !is.na((!!sym(variable))))
    type_frac <- c()
    type_var <- c()
    for (i in 1:replicates) {
        shuffle_plot_data <- plot_data %>%
            sample_n(size=plot_data %>% dim() %>% .[[1]], replace=F)
        type_frac <- append(type_frac, shuffle_plot_data[,type] %>% .[[1]])
        type_var <- append(type_var, plot_data[,variable] %>% .[[1]])
        plot_data$shuffled <- shuffle_plot_data %>% .[[2]]
    }
    mean_shuffle <- data.frame(type_frac, type_var)

    g <- g +
        stat_bin(
            data = mean_shuffle,
            mapping=aes_string(x="type_var", weight="type_frac", y=paste("..", count_method, "..", sep="")),
            geom="step",
            center=0,
            size=0.75,
            color=ns_col,
            bins=bin_count,
            alpha=0.7
        )

    type = 'pS_popular_consensus'
    plot_data <- scvs %>%
        select((!!sym(variable)), (!!sym(type))) %>%
        filter(!is.na((!!sym(type))), !is.na((!!sym(variable))))
    type_frac <- c()
    type_var <- c()
    for (i in 1:replicates) {
        shuffle_plot_data <- plot_data %>%
            sample_n(size=plot_data %>% dim() %>% .[[1]], replace=F)
        type_frac <- append(type_frac, shuffle_plot_data[,type] %>% .[[1]])
        type_var <- append(type_var, plot_data[,variable] %>% .[[1]])
        plot_data$shuffled <- shuffle_plot_data %>% .[[2]]
    }
    mean_shuffle <- data.frame(type_frac, type_var)

    g <- g +
        stat_bin(
            data = mean_shuffle,
            mapping=aes_string(x="type_var", weight="type_frac", y=paste("..", count_method, "..", sep="")),
            geom="step",
            origin=0,
            size=0.75,
            color=s_col,
            bins=bin_count,
            alpha=0.7
        )

    g <- g +
        labs(y=count_method, x=variable) +
        theme_classic() +
        theme(
            text=element_text(size=11, family="Helvetica", face="bold")
        ) +
        scale_x_continuous(limits = xlim, expand = c(0, 0))
    if (sqrty) {
        g <- g + scale_y_continuous(limits = c(NA, ylim), trans='sqrt', expand = c(0.005, 0))
    } else {
        g <- g + scale_y_continuous(limits = c(NA, ylim), expand = c(0.005, 0))
    }
    g
}

N <- 10
DTL <- get_plot(variable='ANY_dist', xlim=c(0,40), sqrty=FALSE, replicates=N) +
    labs(y="Density", x="DTL (Å)")
RSA <- get_plot(variable='rel_solvent_acc', xlim=c(0,1), sqrty=TRUE, replicates=N) +
    labs(y="Density", x="RSA")

# -----------------------------------------------------------------------------
# Put it all together
# -----------------------------------------------------------------------------

display(
    plot_grid(RSA, DTL, ncol=2, align='v'),
    output = file.path(args$output, 'fig.png'),
    w = 7.5,
    h = 2.5
)
```
</details> 

Running the following creates Figure S_SHUFF_COMP under the filename `YY_PLOTS/FIG_S_SHUFF_COMP/fig.png`:


<div class="extra-info" style="{{ command_style  }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
cd ZZ_SCRIPTS/
Rscript figure_s_shuff_comp.R
cd -
```
‣ **Time:** ~1 hour  
‣ **Storage:** Minimal  
</div> 

[![shuff_comp]({{images}}/shuff_comp.png)]({{images}}/shuff_comp.png){:.center-img .width-70}

## Analysis X: An alternative 1D definition for DTL

Besides our Euclidean distance definition of DTL, we also considered a much more primitive distance metric, which was defined not in 3D space but by the distance in sequence. For example, if a gene had only one ligand-binding residue, which occurred at the fifth residue, then the 25th residue would have a DTL of 20. In Figure S_1D_DTL we demonstrate how this metric performs.

[![dtl_1d]({{images}}/dtl_1d.png)]({{images}}/dtl_1d.png){:.center-img .width-70} 

Figure S_1D_DTL is generated with the script `ZZ_SCRIPTS/figure_s_1d_DTL.R`

<details markdown="1"><summary>Show/Hide Script</summary>
```R
#! /usr/bin/env Rscript

request_scvs <- TRUE
request_regs <- FALSE
withRestarts(source("load_data.R"), terminate=function() message('load_data.R: data already loaded. Nice.'))

library(tidyverse)
library(cowplot)

args <- list()
args$output <- "../YY_PLOTS/FIG_S_1D_DTL"

# Create directory
dir.create(args$output, showWarnings=F, recursive=T)

# -----------------------------------------------------------------------------
# per-site D histograms
# -----------------------------------------------------------------------------

get_D_hist <- function(variable, type, col="#9D082D", col_shuffle="#333333", xlim=NA, ylim=NA,
                       sqrty=FALSE, replicates=10, shuffle=T, shuffle_trace=F) {
    plot_data <- scvs %>%
        select((!!sym(variable)), (!!sym(type))) %>%
        filter(!is.na((!!sym(type))), !is.na((!!sym(variable))))
    count_method <- "density"
    g <- ggplot()
    binwidth <- 1
    type_frac <- c()
    type_var <- c()
    for (i in 1:replicates) {
        shuffle_plot_data <- plot_data %>%
            sample_n(size=plot_data %>% dim() %>% .[[1]], replace=F)
        type_frac <- append(type_frac, shuffle_plot_data[,type] %>% .[[1]])
        type_var <- append(type_var, plot_data[,variable] %>% .[[1]])
        plot_data$shuffled <- shuffle_plot_data %>% .[[2]]
        if (shuffle_trace) {
            g <- g +
                stat_bin(
                    data = plot_data,
                    mapping=aes_string(x=variable, weight="shuffled", y=paste("..", count_method, "..", sep="")),
                    geom="step",
                    origin=0,
                    color=col_shuffle,
                    size=0.3,
                    alpha=0.4,
                    bins=85
                )
        }
    }
    mean_shuffle <- data.frame(type_frac, type_var)
    g <- g +
        geom_histogram(
            data = mean_shuffle,
            mapping=aes_string(x="type_var", weight="type_frac", y=paste("..", count_method, "..", sep="")),
            origin=0,
            fill=col_shuffle,
            alpha=0.7,
            binwidth=binwidth
        )
    g <- g +
        stat_bin(
            data = plot_data,
            mapping=aes_string(x=variable, weight=type, y=paste("..", count_method, "..", sep="")),
            geom="step",
            center=0.0,
            color=col,
            size=1.0,
            binwidth=binwidth
        )
    g <- g +
        labs(y=count_method, x=variable) +
        theme_classic() +
        theme(
            text=element_text(size=11, family="Helvetica", face="bold")
        ) +
        scale_x_continuous(limits = xlim, expand = c(0, 0))
    if (sqrty) {
        g <- g + scale_y_continuous(limits = c(NA, ylim), trans='sqrt', expand = c(0.005, 0))
    } else {
        g <- g + scale_y_continuous(limits = c(NA, ylim), expand = c(0.005, 0))
    }
    g
}

N <- 10
pn_D <- get_D_hist(variable='ANY_D', type='pN_popular_consensus', col=ns_col, col_shuffle=ns_col, xlim=c(0,350), ylim=0.045, sqrty=F, replicates=N) +
    labs(y="Density", x="1D distance-to-ligand (# residues)")
ps_D <- get_D_hist(variable='ANY_D', type='pS_popular_consensus', col=s_col, col_shuffle=s_col, xlim=c(0,350), ylim=0.045, sqrty=F, replicates=N) +
    labs(y="Density", x="1D distance-to-ligand (# residues)")

pn_D_sub <- pn_D +
    scale_x_continuous(limits = c(0,20), expand = c(0, 0)) +
    scale_y_continuous(limits = c(NA, NA), expand = c(0.005, 0)) +
    labs(y="", x="") +
    theme(text=element_text(size=10, family="Helvetica", face="bold"))
ps_D_sub <- ps_D +
    scale_x_continuous(limits = c(0,20), expand = c(0, 0)) +
    scale_y_continuous(limits = c(NA, NA), expand = c(0.005, 0)) +
    labs(y="", x="") +
    theme(text=element_text(size=10, family="Helvetica", face="bold"))

# -----------------------------------------------------------------------------
# Put it all together
# -----------------------------------------------------------------------------

plot1 <- ggdraw(pn_D) +
    draw_plot(pn_D_sub, 0.45, 0.45, 0.5, 0.5)
plot2 <- ggdraw(ps_D) +
    draw_plot(ps_D_sub, 0.45, 0.45, 0.5, 0.5)

display(
    plot_grid(
        plot1, plot2, ncol=2, align='v'),
    output = file.path(args$output, 'fig.pdf'),
    w = 7.5,
    h = 2.5
)
display(
    plot_grid(
        plot1, plot2, ncol=2, align='v'),
    output = file.path(args$output, 'fig.png'),
    w = 7.5,
    h = 2.5
)

```
</details> 

which can be ran like so:

<div class="extra-info" style="{{ command_style  }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
cd ZZ_SCRIPTS/
Rscript figure_s_1d_DTL.R
cd -
```
‣ **Time:** FIXME  
‣ **Storage:** Minimal  
</div> 

## Analysis X: Comparison to BioLiP and DTL cutoff

### DTL cutoff

```bash
mkdir 20_BIOLIP
cd 20_BIOLIP
wget http://zhanglab.ccmb.med.umich.edu/BioLiP/download/BioLiP.tar.bz2
tar -zxvf BioLiP.tar.bz2
cd -
python ZZ_SCRIPTS/biolip_dtl_dist.py
cd ZZ_SCRIPTS/
Rscript figure_s_biolip.R
cd -
# Output is YY_PLOTS/FIG_S_BIOLIP
```

## Analysis X: Big linear models

<div class="extra-info" style="{{ command_style  }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
cd ZZ_SCRIPTS
Rscript analysis_big_linear_models.R
cd ..
```
‣ **Time:** 8 hours  
‣ **Storage:** 10 Gb  
‣ **Memory:** Minimal  
</div> 

Creates `table_models.txt`.

## Analysis X: Gene-sample paired models

In Figure 1 we show distributions of Pearson coefficients thousands of different models. What are these models and how were they created? Let's dive into that.

[![fig1cd]({{images}}/fig1cd.png)]({{images}}/fig1cd.png){:.center-img .width-70} 

## Analysis X: Per-group models

<div class="extra-info" style="{{ command_style  }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
for bin_number in 2 4 6 8 10 14 18 22 26 30 34 38; do
    python ZZ_SCRIPTS/analysis_pnps_d_and_rsa.py -b $bin_number \
                                                 -o 17_PNPS_RSA_AND_DTL
done
```
‣ **Time:** 30 min per bin  
‣ **Storage:** Minimal  
</div> 

## Analysis X: dN/dS$^{(gene)}$ between HIMB83 and HIMB122

<div class="extra-info" markdown="1">
<span class="extra-info-header">Step X Info</span>
‣ **Prerequisite steps:** Aux. Step X  
</div> 

This section requires you to first complete Auxiliary Step X, please complete that before continuing.

In this section we calculate the dN/dS of homologous genes shared between HIMB83 and HIMB122, a related SAR11 genome. In the paper, the rationale for doing this was to validate our approach of pN/pS$^{(gene)}$ by averaging across samples and comparing the sample-averaged pN/pS$^{(gene)}$ values to dN/dS$^{(gene)}$ between HIMB83 and HIMB122, with the expectation that these should be qualitatively similar to one another given the evolutionary relatedness of HIMB122 to HIMB83.

## Analysis X: Glutamine synthetase (GS)

<div class="extra-info" style="{{ analysis_style  }}" markdown="1">
<span class="extra-info-header">Analysis X Info</span>
‣ **Prerequisite steps/analyses:** None  
‣ **Checkpoint datapack:** None  
</div> 

### Dodecameric RSA & DTL

In the study we focus on glutamine synthetase for a case study, and when dealing with its structure we make the following point in the text:

<blockquote>
Since the native quaternary structure of GS is a dodecameric complex (12 monomers), our monomeric estimates of RSA and DTL are unrepresentative of the active state of GS. We addressed this by aligning 12 copies of the predicted structure to a solved dodecameric complex of GS in Salmonella typhimurium (PDB ID 1FPY), which HIMB83 GS shares 61% amino acid similarity with (Figure 3a). From this stitched quaternary structure we recalculated RSA and DTL, and as expected, this yielded lower average RSA and DTL estimates due to the presence of adjacent monomers (0.17 versus 0.24 for RSA and 17.8Å versus 21.2Å for DTL).
<div class="blockquote-author">
  <b>Kiefl et al. 2022, pre-print</b>
</div>
</blockquote>

FIXME

<div class="extra-info" style="{{ command_style  }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
mkdir -p 21_GS_COMPLEX
pymol -c ZZ_SCRIPTS/GS_neighor_complex.pml
python ZZ_SCRIPTS/GS_complex_DTL.py
python ZZ_SCRIPTS/GS_complex_RSA.py
```
‣ **Time:** <1 min  
‣ **Storage:** 4 Mb  
</div> 

### How similar to HIMB122 to HIMB83?

First, let's determine how similar HIMB122 is to HIMB83. Here is a genome content comparison between HIMB83 and HIMB122

[![himb122_pan]({{images}}/himb122_pan.png)]({{images}}/himb122_pan.png){:.center-img .width-70}

You are a few seconds away from playing with this data yourself:

```bash
anvi-display-pan -p 07_PANGENOME_COMP_TO_HIMB122/PANGENOME/SAR11-PAN.db -g 07_PANGENOME_COMP_TO_HIMB122/SAR11-GENOMES.db
```

Clearly there is a lot of overlapping gene content. We consider HIMB122 to share a homolog with HIMB83 based on the clustering results of this pangenome. If two genes are clustered into the same gene cluster, they are considered homologs. For example, here is a randomly chosen gene cluster, where we can see the corresponding protein sequences of the two genes:

[![himb122_gc]({{images}}/himb122_gc.png)]({{images}}/himb122_gc.png){:.center-img .width-100} 
On a nucleotide level, we can assess the overall similarity between HIMB122 and HIMB83 using average nucleotide identity (ANI) calculated by pyANI using {% include PROGRAM name="anvi-compute-genome-similarity" %}.

<div class="extra-info" style="{{ command_style  }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
anvi-compute-genome-similarity -e 07_EXTERNAL_GENOMES_COMP_TO_HIMB122.txt \
                               -p 07_PANGENOME_COMP_TO_HIMB122/PANGENOME/SAR11-PAN.db \
                               --program pyANI \
                               --num-threads 6 \
                               -o 07_ANI_HIMB122
```
‣ **Time:** <1 min  
‣ **Storage:** Minimal  
</div> 

Looking at the output `07_ANI_HIMB122/ANIb_percentage_identity.txt`, we can see that the genome ANI of HIMB122 to HIMB83 is **82.6%**. Note that this only considers aligned segments of the genome. If the full sequences are considered this value drops to 63.7% similarity (`07_ANI_HIMB122/ANIb_full_percentage_identity.txt`), however since we are interested in calculating dN/dS$^{(gene)}$ between homologs, it is the former metric that is more relevant.

### Calculating dN/dS$^{(gene)}$ for 1 gene

To calculate dN/dS$^{(gene)}$, I opted to use PAML's `yn00` program, which utilizes the [Yang and Nielson (2000)](https://pubmed.ncbi.nlm.nih.gov/10666704/) counting-based method for dN/dS. While dN/dS is a parameter that is estimated during the process of maximum-likelihood phylogenetic tree estimation, and while such estimates are typically considered higher accuracy and offer greater flexibility than counting-based methods since one can calculate branch-specific dN/dS values, in this simple pairwise case between HIMB83 and HIMB122 I did not see the necessity of such complexity.

First I will walk through the steps required for calculating dN/dS using `yn00` for one homologous pair, and then I present the script that calculates them for all homologs shared between HIMB83 and HIMB122.

For the example I will use gene cluster `GC_00000017`, which has the HIMB83 gene with ID `2003`. This happens to be ribosomal protein L35 (PF01632), which I found out by opening up `functions.txt` and searching for gene `1324`.

As a first step, I need the nucleotide and amino acid sequences for each homolog. For this, I use the program {% include PROGRAM name="anvi-get-sequences-for-gene-clusters" %}:

```bash
anvi-get-sequences-for-gene-clusters -p 07_PANGENOME_COMP_TO_HIMB122/PANGENOME/SAR11-PAN.db \
                                     -g 07_PANGENOME_COMP_TO_HIMB122/SAR11-GENOMES.db \
                                     --gene-cluster-id GC_00000017 \
                                     --max-num-genes-from-each-genome 1 \
                                     -o GC_00000017.aln.faa
anvi-get-sequences-for-gene-clusters -p 07_PANGENOME_COMP_TO_HIMB122/PANGENOME/SAR11-PAN.db \
                                     -g 07_PANGENOME_COMP_TO_HIMB122/SAR11-GENOMES.db \
                                     --gene-cluster-id GC_00000017 \
                                     --max-num-genes-from-each-genome 1 \
                                     -o GC_00000017.aln.fna \
                                     --report-DNA
```

The FASTA deflines of these files look something like this

```
>00000001|gene_cluster:GC_00000017|genome_name:HIMB122|gene_callers_id:4577
ACTGC....
>00000000|gene_cluster:GC_00000017|genome_name:HIMB083|gene_callers_id:2003
ACTGC....
```

which PAML insists on truncating without warning, so we have to modify the deflines with the script `ZZ_SCRIPTS/rename_for_paml.py`.

<details markdown="1"><summary>Show/Hide rename_for_paml.py</summary>
```python
#! /usr/bin/env python

import argparse
import shutil
import anvio.fastalib as u

ap = argparse.ArgumentParser()
ap.add_argument("--fasta", "-f", required=True)
args = ap.parse_args()

output = u.FastaOutput(args.fasta + ".temp")
fasta = u.SequenceSource(args.fasta)

while next(fasta):
    defline = fasta.id.split('|')[2].split(':')[-1]
    output.write_id(defline)
    output.write_seq(fasta.seq, split = False)

fasta.close()
output.close()

shutil.move(args.fasta + '.temp', args.fasta)
```
</details> 

```bash
python ZZ_SCRIPTS/rename_for_paml.py -f GC_00000017.fna
python ZZ_SCRIPTS/rename_for_paml.py -f GC_00000017.aln.faa
```

You probably have noticed that the amino acid FASTA is labeled with `aln` in its name, that's because {% include PROGRAM name="anvi-get-sequences-for-gene-clusters" %} reports amino acid alignments of the protein sequence (_i.e._ there potentially exists gap characters so that the sequences align even if they are of different lengths). However, this alignment does not yet exist for the nucleotide sequence. To create this alignment, we use a program `ZZ_SCRIPTS/pal2nal.pl`, which has been floating around the internet since this [2006 paper](https://academic.oup.com/nar/article/34/suppl_2/W609/2505720) that converts protein alignments to codon alignments.

Using this ancient relic, we can create the codon alignment `GC_00000017.aln.fna` that PAML (`yn00`) so desires.

```
   2    204
HIMB083
ATGCCAAAACTAAAAACAAAAAGCTCAGCAAAAAAAAGATTTAAAATATCTGCTAAGGGT
AAAGTTATTTCTGCTCAGGCAGGTAAGAGACATGGTATGATTAAAAGAACGAATTCACAA
ATAAGAAAACTAAGAGGCACCACAACATTGTCTAAACAAGATGGAAAAATTGTTAAGTCA
TACATGCCTTACAGTTTAAGAGGT
HIMB122
ATGCCAAAATTAAAAACAAAAAGCTCTGCAAAAAAAAGATTTAAAATATCAGCTAAGGGT
AAAGTTATATCTGCTCAAGCAGGTAAAAGACACGGTATGATTAAACGAACGAATTCACAA
ATTAGAAAATTAAGAGGCACTACAACTTTGTCTAAACAAGATGGCAAAATTGTTAAATCA
TACATGCCTTACAGTTTAAGAGGG
```

Nevertheless, the final thing we need to do is make a control file that specifies how `yn00` should be run. Basically, this is where you provide inputs, outputs, and any parameters to tune `yn00`'s behavior. This is the vanilla file we want:

```
seqfile = GC_00000017.aln.fna * sequence data file name
outfile = GC_00000017_output * main result file
verbose = 0 * 1: detailed output (list sequences), 0: concise output
icode = 0 * 0:universal code; 1:mammalian mt; 2-10:see below
weighting = 0 * weighting pathways between codons (0/1)?
commonf3x4 = 0 * use one set of codon freqs for all pairs (0/1)?
```

Create this file, name it `GC_00000017.ctl` and then run `yn00` with `yn00 GC_00000017.ctl`. Things should finish nearly instantly, and the newly created `GC_00000017_output` file should now exist in your directory. It is a difficult file to parse, however the dN/dS value we are after can be parsed with the following command.

```bash
echo $(grep 'omega' GC_00000017_output -A 2 | tail -n 1 | awk '{print $7}')
```

Which yields the answer `0.00000`. Ok so this was a bit underwhelming, because in retrospect these homologs have 0 amino acid substitutions relative to one another, so by definition dN = 0. Regardless, hopefully this has illustrated the workflow that is about to occur.

If you followed this exercise, make sure you clean up all of the files you made with the following commands:

```bash
rm GC_00000017.aln.faa GC_00000017.aln.fna GC_00000017.ctl GC_00000017.fna GC_00000017_output
```

### Calculating dN/dS$^{(gene)}$ for all

That is quite a bit of tedium to calculate dN/dS for one homologous pair, so I wrote a script that can calculate dN/dS for all of them called `ZZ_SCRIPTS/calculate_dnds.sh`

<details markdown="1"><summary>Show/Hide Script</summary>
```bash
#! /usr/bin/env bash

GOI_GCOI=$1
PAN=$2
STORAGE=$3
OUTPUT=$4

rm -rf $OUTPUT
mkdir $OUTPUT
cd $OUTPUT

cat ../$GOI_GCOI | while read g gcoi; do

    mkdir -p $g
    cd $g

    echo $gcoi

    # Get the protein alignment
    anvi-get-sequences-for-gene-clusters -p ../../$PAN -g ../../$STORAGE -o $gcoi.aln.faa --gene-cluster-id $gcoi --max-num-genes-from-each-genome 1

    FILE=$gcoi.aln.faa
    if [ -f $FILE ]; then
        # Get the unaligned nucleotie sequences
        anvi-get-sequences-for-gene-clusters -p ../../$PAN -g ../../$STORAGE -o $gcoi.fna --gene-cluster-id $gcoi --max-num-genes-from-each-genome 1 --report-DNA

        # Rename the FASTA deflines so PAML likes them
        python ../../ZZ_SCRIPTS/rename_for_paml.py -f $gcoi.fna
        python ../../ZZ_SCRIPTS/rename_for_paml.py -f $gcoi.aln.faa

        # Create a codon alignment from the protein alignment
        ../../ZZ_SCRIPTS/pal2nal.pl $gcoi.aln.faa $gcoi.fna -output paml -nogap > $gcoi.aln.fna

        # Create the control file that yn00 wants
        python ../../ZZ_SCRIPTS/_gen_yn00_ctl_file.py -a $gcoi.aln.fna --codeml-output ${gcoi}_output --output $gcoi.ctl

        # Run yn00
        yn00 $gcoi.ctl

        # Parse the yn00 output to yield dN/dS (omega)
        omega=$(grep 'omega' ${gcoi}_output -A 2 | tail -n 1 | awk '{print $7}')
    else
        omega=""
    fi

    echo -e "$g\t$omega" >> ../dnds.txt
    cd ..

done
```
</details>

This script does exactly what we just did manually, except it wraps it up into a loop thats iterated for each homologous pair. It takes as input a file you generated earlier in this analysis, `goi_gcoi_COMP_TO_HIMB122`, which corresponds each gene cluster to the gene ID of HIMB83. The script should be ran like so:

<div class="extra-info" style="{{ command_style  }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
bash ZZ_SCRIPTS/calculate_dnds.sh goi_gcoi_COMP_TO_HIMB122 \
                                  07_PANGENOME_COMP_TO_HIMB122/PANGENOME/SAR11-PAN.db \
                                  07_PANGENOME_COMP_TO_HIMB122/SAR11-GENOMES.db \
                                  19_DNDS_HIMB122
```
‣ **Time:** 3 hours  
‣ **Storage:** 30 Mb  
</div> 

For many genes, you may see the following error message coming from anvi'o:

```
Config Error: Bad news: the combination of your filters resulted in zero gene clusters :/
              These are the filtesr anvi'o used: --min-num-genomes-gene-cluster-occurs 0,
              --max-num-genomes-gene-cluster-occurs 2, --min-num-genes-from-each-genome 0,
              --max-num-genes-from-each-genome 1, --min-functional-homogeneity-index
              -1.000000, --max-functional-homogeneity-index 1.000000, --min-geometric-
              homogeneity-index -1.000000, --max-geometric-homogeneity-index 1.000000, --min-
              combined-homogeneity-index -1.000000, and --max-combined-homogeneity-index
              1.000000. None of your 1 gene clusters in your 2 genomes that were included this
              analysis matched to this combination (please note that number of genomes may be
              smaller than the actual number of genomes in the original pan genome if other
              filters were applied to the gene clusters dictionary prior).
```

This is actually desired behavior. We have set the parameter `--max-num-genes-from-each-genome 1` in {% include PROGRAM name="anvi-get-sequences-for-gene-clusters" %} in order to prevent circumstances where one gene cluster had multiple genes from a single genome. In such a scenario, we opt to ignore the cluster altogether, since it is ambiguous which homolog should be used for the comparison.

This will take some time to run, but when it has finished the directory `19_DNDS_HIMB122` will be populated with all of the alignment information of each homologous comparison. But most important is the file `19_DNDS_HIMB122/dnds.txt`, which holds all of the dN/dS$^{(gene)}$ for each HIMB83 gene with respect to the homologous HIMB122 gene.

### Visualizing

To summarize the results, I created a scatter plot between sample-averaged pN/pS and dN/dS which ended up being Figure S_DNDS.

<details markdown="1"><summary>Show/Hide Script</summary>
```R
#! /usr/bin/env Rscript

request_scvs <- FALSE
request_regs <- FALSE
withRestarts(source("load_data.R"), terminate=function() message('load_data.R: data already loaded. Nice.'))

library(latex2exp)
library(optparse)
library(tidyverse)
library(cowplot)

args <- list()
args$output <- "../YY_PLOTS/FIG_S_DNDS"

# Create directory
dir.create(args$output, showWarnings=F, recursive=T)

# -----------------------------------------------------------------------------
# Load and merge dN/dS data with pN/pS
# -----------------------------------------------------------------------------

# pN/pS exists, but needs to be sample-averaged
pnps_sample_averaged <- pnps %>%
    group_by(gene_callers_id) %>%
    summarize(pnps = mean(pnps, na.rm=T))

# Load up dN/dS
dnds <- read_tsv("../19_DNDS_HIMB122/dnds.txt", col_names=F) %>%
    rename(gene_callers_id=X1, dnds=X2) %>%
    filter(!is.na(dnds))

# Merge them together
df <- full_join(pnps_sample_averaged, dnds) %>%
    mutate(log_pnps = log10(pnps), log_dnds = log10(dnds))

# -----------------------------------------------------------------------------
# Make some plots
# -----------------------------------------------------------------------------

g1 <- ggplot(data=df) +
    geom_histogram(aes(dnds)) +
    scale_y_continuous(expand=c(0,0))
display(g1, file.path(args$output, "dnds_histogram.png"), width=3.2, height=2.8)

overlap <- df %>% filter(!is.na(pnps), !is.na(dnds))
R2 <- overlap %>% lm(pnps ~ dnds, data=.) %>% summary() %>% .$r.squared %>% round(2)
g2 <- ggplot(data=overlap, aes(log_dnds, log_pnps)) +
    geom_abline(slope=1, size=1.25, alpha=0.7) +
    geom_text_repel(aes(label=gene_callers_id), alpha=0.8, nudge_x=0.4) +
    annotate('text', label=paste("R² =", R2), x=-3.7, y=-0.7, size=3, bold=T) +
    my_theme(8) +
    labs(x=TeX("$\\log_{10}$ dN/dS$^{(gene)}$", bold=T), y=TeX("$\\log_{10}$ sample-averaged pN/pS$^{(gene)}$", bold=T))
display(g2, file.path(args$output, "dnds_vs_pnps.png"), width=3.2, height=2.8)
```
</details> 

Running the following yields Figure S_DNDS under the filename `YY_PLOTS/FIG_S_DNDS/dnds_vs_pnps.png`

<div class="extra-info" style="{{ command_style  }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
cd ZZ_SCRIPTS/
Rscript figure_s_dnds.R
cd ..
```
‣ **Time:** <1 min  
‣ **Storage:** Minimal  
</div> 

[![dnds]({{images}}/dnds.png)]({{images}}/dnds.png){:.center-img .width-70} 

Since it is much more common for slightly negative polymorphisms to drift to observable frequencies than it is for them to fixate, it is expected that in the majority of cases sample-averaged pN/pS exceeds dN/dS, and it is indeed what we see (most genes are above the black line $y = x$). Interestingly, there are a select number of genes that have quite high rates of polymorphim, despite dN/dS being very low. I think there is probably an interesting story involving all of the genes with high polymorphism rates but low substitution rates (left side of plot).

## Aux. Step X: Pangenome detour

<div class="extra-info" markdown="1">
<span class="extra-info-header">Step X Info</span>
‣ **Prerequisite steps:** None  
‣ **Checkpoint datapack:** None  
‣ **Central:** No  
</div> 

This is the first **auxiliary** step, meaning it is not required for the central analyses of this paper. For this reason, you can skip this step if you want.

But if you're reading this, you've more than likely beeen redirected to this section of the workflow, because the analysis you're interested in requires this step to be completed. In that case, you're in the right place.

How does the gene content differ between HIMB83 and the other 20 genomes? What genes are present and absent, what is the percent similarity of the homologs? What is the average nucleotide identity (ANI) of these genomes to HIMB83?

These questions all fall under the umbrella of [pangenomics](https://merenlab.org/momics/#pangenomics), and since we now have {% include ARTIFACT name="contigs-db" text="contigs-dbs" %} for each genome (in `07_SPLIT/`), we are in the perfect position for a quick detour into the land of pangenomics. While a detour, this is required for some supplemental figures down the road.

In this paper I made 2 pangenomes. One compares all 21 SAR11 genomes to one another, and the other compares HIMB83 to its closest relative in this genome collection, HIMB122.

To create this 2 pangenomes, I created a bash script called `ZZ_SCRIPTS/make_pangenomes.sh`.

<details markdown="1"><summary>Show/Hide Script</summary>
```bash
#! /usr/bin/env bash

# Delete all previous pangenomes
rm -rf 07_PANGENOME*
rm -rf 07_SUMMARY_PAN*

mkdir 07_PANGENOME
mkdir 07_PANGENOME_COMP_TO_HIMB122

python ZZ_SCRIPTS/gen_external_genomes_file.py

# Full pangeome
anvi-gen-genomes-storage -e 07_EXTERNAL_GENOMES.txt -o 07_PANGENOME/SAR11-GENOMES.db
anvi-pan-genome -g 07_PANGENOME/SAR11-GENOMES.db -n SAR11 -o 07_PANGENOME/PANGENOME -T $1
anvi-script-add-default-collection -p 07_PANGENOME/PANGENOME/SAR11-PAN.db
anvi-summarize -p 07_PANGENOME/PANGENOME/SAR11-PAN.db -g 07_PANGENOME/SAR11-GENOMES.db -C DEFAULT -o 07_SUMMARY_PAN
gzip -d 07_SUMMARY_PAN/SAR11_gene_clusters_summary.txt.gz

# HIMB122
anvi-gen-genomes-storage -e 07_EXTERNAL_GENOMES_COMP_TO_HIMB122.txt -o 07_PANGENOME_COMP_TO_HIMB122/SAR11-GENOMES.db
anvi-pan-genome -g 07_PANGENOME_COMP_TO_HIMB122/SAR11-GENOMES.db -n SAR11 -o 07_PANGENOME_COMP_TO_HIMB122/PANGENOME -T $1
anvi-script-add-default-collection -p 07_PANGENOME_COMP_TO_HIMB122/PANGENOME/SAR11-PAN.db
anvi-summarize -p 07_PANGENOME_COMP_TO_HIMB122/PANGENOME/SAR11-PAN.db -g 07_PANGENOME_COMP_TO_HIMB122/SAR11-GENOMES.db -C DEFAULT -o 07_SUMMARY_PAN_COMP_TO_HIMB122
gzip -d 07_SUMMARY_PAN_COMP_TO_HIMB122/SAR11_gene_clusters_summary.txt.gz
```
</details> 

First, a pair of {% include ARTIFACT name="external-genomes" %} files are made: `07_EXTERNAL_GENOMES.txt` and `07_EXTERNAL_GENOMES_COMP_TO_HIMB122.txt`. These are used by the program {% include PROGRAM name="anvi-gen-genomes-storage" %} to create {% include ARTIFACT name="genomes-storage-db" text="genomes-storage-dbs" %}, which are in turn used by the program {% include PROGRAM name="anvi-pan-genome" %} to create {% include ARTIFACT name="pan-db" text="pan-dbs" %}. The data within these databases are summarized with the program {% include PROGRAM name="anvi-summarize" %}, which creates two output directories: `07_SUMMARY_PAN_COMP_TO_HIMB122` and `07_SUMMARY_PAN`.

When you're ready, run it:

<div class="extra-info" style="{{ command_style  }}" markdown="1">
<span class="extra-info-header">Command #X</span>
```bash
bash ZZ_SCRIPTS/make_pangenomes.sh <NUM_THREADS>
```
‣ **Time:** ~(25/`<NUM_THREADS>`) min  
‣ **Storage:** 259 Mb  
</div> 

Once its finished, feel free to explore whatever you want. For example, the `07_SUMMARY_PAN*` directories have tabular data you can open in Excel, as well as HTML files that can be opened and explored in your browser. You could also load up the interactive interface and explore this data interactively with {% include PROGRAM name="anvi-display-pan" %}:

```
anvi-display-pan -p 07_PANGENOME/PANGENOME/SAR11-PAN.db \
                 -g 07_PANGENOME/SAR11-GENOMES.db
```

Which would present you with a plot like this, which you can interactively explore:

[![pan_1]({{images}}/pan_1.png)]( {{images}}/pan_1.png){:.center-img .width-70}

However, for the purposes of this study, we are done here. These pangenomes will be useful for creating some supplemental figures you're probably already aware of.

## Reproducing numbers in the text

FIXME Introduce

----------------------------

<blockquote>
(...) quote goes here (...)
</blockquote>

```R
code goes here
```

----------------------------


<blockquote>
(...) quote goes here (...)
</blockquote>

```R
code goes here
```

----------------------------





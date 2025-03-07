---
layout: page
title: A reproducible workflow for Meyer & Eren manuscript 2024
modified: 2024-08-30
excerpt: "A bioinformatics workflow for our metagenomics data integration across observatories"
comments: true
authors: [raissa]
---

{% capture images %}{{site.url}}/data/micr-obs-future/images{% endcapture %}


<div class="extra-info" markdown="1">

<span class="extra-info-header">Summary</span>

**The purpose of this page** is to provide access to reproducible data products and analyses for the manuscript "**Marine Microbial Observatories for the Future: From Samples to Data to Legacy using integrated omics strategies**" by Meyer & Eren


In this study, we 
-  integrated metagenomics metadata and data from observatories (Hawaii Ocean Time-Series and Bermuda Atlantic Time-series Study), sampling expeditions (Bio-GO-SHIP, bioGEOTRACES, Malaspina, and Tara Oceans), and citizen science initiatives (Ocean Sampling Day)
-  generated an anvi'o contigs database describing 51 SAR11 isolate genomes (reference genomes)
-  competitively recruited reads from the above-listed projects' metagenomes to the SAR11 reference genomes, and profile the recruitment results
-  investigated the patterns in genes across metagenomes recruited to the individual SAR11 reference genomes

Sections in this document will detail all the steps of downloading and processing SAR11 genomes and metagenomes, mapping metagenomic reads onto the SAR11 genomes, as well as analyzing and visualizing the outcomes.

For the curation of metadata, please consult the [public-marine-omics-metadata](https://github.com/merenlab/public-marine-omics-metadata/tree/main) GitHub repository, where all steps of metadata gathering, curation, and standardization are described in detail.

{:.notice}
If you have any questions, notice an issue, and/or are unable to find an important piece of information here, please feel free to leave a comment down below, send an e-mail to [us](/people/), or get in touch with us through Discord:

{% include _join-anvio-discord.html %}

In the following table, we are giving a summary of the samples and their projects included in this analysis (information for depths up to 100 m). In total, we started this analysis with 1825 samples, so buckle up! 💺

| Project | Acronym | Accession | Years sampled | # Samples |
|:--|:--|:--|:--|:--|
| [Bermuda Atlantic Time-series Study](https://bats.bios.asu.edu/about/) | BATS | PRJNA385855 | 2003 - 2004, 2009 | 40 |
| [bioGEOTRACES](https://www.nature.com/articles/sdata2018176) | BGT | PRJNA385854 | 2010, 2011 | 323 |
| [Bio-GO-SHIP](https://biogoship.org) | BGS | PRJNA656268 | 2011, 2013, 2014, 2016 - 2018 2020 | 969 |
| [Hawaii Ocean Time-series](http://hahana.soest.hawaii.edu/hot/hot_jgofs.html) | HOT1 \| HOT3 |  PRJNA385855 \| PRJNA352737 | 2003, 2004 \| 2014 - 2017 | 28 \| 230 |
| [Malaspina](https://www.nature.com/articles/s41597-024-02974-1) | MAL | PRJEB52452 | 2011 | 16 |
| [Ocean Sampling Day 2014](https://doi.org/10.1186/s13742-015-0066-5) | OSD | PRJEB8682 | 2014 | 127 |
| [Tara Oceans](https://fondationtaraocean.org/en/home/) | TARA | PRJEB1787 | 2009 - 2012 | 92 |


</div>

# General notices

{:.notice}
All anvi'o analyses in this document are performed using the anvi'o developer version during the era of `v8`. Please see [the installation notes](https://anvio.org/install/) to download the appropriate version through PyPI, Docker, or GitHub.


To check your version, you can use the `anvi-help -v` command. Here is mine:
```bash
$ anvi-help -v
Anvi'o .......................................: marie (v8-dev)
Python .......................................: 3.10.13

Profile database .............................: 40
Contigs database .............................: 24
Pan database .................................: 20
Genome data storage ..........................: 7
Auxiliary data storage .......................: 2
Structure database ...........................: 2
Metabolic modules database ...................: 4
tRNA-seq database ............................: 2
```

{:.notice}
You will see that many of the commands utilize the SLURM wrapper `clusterize` to submit jobs to our HPC clusters. `Clusterize` was developed by Evan Kiefl and is described here: https://github.com/ekiefl/clusterize. Thank you, Evan! 

## Metadata
For the curation of metadata, please consult the [public-marine-omics-metadata](https://github.com/merenlab/public-marine-omics-metadata/tree/main) GitHub repository, where all steps of metadata gathering, curation, and standardization are described in detail.

Remember, team: Data without metadata is like a mystery dish - sure, it might look good, but you don't know what's in it, you don't know if you can and are allowed to eat it, you don't know where it came from, you don't know who made it, you don't even really know what it is, ... 😖. Our data is often only worth as much as our metadata can support it.

## SAR11 cultivar genomes
This section explains how to prepare the set of 99 SAR11 isolate genomes, ending up with 51 quality-checked and dereplicated reference genomes. 

The isolate genomes available at the time of this analysis are included in a file called `SAR11_June2024_bycontig.fa`. We cannot make the file itself public yet, as it includes some unpublished isolate genomes (Thank you, Freel et al. [in prep] for sharing those with us!). To see more information and the source of each isolate genome, expand the section below. We would like to thank all researchers who have provided these genomes.

---

<details markdown="1"><summary>Click here to show/hide an overview of the reference genomes used</summary>

| Genome name                         | Genus     | Species        | Type | Study                        | Source.Location             |
|-------------------------------------|-----------|----------------|------|------------------------------|-----------------------------|
| HIMB1321                            | Ia.3.III  | Ia.3.III       | No   | Brandon 2006                 | coastal, Oahu, Hawaii, USA  |
| HIMB122                             | Ia.3.VI   | Ia.3.VI        | No   | Brandon 2006                 | coastal, Oahu, Hawaii, USA  |
| HIMB140                             | Ia.3.VI   | Ia.3.VI        | No   | Brandon 2006                 | coastal, Oahu, Hawaii, USA  |
| HIMB5                               | Ia.3.II   | Ia.3.II        | Type | Grote et al., 2012           | coastal, Oahu, Hawaii, USA  |
| HIMB4                               | Ia.3.III  | Ia.3.III       | Type | Grote et al., 2012           | coastal, Oahu, Hawaii, USA  |
| HIMB83 (might be listed as HIMB083) | Ia.3.V    | Ia.3.V         | No   | Grote et al., 2012           | coastal, Oahu, Hawaii, USA  |
| RS39                                | Ia.4.RS39 | Ia.4.RS39      | Type | Jimenez-Infante et al., 2017 | central Red Sea             |
| RS40                                | Ib.1      | Ib.1           | Type | Jimenez-Infante et al., 2017 | central Red Sea             |
| NP1                                 | Ia.1.I    | P. giovannonii | Type | Morris et al., (2020)        | open ocean, NE Pacific      |
| HTCC9565                            | Ia.1.I    | Ia.1           | Type | NA                           | open ocean, NE Pacific      |
| HTCC9022                            | Ia.3.IV   | Ia.3.IV        | No   | NA                           | coastal Oregon, USA         |
| HTCC1002                            | Ia.1.I    | Ia.1.I         | No   | Rappe et al., 2002           | coastal Oregon, USA         |
| HTCC1013                            | Ia.1.I    | Ia.1.I         | No   | Rappe et al., 2002           | coastal Oregon, USA         |
| HTCC1016                            | Ia.1.I    | Ia.1.I         | No   | Rappe et al., 2002           | coastal Oregon, USA         |
| HTCC1040                            | Ia.1.I    | Ia.1.I         | No   | Rappe et al., 2002           | coastal Oregon, USA         |
| HTCC1062                            | Ia.1.I    | Ia.1.I         | Type | Rappe et al., 2002           | coastal Oregon, USA         |
| HTCC7211                            | Ia.3.I    | Ia.3.I         | Type | Stingl et al., 2007          | Sargasso Sea, BATS          |
| HTCC7214                            | Ia.3.I    | Ia.3.I         | No   | Stingl et al., 2007          | Sargasso Sea, BATS          |
| HTCC7217                            | Ia.3.I    | Ia.3.I         | No   | Stingl et al., 2007          | Sargasso Sea, BATS          |
| HTCC8051                            | Ia.3.IV   | Ia.3.IV        | Type | Stingl et al., 2007          | coastal Oregon, USA         |
| FZCC0015                            | Ia.3.V    | Ia.3.V         | No   | Zhao et al., 2019            | coastal Pingtang, China     |
| HIMB123                             | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1420                            | -  | -     | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB1427                            | -   | -        | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB1520                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1564                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1412                            | -   | -   | Type | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB1505                            | -  | -       | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1527                            | -  | -       | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1549                            | -  | -       | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1409                            | -    | -         | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB1413                            | -   | -         | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB1444                            | -   | -         | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB1456                            | -    | -         | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB2201                            | -   | -         | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB2250                            | -    | -        | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB1430                            | -   | -        | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB1485                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1488                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1490                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1491                            | -  | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1493                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1494                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1495                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1506                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1507                            | -  | -       | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1509                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1513                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1518                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1521                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1526                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1542                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1552                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1556                            | -   | -        | Type | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1559                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1573                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1577                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1587                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1593                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1597                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1611                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1623                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1631                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1636                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1641                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1662                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1685                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1695                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1701                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1702                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1709                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1710                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1715                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1723                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1746                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1748                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1758                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1765                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1770                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1782                            | -   | -        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB2211                            | -   | -        | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB1483                            | - | -      | Type | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1402                            | -   | -        | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB1437                            | -   | -        | Type | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB1863                            | -   | -        | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB2104                            | -   | -        | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB2200                            | -   | -        | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB2204                            | -   | -        | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB2215                            |    | -        | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB1436                            | - | -      | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB2187                            | - | -      | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB2226                            | - | -      | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB2305                            | -      | -           | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB114                             | -      | -           | Type | Freel et al., in prep        |                             |
| HIMB1517                            | -      | -           | No   | Freel et al., in prep        |                             |
| HIMB1565                            | -      | -           | No   | Freel et al., in prep        |                             |
| IMCC9063                            | IIIa      | IIIa           | No   | Oh et al., 2011              |                             |
| LSUCC0530                           | IIIb      | IIIb           | Type | Henson et al., 2018          |                             |
| HIMB58                              | IIa       | IIa            | Type | Brandon M.S. thesis (2006)   |                             |

</details>

--- 

The first thing we did with these reference genomes in `SAR11_June2024_bycontig.fa`, is to separate them into individual .fa files - one per reference genome. We are doing this because we would like to assess the quality of each reference genome and dereplicate them in the next steps, and that requires them to be in separate files.

To perform that separation, we wrote the following script and ran it in the same directory that we have the `SAR11_June2024_bycontig.fa` file in (`fastaOriginal/`).
``` bash
nano separateFasta.py
```

content
``` python
from collections import defaultdict

def parse_fasta(file_path):
    """
    Parses the input FASTA file and organizes sequences by their ID.
    
    Args:
        file_path (str): Path to the input FASTA file.
        
    Returns:
        dict: A dictionary where the keys are sequence IDs and the values are lists of sequences.
    """
    with open(file_path, 'r') as file:
        sequences = defaultdict(list)  # Initialize a dictionary to hold sequences by their IDs
        header = None  # Variable to hold the current header
        for line in file:
            line = line.strip()  # Remove any leading/trailing whitespace
            if line.startswith('>'):  # Check if the line is a header (would start with ">")
                # Extract the sequence ID (part before the '-' character)
                header = line.split('-')[0][1:]
            else:
                # Append the sequence line to the corresponding sequence ID in the dictionary
                sequences[header].append(line)
    return sequences

def write_fasta(sequences):
    """
    Writes sequences to separate FASTA files based on their IDs.
    
    Args:
        sequences (dict): A dictionary where the keys are sequence IDs and the values are lists of sequences.
    """
    for seq_id, seq_list in sequences.items():
        # Create a file name based on the sequence ID
        file_name = f"{seq_id}.fa"
        with open(file_name, 'w') as file:
            for i, seq in enumerate(seq_list, start=1):
                # Write the header and sequence to the file
                file.write(f">{seq_id}-{i:09}\n")
                file.write(f"{seq}\n")

def main():
    # Path to the input FASTA file
    input_file = "SAR11_June2024_bycontig.fa"
    # Parse the FASTA file to get sequences organized by their IDs
    sequences = parse_fasta(input_file)
    # Write the sequences to separate FASTA files
    write_fasta(sequences)

# If this script is executed directly, run the main function
if __name__ == "__main__":
    main()

```

run
``` bash
python3 separateFasta.py
```

This resulted in 99 individual FASTA files, ready to be evaluated.

### Using checkM to evaluate the quality of isolate genomes

Before moving on with any of the reference genomes, we used `CheckM` to evaluate the completeness and contamination of the isolate genomes. These metrics are sometimes also given in the databases hosting the publicly available reference genomes (e.g., see [https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_900177485.1/](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_900177485.1/)), but it is always better to double-check. Besides, we need the completeness and contamination metrics for the unpublished reference genomes anyhow. 

Going into this, know that, since all these genomes are isolate genomes, we expect them to be of quite high completeness and low contamination (rightly so, as you'll see).

{:.notice}
There is also an option to include `checkM` in the `dRep` step that is following, however, that does not seem to work for everyone, so we are doing it separately.

<div class="extra-info" markdown="1">

<span class="extra-info-header">Citations and version</span>

CheckM citation:

- Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW. 2015. [CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes](http://genome.cshlp.org/content/25/7/1043.short). Genome Research, 25: 1043–1055.

CheckM relies on several other software packages:

- [pplacer](http://matsen.fhcrc.org/pplacer/): Matsen FA, Kodner RB, Armbrust EV. 2010. pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree. BMC Bioinformatics 11: doi:10.1186/1471-2105-11-538.
- [prodigal](http://prodigal.ornl.gov/): Hyatt D, Locascio PF, Hauser LJ, Uberbacher EC. 2012. Gene and translation initiation site prediction in metagenomic sequences. Bioinformatics 28: 2223–2230.
- [HMMER](http://hmmer.org/): http://hmmer.org/

We used CheckM version 1.2.2

</div>


We ran `checkM` on all files in the directory `fastaOriginal/` that end on `.fa` (so on the output files of the splitting we did above) and added the output into a directory called `check output/`. 

``` bash
clusterize -j checkM -o checkm.log -n 1 "checkm lineage_wf -t 40 -x fa ./fastaOriginal ./checkMoutput/ -f out_checkM.tab --tab_table"
```

{:.notice}
This was our first use of `clusterize`. It's that easy!

If you would like to inspect the output file, feel free to grab it [here](files/out_checkM.tsv). Otherwise, we have also included its content in the expandable section below.

---

<details markdown="1"><summary>Click here to show/hide the checkM output table</summary>


 
| Bin Id    | Marker lineage             | # genomes | # markers | # marker sets | 0  | 1   | 2  | 3 | 4 | 5 | Completeness | Contamination | Strain heterogeneity |
|-----------|----------------------------|-----------|-----------|---------------|----|-----|----|---|---|---|--------------|---------------|----------------------|
| FZCC0015  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 324 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |
| HIMB058   | k__Bacteria (UID2495)      | 2993      | 139       | 83            | 1  | 138 | 0  | 0 | 0 | 0 | 98.80        | 0.00          | 0.00                 |
| HIMB114   | k__Bacteria (UID2495)      | 2993      | 140       | 84            | 1  | 139 | 0  | 0 | 0 | 0 | 98.81        | 0.00          | 0.00                 |
| HIMB122   | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 322 | 2  | 0 | 0 | 0 | 100.00       | 0.95          | 0.00                 |
| HIMB123   | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 33 | 285 | 6  | 0 | 0 | 0 | 86.73        | 1.90          | 100.00               |
| HIMB1321  | o__Rickettsiales (UID3809) | 83        | 323       | 210           | 0  | 320 | 3  | 0 | 0 | 0 | 100.00       | 0.95          | 0.00                 |
| HIMB140   | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 1  | 323 | 0  | 0 | 0 | 0 | 99.53        | 0.00          | 0.00                 |
| HIMB1402  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 1  | 323 | 0  | 0 | 0 | 0 | 99.53        | 0.00          | 0.00                 |
| HIMB1409  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 324 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |
| HIMB1412  | o__Rickettsiales (UID3809) | 83        | 323       | 210           | 0  | 323 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |
| HIMB1413  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 324 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |
| HIMB1420  | o__Rickettsiales (UID3809) | 83        | 323       | 210           | 0  | 323 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |
| HIMB1427  | o__Rickettsiales (UID3809) | 83        | 323       | 210           | 0  | 311 | 12 | 0 | 0 | 0 | 100.00       | 4.29          | 91.67                |
| HIMB1430  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 5  | 318 | 1  | 0 | 0 | 0 | 97.63        | 0.47          | 0.00                 |
| HIMB1436  | o__Rickettsiales (UID3809) | 83        | 323       | 210           | 0  | 322 | 1  | 0 | 0 | 0 | 100.00       | 0.48          | 0.00                 |
| HIMB1437  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 1  | 323 | 0  | 0 | 0 | 0 | 99.53        | 0.00          | 0.00                 |
| HIMB1444  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 1  | 321 | 2  | 0 | 0 | 0 | 99.76        | 0.63          | 0.00                 |
| HIMB1456  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 3  | 321 | 0  | 0 | 0 | 0 | 99.05        | 0.00          | 0.00                 |
| HIMB1483  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 1  | 322 | 1  | 0 | 0 | 0 | 99.53        | 0.47          | 0.00                 |
| HIMB1485  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 323 | 0  | 0 | 1 | 0 | 100.00       | 1.42          | 0.00                 |
| HIMB1488  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 321 | 2  | 1 | 0 | 0 | 100.00       | 1.90          | 0.00                 |
| HIMB1490  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 8  | 315 | 1  | 0 | 0 | 0 | 96.21        | 0.24          | 0.00                 |
| HIMB1491  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 323 | 1  | 0 | 0 | 0 | 100.00       | 0.24          | 0.00                 |
| HIMB1493  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 324 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |
| HIMB1494  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 322 | 2  | 0 | 0 | 0 | 100.00       | 0.95          | 0.00                 |
| HIMB1495  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 322 | 2  | 0 | 0 | 0 | 100.00       | 0.95          | 0.00                 |
| HIMB1505  | o__Rickettsiales (UID3809) | 83        | 323       | 210           | 0  | 320 | 2  | 1 | 0 | 0 | 100.00       | 1.43          | 20.00                |
| HIMB1506  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 1  | 320 | 3  | 0 | 0 | 0 | 99.53        | 1.42          | 0.00                 |
| HIMB1507  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 1  | 322 | 1  | 0 | 0 | 0 | 99.53        | 0.16          | 0.00                 |
| HIMB1509  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 321 | 3  | 0 | 0 | 0 | 100.00       | 1.42          | 0.00                 |
| HIMB1513  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 323 | 1  | 0 | 0 | 0 | 100.00       | 0.47          | 0.00                 |
| HIMB1517  | k__Bacteria (UID2495)      | 2993      | 140       | 84            | 1  | 137 | 2  | 0 | 0 | 0 | 98.81        | 2.38          | 100.00               |
| HIMB1518  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 324 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |
| HIMB1520  | o__Rickettsiales (UID3809) | 83        | 323       | 210           | 0  | 323 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |
| HIMB1521  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 323 | 1  | 0 | 0 | 0 | 100.00       | 0.47          | 0.00                 |
| HIMB1526  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 323 | 1  | 0 | 0 | 0 | 100.00       | 0.47          | 0.00                 |
| HIMB1527  | o__Rickettsiales (UID3809) | 83        | 323       | 210           | 2  | 318 | 3  | 0 | 0 | 0 | 99.29        | 0.95          | 0.00                 |
| HIMB1542  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 324 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |
| HIMB1549  | o__Rickettsiales (UID3809) | 83        | 323       | 210           | 1  | 319 | 3  | 0 | 0 | 0 | 99.52        | 0.71          | 0.00                 |
| HIMB1552  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 323 | 1  | 0 | 0 | 0 | 100.00       | 0.47          | 0.00                 |
| HIMB1556  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 322 | 2  | 0 | 0 | 0 | 100.00       | 0.63          | 0.00                 |
| HIMB1559  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 1  | 321 | 2  | 0 | 0 | 0 | 99.53        | 0.95          | 0.00                 |
| HIMB1564  | o__Rickettsiales (UID3809) | 83        | 323       | 210           | 0  | 323 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |
| HIMB1565  | k__Bacteria (UID2495)      | 2993      | 140       | 84            | 2  | 138 | 0  | 0 | 0 | 0 | 97.62        | 0.00          | 0.00                 |
| HIMB1573  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 324 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |
| HIMB1577  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 1  | 322 | 1  | 0 | 0 | 0 | 99.53        | 0.47          | 0.00                 |
| HIMB1587  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 1  | 321 | 2  | 0 | 0 | 0 | 99.53        | 0.95          | 0.00                 |
| HIMB1593  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 324 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |
| HIMB1597  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 1  | 321 | 2  | 0 | 0 | 0 | 99.53        | 0.95          | 0.00                 |
| HIMB1611  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 323 | 1  | 0 | 0 | 0 | 100.00       | 0.47          | 0.00                 |
| HIMB1623  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 1  | 323 | 0  | 0 | 0 | 0 | 99.53        | 0.00          | 0.00                 |
| HIMB1631  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 1  | 322 | 1  | 0 | 0 | 0 | 99.53        | 0.47          | 0.00                 |
| HIMB1636  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 323 | 1  | 0 | 0 | 0 | 100.00       | 0.47          | 0.00                 |
| HIMB1641  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 323 | 1  | 0 | 0 | 0 | 100.00       | 0.47          | 0.00                 |
| HIMB1662  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 1  | 323 | 0  | 0 | 0 | 0 | 99.53        | 0.00          | 0.00                 |
| HIMB1685  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 324 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |
| HIMB1695  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 322 | 2  | 0 | 0 | 0 | 100.00       | 0.95          | 0.00                 |
| HIMB1701  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 1  | 322 | 1  | 0 | 0 | 0 | 99.53        | 0.47          | 0.00                 |
| HIMB1702  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 323 | 1  | 0 | 0 | 0 | 100.00       | 0.47          | 0.00                 |
| HIMB1709  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 324 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |
| HIMB1710  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 324 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |
| HIMB1715  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 5  | 319 | 0  | 0 | 0 | 0 | 98.10        | 0.00          | 0.00                 |
| HIMB1723  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 324 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |
| HIMB1746  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 323 | 1  | 0 | 0 | 0 | 100.00       | 0.47          | 0.00                 |
| HIMB1748  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 319 | 5  | 0 | 0 | 0 | 100.00       | 1.90          | 100.00               |
| HIMB1758  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 323 | 1  | 0 | 0 | 0 | 100.00       | 0.47          | 0.00                 |
| HIMB1765  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 324 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |
| HIMB1770  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 324 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |
| HIMB1782  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 322 | 2  | 0 | 0 | 0 | 100.00       | 0.95          | 0.00                 |
| HIMB1863  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 315 | 9  | 0 | 0 | 0 | 100.00       | 3.32          | 100.00               |
| HIMB2104  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 324 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |
| HIMB2187  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 324 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |
| HIMB2200  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 323 | 1  | 0 | 0 | 0 | 100.00       | 0.47          | 0.00                 |
| HIMB2201  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 1  | 323 | 0  | 0 | 0 | 0 | 99.76        | 0.00          | 0.00                 |
| HIMB2204  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 1  | 323 | 0  | 0 | 0 | 0 | 99.53        | 0.00          | 0.00                 |
| HIMB2211  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 324 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |
| HIMB2215  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 8  | 316 | 0  | 0 | 0 | 0 | 98.10        | 0.00          | 0.00                 |
| HIMB2226  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 324 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |
| HIMB2250  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 1  | 323 | 0  | 0 | 0 | 0 | 99.53        | 0.00          | 0.00                 |
| HIMB2305  | k__Bacteria (UID2495)      | 2993      | 139       | 83            | 0  | 139 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |
| HIMB4     | o__Rickettsiales (UID3809) | 83        | 323       | 210           | 0  | 317 | 5  | 1 | 0 | 0 | 100.00       | 2.30          | 12.50                |
| HIMB5     | o__Rickettsiales (UID3809) | 83        | 323       | 210           | 0  | 322 | 1  | 0 | 0 | 0 | 100.00       | 0.48          | 0.00                 |
| HIMB83    | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 1  | 323 | 0  | 0 | 0 | 0 | 99.53        | 0.00          | 0.00                 |
| HTCC1002  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 323 | 1  | 0 | 0 | 0 | 100.00       | 0.47          | 0.00                 |
| HTCC1013  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 323 | 1  | 0 | 0 | 0 | 100.00       | 0.47          | 0.00                 |
| HTCC1016  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 323 | 1  | 0 | 0 | 0 | 100.00       | 0.47          | 0.00                 |
| HTCC1040  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 323 | 1  | 0 | 0 | 0 | 100.00       | 0.47          | 0.00                 |
| HTCC1062  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 323 | 1  | 0 | 0 | 0 | 100.00       | 0.47          | 0.00                 |
| HTCC7211  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 323 | 1  | 0 | 0 | 0 | 100.00       | 0.47          | 0.00                 |
| HTCC7214  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 2  | 321 | 1  | 0 | 0 | 0 | 99.05        | 0.47          | 0.00                 |
| HTCC7217  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 322 | 2  | 0 | 0 | 0 | 100.00       | 0.95          | 0.00                 |
| HTCC8051  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 1  | 323 | 0  | 0 | 0 | 0 | 99.76        | 0.00          | 0.00                 |
| HTCC9022  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 324 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |
| HTCC9565  | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 2  | 322 | 0  | 0 | 0 | 0 | 99.05        | 0.00          | 0.00                 |
| IMCC9063  | k__Bacteria (UID2495)      | 2993      | 140       | 84            | 1  | 139 | 0  | 0 | 0 | 0 | 98.81        | 0.00          | 0.00                 |
| LSUCC0530 | k__Bacteria (UID2495)      | 2993      | 139       | 83            | 0  | 139 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |
| NP1       | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 2  | 320 | 2  | 0 | 0 | 0 | 99.05        | 0.63          | 0.00                 |
| RS39      | o__Rickettsiales (UID3809) | 83        | 324       | 211           | 0  | 324 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |
| RS40      | k__Bacteria (UID2495)      | 2993      | 139       | 83            | 0  | 139 | 0  | 0 | 0 | 0 | 100.00       | 0.00          | 0.00                 |

</details>

---

Looking at the output, this is what we see: 
- Besides HIMB123 (completeness: 86.73), all isolate genomes have a completeness of ≥96.00. However, even 86.73 is sufficient for our purposes.
- The highest contamination value is HIMB1427 (4.29%), followed by HIMB1863 (3.32%), HIMB1517 (2,38%), HIMB4 (2.30%), and HIMB1748, HIMB1488, and HIMB123 (all 1.90%), HIMB1505 (1.43%) HIMB1509 and HIMB1506, HIMB1485 (all 1.42 %). All others are below 1% contamination.

We kept all isolate genomes since they are all above 80% completeness and below 5% contamination.

### Dereplicating the isolate genomes

In the following section, we describe how we dereplicated the 99 reference genomes to retain a single representative genome from each group of highly similar genomes. Dereplication is useful for reducing redundancy, optimizing computational efficiency, and ensuring that downstream analyses focus on the unique diversity of the genome set. 

To accomplish this, we used `dRep`, a Python-based tool designed for rapid pair-wise comparison and clustering of genomes. By utilizing `dRep`, one can identify and retain one representative genome per cluster based on similarity thresholds, ensuring that only distinct genomes are preserved for further analysis.

<div class="extra-info" markdown="1">

<span class="extra-info-header">dRep citation and version</span>

Olm, M., Brown, C., Brooks, B. et al. dRep: a tool for fast and accurate genomic comparisons that enables improved genome recovery from metagenomes through de-replication. ISME J 11, 2864–2868 (2017). https://doi.org/10.1038/ismej.2017.126

We used dRep v3.5.0

</div>

To use `dRep`, we wrote a bash script, telling it to dereplicate at 95% ANI (preclustering at 90% ANI) and then submit it with `clusterize`. 

{:.notice}
ANI stands for Average Nucleotide Identity. `dRep` uses this measure to cluster genomes before selecting a representative for each cluster. https://drep.readthedocs.io/en/latest/choosing_parameters.html

``` bash
nano drep_clusterize.sh
```
content
``` bash
#!/bin/bash 

# Load necessary modules 
# Make sure to load the dRep module if your environment uses module systems
module load dRep

# Set your input directory and output directory
INPUT_DIR=fastaOriginal  # Directory containing the input FASTA files
OUTPUT_DIR=fastaDrep     # Directory where dRep output will be saved

# Step 1: Dereplicate genomes at 95% ANI
# - OUTPUT_DIR: specifies the directory where the output will be saved
# - -g $INPUT_DIR/*.fa -g $INPUT_DIR/*.fna: includes all .fa and .fna files from the input directory
# - -sa 0.95: sets the secondary ANI threshold for dereplication to 95%
# - -pa 0.9: sets the primary ANI threshold for initial clustering to 90%
# - --ignoreGenomeQuality: ignores genome quality. Usually, one could use -comp to set the minimum genome completeness for dereplication and -con to set the maximum genome contamination for dereplication
dRep dereplicate $OUTPUT_DIR -g $INPUT_DIR/*.fa -sa 0.95 -pa 0.9 --ignoreGenomeQuality
```

Submit the job with `clusterize`.
``` bash
# clusterize command explanation:
# -j dRep_workflow: Sets the job name to 'dRep_workflow'
# -o dRep_workflow.log: Specifies the log file for the job output
# -n 16: Requests 16 CPU cores for the job
# "bash drep_clusterize.sh": Runs the 'drep_clusterize.sh' script
clusterize -j dRep_workflow \
           -o dRep_workflow.log \
           -n 1 \
           "bash drep_clusterize.sh"
```

**51 genomes passed.** These are the ones we will continue with.

If you would like to see all `dRep` output files, feel free to do so by clicking [here](files/output_dRep.zip) (zipped). 

---
<details markdown="1"><summary>Click to show/hide primary clustering dendrogram</summary>

We are showing the primary clustering dendrogram here because it gives a more digestible overview, for more accurate clustering information, please consult the secondary clustering dendrogram included in the dRep output folder linked above.

blue and purple stars: representatives after dereplicating. These are the reference genomes we will continue with.

[![dRep]({{images}}/Primary_clustering_dendrogram.pdf)]({{images}}/Primary_clustering_dendrogram.pdf){:.center-img .width-70}


</details>

---

### Concatenating FASTA files of dereplicated reference genomes into one and simplifying deflines

We combined all the FASTA files of the dereplicated reference genomes into a single file because, later on, we performed competitive read recruitment. For this process to work correctly, all reference genomes needed to be in one unified FASTA file to ensure that the reads could be compared against all genomes in a competitive manner.

To use anvi'o (which does not allow special characters in the deflines of FASTA files) and to easily identify each reference genome throughout the analysis, we simplified the deflines using anvi'o's [`anvi-script-reformat-fasta`](https://anvio.org/help/main/programs/anvi-script-reformat-fasta/) program. We used the `--prefix` flag to include the name of the reference genome, and the `--simplify-names` flag to remove any special characters and standardize the names. Additionally, the `--report-file` flag was used to generate a file that maps the original deflines in the FASTA files to the reformatted ones. Finally, we concatenated all FASTA files with these reformatted deflines into a single file named `all_fasta.fa`.

We are running the following commands in the directory in which the FASTA files are. 

```bash
# assuming each .fa file is named according to genome name

for g in *.fa; do
  name=${g%.fa}
  anvi-script-reformat-fasta --prefix $name --simplify-names --report-file ${name}-reformat-report.txt -o ${name}_reformatted.fa $g
done

# afterwards, concatenate
cat *_reformatted.fa > all_fasta.fa
```

To check if it worked, we `grep`ed for the number of carrots (>) across input files VS in the all_fasta.fa file. The values should match exactly.

```bash
grep -c ">" *reformatted.fa | awk -F: '{s+=$2} END {print s}'
grep -c ">" all_fasta.fa
```




## Metagenomes

This section explains how to download and quality filter short metagenomic reads from the Hawaii Ocean Time-Series ([Biller et al., 2018](https://doi.org/10.1038/sdata.2018.176), [Mende et al., 2017](https://doi.org/10.1038/s41564-017-0008-3)), the Bermuda Atlantic Time-series Study ([Biller et al. 2018](https://doi.org/10.1038/sdata.2018.176)) OceaTARA Oceans project ([Sunagawa et al., 2015](https://doi.org/10.1126/science.1261359)), the Malaspina Expedition ([Sánchez et al., 2024](https://doi.org/10.1038/s41597-024-02974-1)), the Bio-GO-SHIP project ([Larkin et al., 2021](https://doi.org/10.1038/s41597-021-00889-9)), the bioGEOTRACES project ([Biller et al. 2018](https://doi.org/10.1038/sdata.2018.176)), and the Ocean Sampling Day project ([Kopf et al., 2015](https://doi.org/10.1186/s13742-015-0066-5)), as used in this analysis.

{:.notice}
All metagenomes we analyzed are publicly available through the European Nucleotide Archive ([ENA](https://www.ebi.ac.uk/ena/browser/home)) and [NCBI](https://www.ncbi.nlm.nih.gov). The accession numbers to each project are given in the Summary table at the top of this workflow.

{:.notice}
If you are reproducing this workflow, please note that we are using metagenomes from a total of **1826** samples. Make sure you have the storage and computational resources to handle them.

### Downloading the metagenomes

To download the metagenomes, we used anvi'o's [`sra_download`](https://anvio.org/help/main/workflows/sra-download/). For this, we needed a `SRA_accession_list.txt` for each project and a `download_config.json` config file.

#### Prep the download

The `SRA_accession_list.txt` artifact should look like this: no headers, only a list of run accession numbers.
```bash
$ cat SRA_accession_list.txt
ERR6450080
ERR6450081
SRR5965623
```

We created multiple `SRA_accession_list.txt`s, one per project. Or more accurately, they were created as part of our work in the [public-marine-omics-metadata](https://github.com/merenlab/public-marine-omics-metadata/tree/main) GitHub repository. That way, they only include the accession numbers to the samples that are in accordance with our metadata requirements. 

To get the files, click [here](files/SRA_accession-txt.zip) for a zipped folder including them.

For the different projects, they are called: 
```bash
SRA_accession_PRJEB1787_TARA.txt
SRA_accession_PRJEB8682_OSD.txt
SRA_accession_PRJEB52452_MAL.txt
SRA_accession_PRJNA385854_BGT.txt
SRA_accession_PRJNA385855_BATS.txt
SRA_accession_PRJNA385855_PRJNA352737_HOT_combined.txt
SRA_accession_PRJNA656268_BGS.txt
```

FYI, we created a directory for each of the observatories, with a `00_WORKFLOW_FILES` directory, to which these files were added.

In addition to the `SRA_accession_list.txt`s, we needed a `download_config.json` for the [`sra_download`](https://anvio.org/help/main/workflows/sra-download/) function. To get the config file, we ran:
```bash
anvi-run-workflow -w sra_download --get-default-config download_config.json
```

In each of the `donwload_config.json` files (again, one per project), we exchanged the `"SRA_accession_list": SRA_accession_list.txt",` with `"SRA_accession_list": "00_WORKFLOW_FILES/SRA_accession_[PROJECT_ACCESSION_NUMBER]_[PROJECT_ACRONYM].txt",` (of course adding the respective project accession number and project acronym).
```bash
{
    "SRA_accession_list": "00_WORKFLOW_FILES/SRA_accession_[PROJECT_ACCESSION_NUMER]_[PROJECT_ACRONYM].txt",
    "prefetch": {
        "--max-size": "40g",
        "threads": 2
    },
    "fasterq_dump": {
        "threads": 6
    },
    "pigz": {
        "threads": 8
    },
    "output_dirs": {
        "SRA_prefetch": "01_NCBI_SRA",
        "FASTAS": "02_FASTA",
        "LOGS_DIR": "00_LOGS"
    },
    "max_threads": "",
    "config_version": "3",
    "workflow_name": "sra_download"
}
```

#### Do the download

We then were ready to use this `download_config.json` in combination with anvi'o's [`anvi-run-workflow`](https://anvio.org/help/main/programs/anvi-run-workflow/) and the workflow [`sra_download`](https://anvio.org/help/main/workflows/sra-download/).

We did a dry run ...
```bash
anvi-run-workflow -w sra_download -c 00_WORKFLOW_FILES/download_config.json -A -n -q
```
... and then sent the real deal off to get us the metagenomes!
```bash
clusterize -j sra_download_workflow \
                 -o sra_download.log \
                 -n 1 \
                 --mail-user raissa.meyer@awi.de "anvi-run-workflow -w sra_download \
    -c 00_WORKFLOW_FILES/download_config.json  -A \
   --cluster 'clusterize -o {log} -n {threads}' \
   --resources nodes=100 \
   --jobs 100 \
   --rerun-incomplete \
   --keep-going"
```

Following the download, we searched for the mention of any errors in the `.log` files to see if everything went as expected.
```bash
grep -i ERROR sra_download_*.log
```

We additionally used the following script to check if there is a file for each accession number in the `SRA_accession_list.txt`s.
```bash
nano check_fastq_files.sh
```
script
```bash
#!/bin/bash

# Path to the file containing accession numbers
accession_file="00_WORKFLOW_FILES/SRA_accession_PRJNA385855_BATS.txt"

# Directory where the FASTQ files are stored
fastq_dir="02_FASTA"

# Flag to track if all files are present
all_files_exist=true

# Read each line (accession number) from the file
while IFS= read -r accession; do
  # Construct the expected file names for both read pairs
  file1="${fastq_dir}/${accession}_1.fastq.gz"
  file2="${fastq_dir}/${accession}_2.fastq.gz"
  
  # Check if the first file exists
  if [[ ! -f "$file1" ]]; then
    echo "Missing $file1"
    all_files_exist=false
  else
    echo "Found $file1"
  fi
  
  # Check if the second file exists
  if [[ ! -f "$file2" ]]; then
    echo "Missing $file2"
    all_files_exist=false
  else
    echo "Found $file2"
  fi
done < "$accession_file"

# Final feedback
if $all_files_exist; then
  echo "All files are present."
else
  echo "Some files are missing."
fi

```
give permissions for executing the script
```bash
chmod +x check_fastq_files.sh
```
run 
```bash
./check_fastq_files.sh
```
Everything was well.


### QC

For QC, we used part of the anvi'o [`metagenomics` workflow](https://anvio.org/help/main/workflows/metagenomics/) to remove noise from raw reads prior to mapping. Namely, the `illumina-utils` program ([Eren et al., 2013](https://doi.org/10.1371/journal.pone.0066643)).

For that, we needed a [`samples-txt`](https://anvio.org/help/main/artifacts/samples-txt/) file following this structure (the group column is added as a bonus if one needs it, but not necessary):
```bash
column -t samples.txt
sample     group  r1                                           r2
sample_01  G01    three_samples_example/sample-01-R1.fastq.gz  three_samples_example/sample-01-R2.fastq.gz
sample_02  G02    three_samples_example/sample-02-R1.fastq.gz  three_samples_example/sample-02-R2.fastq.gz
sample_03  G02    three_samples_example/sample-03-R1.fastq.gz  three_samples_example/sample-03-R2.fastq.gz
```
Generally, the file needs to include the sample name and the locations of raw paired-end reads for each sample.


#### Prep QC

In the following, we describe how we created one of those `samples.txt` files. We did it in two steps.

**1. Create a more extensive mapping of biosample to run to project to custom_sample_name**

We created a detailed mapping (which also has the Project accession, real Sample accession, and Run accession, as well as a custom_sample_name [created based on this schema: prefix: [PROJECT_ACRONYM], then the value in `sample_accession`, then the value in `depth`, then the value in `collection_date`, all separated by underscores]) for each project. These mapping files are accessible [here](files/BioSample_to_SRA_accession-csv.zip) (zipped).

{:.notice}
Note, that we are using files created as part of the [public-marine-omics-metadata](https://github.com/merenlab/public-marine-omics-metadata/tree/main) GitHub repository: the `_patch2.csv` files. 

From the `_patch2.csv` files, which contain metadata for each project, we used the following script to get the extensive mapping files.

```bash
nano sampleMapping.py
```
script
```python
import pandas as pd
import re
import glob
import os

# Get all files matching the pattern *_patch2.csv in the directory ../../data/
files = glob.glob("../data/*_patch2.csv")

def extract_numeric_depth(depth_value):
    match = re.match(r'(\d+)', str(depth_value))
    return int(match.group(1)) if match else None

def process_file(file):
    # Load the metadata file
    df = pd.read_csv(file)

    # Create a new DataFrame with unique biosamples
    biosamples = df['biosample'].unique().tolist()
    bios = pd.DataFrame(index=biosamples, columns=["run", "bioproject", "custom_sample_name"])

    # Extract the base filename to use in custom sample names
    base_filename = os.path.basename(file).replace('_patch2.csv', '')

    # Populate the new DataFrame
    for b in biosamples:
        # Get the corresponding run accessions and join them into a comma-separated string
        bios.loc[b, "run"] = ",".join(df[df['biosample'] == b]['run'].tolist())

        # Get the corresponding bioproject (assuming they are the same for all run accessions of a biosample)
        bios.loc[b, "bioproject"] = df[df['biosample'] == b]['bioproject'].iloc[0]

        # Get the depth and collection date for the biosample (assuming they are the same for all run accessions of a biosample)
        depth = df[df['biosample'] == b]['depth'].iloc[0]
        depth_int = extract_numeric_depth(depth)  # Extract numeric part of the depth
        collection_date = df[df['biosample'] == b]['collection_date'].iloc[0]
        collection_date_formatted = pd.to_datetime(collection_date).strftime('%Y-%m-%d')  # Format the date to Y-M-D

        # Create the custom sample name using the base_filename
        bios.loc[b, "custom_sample_name"] = f"{base_filename}_{b}_{depth_int}_{collection_date_formatted}"

    # Reset the index to have 'biosample' as a column
    bios.reset_index(inplace=True)
    bios.rename(columns={'index': 'biosample'}, inplace=True)

    # Define the output filename, removing "_patch2" from the base filename
    output_filename = f"../data/BioSample_to_SRA_accessions_{base_filename}.csv"

    # Save the new DataFrame to a CSV file
    bios.to_csv(output_filename, sep="\t", index=False)

# Process each file
for file in files:
    process_file(file)
```

You will get the following files
```
BioSample_to_SRA_accessions_BATS.csv
BioSample_to_SRA_accessions_BGS.csv
BioSample_to_SRA_accessions_BGT.csv
BioSample_to_SRA_accessions_HOT1.csv
BioSample_to_SRA_accessions_HOT3.csv
BioSample_to_SRA_accessions_MAL.csv
BioSample_to_SRA_accessions_OSD.csv
BioSample_to_SRA_accessions_TARA.csv
```

**2. Create the sample-txt artifact from that and the files in our 02_FASTA directories**

Following Step 1., we created the `samples-txt` artifacts. To get those files, you can click [here](files/samples-txt.zip) (zipped) and select the `samples_raw.txt` files.

To create those `samples-txt` files, we had to account for there being two `BioSample_to_SRA_accession_[PROJECT_ACRONYM].csv` files in the HOT directory. Further, anvi'o does not appreciate the use of '-' in sample names (e.g., for the date portion of the sample name), so those had to be substituted with '_'.

We created and ran this script in the directory above our project directories.
```bash
nano makeSamples-txt.py
```
script
```python
import pandas as pd
import os
import glob

# Directories to process
directories = ["HOT", "BATS", "BIOGEOTRACES", "BIOGOSHIP", "MALASPINA", "OSDay", "TARA"]

# Function to create samples_raw.txt from the CSV files
def create_samples_raw_file(directory):
    # Find the relevant CSV files in the directory
    csv_file_pattern = os.path.join(directory, "BioSample_to_SRA_accessions_*.csv")
    csv_files = glob.glob(csv_file_pattern)

    if not csv_files:
        print(f"No CSV file found in {directory}")
        return
    
    for csv_file in csv_files:
        df = pd.read_csv(csv_file, sep="\t")

        # Check if the required columns are in the DataFrame
        expected_columns = ['custom_sample_name', 'run']
        for col in expected_columns:
            if col not in df.columns:
                print(f"Column '{col}' is missing in {csv_file}. Available columns: {df.columns.tolist()}")
                return
        
        # Directory containing the FASTQ files
        fastq_dir = os.path.join(directory, "02_FASTA")

        # List to hold the rows for the samples_raw.txt file
        samples_raw_data = []

        # Check for each run accession
        for index, row in df.iterrows():
            custom_sample_name = row['custom_sample_name']
            run_accessions = row['run'].split(',')
            
            for run_accession in run_accessions:
                r1_file = os.path.abspath(os.path.join(fastq_dir, f"{run_accession}_1.fastq.gz"))
                r2_file = os.path.abspath(os.path.join(fastq_dir, f"{run_accession}_2.fastq.gz"))
                
                if os.path.exists(r1_file) and os.path.exists(r2_file):
                    samples_raw_data.append({
                        "sample": custom_sample_name,
                        "r1": r1_file,
                        "r2": r2_file
                    })
                    print(f"Files found for run accession {run_accession}:")
                    print(f"  {r1_file}")
                    print(f"  {r2_file}")
                else:
                    print(f"Files missing for run accession {run_accession}:")
                    if not os.path.exists(r1_file):
                        print(f"  Missing: {r1_file}")
                    if not os.path.exists(r2_file):
                        print(f"  Missing: {r2_file}")

        # Create a DataFrame from the collected data
        samples_raw_df = pd.DataFrame(samples_raw_data)

        # Determine the output file name
        if directory == "HOT":
            # Specific filenames for HOT directory
            base_filename = os.path.basename(csv_file).replace("BioSample_to_SRA_accessions_", "").replace(".csv", "")
            output_file_path = os.path.join(directory, f"samples_raw_{base_filename}.txt")
        else:
            # Generic samples_raw.txt for other directories
            output_file_path = os.path.join(directory, "samples_raw.txt")
        
        # Save the DataFrame to a .txt file
        samples_raw_df.to_csv(output_file_path, sep="\t", index=False)

        print(f"{output_file_path} file has been created successfully.")

# Function to replace hyphens with underscores in sample names in the samples_raw.txt files
def process_samples_file(directory):
    # Determine the pattern to search for the appropriate files to process
    if directory == "HOT":
        sample_raw_files = glob.glob(os.path.join(directory, "samples_raw_*.txt"))
    else:
        sample_raw_files = [os.path.join(directory, "samples_raw.txt")]

    for file_path in sample_raw_files:
        # Check if the file exists
        if not os.path.isfile(file_path):
            print(f"No {os.path.basename(file_path)} file found in {directory}")
            continue

        # Read the contents of the file
        with open(file_path, 'r') as file:
            lines = file.readlines()

        # Replace hyphens with underscores in sample names
        modified_lines = []
        for line in lines:
            modified_line = line.replace("-", "_")
            modified_lines.append(modified_line)

        # Write the modified contents back to the file
        with open(file_path, 'w') as file:
            file.writelines(modified_lines)
        
        print(f"Processed {file_path}")

# Process each directory
for directory in directories:
    create_samples_raw_file(directory)
    process_samples_file(directory)

print("Completed processing all directories.")
```

After this step, we also counted how many samples and runs were downloaded.
```bash
nano countSamplesNruns.py
```
script
```python
import pandas as pd
import os
import glob

# Directories to process
directories = ["HOT", "BATS", "BIOGEOTRACES", "BIOGOSHIP", "MALASPINA", "OSDay", "TARA"]

# Function to count and report the number of unique samples and runs
def count_samples_and_runs(directory):
    # Determine the pattern to search for the appropriate files to process
    if directory == "HOT":
        sample_raw_files = glob.glob(os.path.join(directory, "samples_raw_*.txt"))
    else:
        sample_raw_files = [os.path.join(directory, "samples_raw.txt")]

    for file_path in sample_raw_files:
        # Check if the file exists
        if not os.path.isfile(file_path):
            print(f"No {os.path.basename(file_path)} file found in {directory}")
            continue
        
        # Check if the file is empty
        if os.path.getsize(file_path) == 0:
            print(f"{file_path} is empty. Skipping.")
            continue

        # Load the samples_raw.txt file
        try:
            df = pd.read_csv(file_path, sep="\t")
        except pd.errors.EmptyDataError:
            print(f"Error: {file_path} is empty or corrupted. Skipping.")
            continue

        # Count the number of unique samples and total runs
        unique_samples = df['sample'].nunique()
        total_runs = df.shape[0]

        # Report the counts
        print(f"{file_path} contains {unique_samples} unique samples and {total_runs} runs.")

# Process each directory
for directory in directories:
    count_samples_and_runs(directory)

print("Completed counting samples and runs in all directories.")
```
run
```
python3 countSamplesNruns.py
```


#### Do QC
With that, we were ready to get the `config file` we needed to perform the QC step.

To get the `config file` needed to perfom the quality control step, we ran
```bash
anvi-run-workflow -w metagenomics --get-default-config QC_config.json 
```

In this `config file`, we
- replaced `"samples_txt": "samples.txt"` with `"samples_txt": "samples_raw.txt"` 
- turned off everything besides 
   - `idba_ud` (needs to be on to trick snakemake, but we're not actually running it, because when submitting the job we tell it to stop after gzipping), 
   - `iu_filter_quality_minoche` (and set threads to 4), and 
   - `gzip_fastqs`
   - `anvi_script_reformat_fasta` (needs to be on for anvi'o not to complain but is not used here)

---

<details markdown="1"><summary>Click here to show/hide the content of the contigs file</summary>

```bash
{
    "anvi_gen_contigs_database": {
        "--project-name": "{group}",
        "--description": "",
        "--skip-gene-calling": "",
        "--ignore-internal-stop-codons": "",
        "--skip-mindful-splitting": "",
        "--contigs-fasta": "",
        "--split-length": "",
        "--kmer-size": "",
        "--skip-predict-frame": "",
        "--prodigal-translation-table": "",
        "threads": ""
    },
    "centrifuge": {
        "threads": 2,
        "run": "",
        "db": ""
    },
    "anvi_run_hmms": {
        "run": false,
        "threads": 5,
        "--also-scan-trnas": true,
        "--installed-hmm-profile": "",
        "--hmm-profile-dir": "",
        "--add-to-functions-table": ""
    },
    "anvi_run_kegg_kofams": {
        "run": false,
        "threads": 4,
        "--kegg-data-dir": "",
        "--hmmer-program": "",
        "--keep-all-hits": "",
        "--log-bitscores": "",
        "--just-do-it": ""
    },
    "anvi_run_ncbi_cogs": {
        "run": false,
        "threads": 5,
        "--cog-data-dir": "",
        "--temporary-dir-path": "",
        "--search-with": ""
    },
    "anvi_run_scg_taxonomy": {
        "run": false,
        "threads": 6,
        "--scgs-taxonomy-data-dir": ""
    },
    "anvi_run_trna_scan": {
        "run": false,
        "threads": 6,
        "--trna-cutoff-score": ""
    },
    "anvi_script_reformat_fasta": {
        "run": true,
        "--prefix": "{group}",
        "--simplify-names": true,
        "--keep-ids": "",
        "--exclude-ids": "",
        "--min-len": "",
        "--seq-type": "",
        "threads": ""
    },
    "emapper": {
        "--database": "bact",
        "--usemem": true,
        "--override": true,
        "path_to_emapper_dir": "",
        "threads": ""
    },
    "anvi_script_run_eggnog_mapper": {
        "--use-version": "0.12.6",
        "run": "",
        "--cog-data-dir": "",
        "--drop-previous-annotations": "",
        "threads": ""
    },
    "samples_txt": "samples_raw.txt",
    "metaspades": {
        "additional_params": "--only-assembler",
        "threads": 7,
        "run": "",
        "use_scaffolds": ""
    },
    "megahit": {
        "--min-contig-len": 1000,
        "--memory": 0.4,
        "threads": 7,
        "run": "",
        "--min-count": "",
        "--k-min": "",
        "--k-max": "",
        "--k-step": "",
        "--k-list": "",
        "--no-mercy": "",
        "--no-bubble": "",
        "--merge-level": "",
        "--prune-level": "",
        "--prune-depth": "",
        "--low-local-ratio": "",
        "--max-tip-len": "",
        "--no-local": "",
        "--kmin-1pass": "",
        "--presets": "",
        "--mem-flag": "",
        "--use-gpu": "",
        "--gpu-mem": "",
        "--keep-tmp-files": "",
        "--tmp-dir": "",
        "--continue": "",
        "--verbose": ""
    },
    "idba_ud": {
        "--min_contig": 1000,
        "threads": 7,
        "run": true,
        "--mink": "",
        "--maxk": "",
        "--step": "",
        "--inner_mink": "",
        "--inner_step": "",
        "--prefix": "",
        "--min_count": "",
        "--min_support": "",
        "--seed_kmer": "",
        "--similar": "",
        "--max_mismatch": "",
        "--min_pairs": "",
        "--no_bubble": "",
        "--no_local": "",
        "--no_coverage": "",
        "--no_correct": "",
        "--pre_correction": "",
        "use_scaffolds": ""
    },
    "iu_filter_quality_minoche": {
        "run": true,
        "--ignore-deflines": true,
        "--visualize-quality-curves": "",
        "--limit-num-pairs": "",
        "--print-qual-scores": "",
        "--store-read-fate": "",
        "threads": 4
    },
    "gzip_fastqs": {
        "run": true,
        "threads": 2
    },
    "bowtie": {
        "additional_params": "--no-unal",
        "threads": 3
    },
    "samtools_view": {
        "additional_params": "-F 4",
        "threads": ""
    },
    "anvi_profile": {
        "threads": 3,
        "--sample-name": "{sample}",
        "--overwrite-output-destinations": true,
        "--report-variability-full": "",
        "--skip-SNV-profiling": "",
        "--profile-SCVs": "",
        "--description": "",
        "--skip-hierarchical-clustering": "",
        "--distance": "",
        "--linkage": "",
        "--min-contig-length": "",
        "--min-mean-coverage": "",
        "--min-coverage-for-variability": "",
        "--cluster-contigs": "",
        "--contigs-of-interest": "",
        "--queue-size": "",
        "--write-buffer-size-per-thread": "",
        "--fetch-filter": "",
        "--min-percent-identity": "",
        "--max-contig-length": ""
    },
    "anvi_merge": {
        "--sample-name": "{group}",
        "--overwrite-output-destinations": true,
        "--description": "",
        "--skip-hierarchical-clustering": "",
        "--enforce-hierarchical-clustering": "",
        "--distance": "",
        "--linkage": "",
        "threads": ""
    },
    "import_percent_of_reads_mapped": {
        "run": true,
        "threads": ""
    },
    "krakenuniq": {
        "threads": 3,
        "--gzip-compressed": true,
        "additional_params": "",
        "run": "",
        "--db": ""
    },
    "remove_short_reads_based_on_references": {
        "delimiter-for-iu-remove-ids-from-fastq": " ",
        "dont_remove_just_map": "",
        "references_for_removal_txt": "",
        "threads": ""
    },
    "anvi_cluster_contigs": {
        "--collection-name": "{driver}",
        "run": "",
        "--driver": "",
        "--just-do-it": "",
        "--additional-params-concoct": "",
        "--additional-params-metabat2": "",
        "--additional-params-maxbin2": "",
        "--additional-params-dastool": "",
        "--additional-params-binsanity": "",
        "threads": ""
    },
    "gen_external_genome_file": {
        "threads": ""
    },
    "export_gene_calls_for_centrifuge": {
        "threads": ""
    },
    "anvi_import_taxonomy_for_genes": {
        "threads": ""
    },
    "annotate_contigs_database": {
        "threads": ""
    },
    "anvi_get_sequences_for_gene_calls": {
        "threads": ""
    },
    "gunzip_fasta": {
        "threads": ""
    },
    "reformat_external_gene_calls_table": {
        "threads": ""
    },
    "reformat_external_functions": {
        "threads": ""
    },
    "import_external_functions": {
        "threads": ""
    },
    "anvi_run_pfams": {
        "run": "",
        "--pfam-data-dir": "",
        "threads": ""
    },
    "iu_gen_configs": {
        "--r1-prefix": "",
        "--r2-prefix": "",
        "threads": ""
    },
    "gen_qc_report": {
        "threads": ""
    },
    "merge_fastqs_for_co_assembly": {
        "threads": ""
    },
    "merge_fastas_for_co_assembly": {
        "threads": ""
    },
    "bowtie_build": {
        "additional_params": "",
        "threads": ""
    },
    "anvi_init_bam": {
        "threads": ""
    },
    "krakenuniq_mpa_report": {
        "threads": ""
    },
    "import_krakenuniq_taxonomy": {
        "--min-abundance": "",
        "threads": ""
    },
    "anvi_summarize": {
        "additional_params": "",
        "run": "",
        "threads": ""
    },
    "anvi_split": {
        "additional_params": "",
        "run": "",
        "threads": ""
    },
    "references_mode": "",
    "all_against_all": "",
    "kraken_txt": "",
    "collections_txt": "",
    "output_dirs": {
        "FASTA_DIR": "02_FASTA",
        "CONTIGS_DIR": "03_CONTIGS",
        "QC_DIR": "01_QC",
        "MAPPING_DIR": "04_MAPPING",
        "PROFILE_DIR": "05_ANVIO_PROFILE",
        "MERGE_DIR": "06_MERGED",
        "TAXONOMY_DIR": "07_TAXONOMY",
        "SUMMARY_DIR": "08_SUMMARY",
        "SPLIT_PROFILES_DIR": "09_SPLIT_PROFILES",
        "LOGS_DIR": "00_LOGS"
    },
    "max_threads": "",
    "config_version": "3",
    "workflow_name": "metagenomics"
}

```

</details>

---

We did a dry run to check the stats before letting it loose.

```
anvi-run-workflow -w metagenomics -c QC_config.json -A --until gzip_fastqs -n -q
```

Then, we submitted the job using `clusterize`.

{:.notice}
Note that we used the flag `--unil gzip-fastqs` so that we would not run the entire `metagenomics` workflow but only the quality control portion.

```bash
clusterize -j QC_workflow \
                 -o QC_workflow.log \
                 -n 1 \
                 "anvi-run-workflow -w metagenomics \
    -c QC_config.json   -A \
   --until gzip_fastqs \
   --cluster 'clusterize -o {log} -n {threads} -j {rule} ' \
   --resources nodes=500 \
   --jobs 500 \
   --rerun-incomplete \
   --keep-going"
```
Once it was done, we `grep`ed for errors in the log file.
```bash
grep -i ERROR QC_workflow_*.log
```
All was well. The outputs of this workflow are available in the [QC_stats.zip](files/QC_stats.zip) archive.

#### Making a samples_qc.txt

Since we split the metagenomics workflow in two (1. QC, 2. read recruitment), we then created `samples_qc.txt`s (again, one per project), which will be concatenated into one in the step after this, and then used as the `samples-txt` artifact in the steps following that. The projects' `samples_qc.txt`s are also included in the [samples-txt folder](files/samples-txt.zip), in case you would like to have a peek. 

```
nano makeQCsample-txt.py
```
script
```python
import os
import pandas as pd

# Read the samples_raw.txt file and extract the first column
samples_df = pd.read_csv('samples_raw.txt', sep='\t')
samples = samples_df['sample'].tolist()

# Get the absolute path of the 01_QC/ directory
qc_dir = os.path.abspath('01_QC')

# List all files in the 01_QC/ directory
qc_files = os.listdir(qc_dir)

# Initialize a list to store the results
qc_data = []

# Initialize a list to store missing samples
missing_samples = []

# Check if each sample has both R1 and R2 files
for sample in samples:
    r1_file = f"{sample}-QUALITY_PASSED_R1.fastq.gz"
    r2_file = f"{sample}-QUALITY_PASSED_R2.fastq.gz"
    
    r1_path = os.path.join(qc_dir, r1_file)
    r2_path = os.path.join(qc_dir, r2_file)
    
    if r1_file in qc_files and r2_file in qc_files:
        qc_data.append({
            "sample": sample,
            "r1": r1_path,
            "r2": r2_path
        })
    else:
        missing_samples.append({
            "sample": sample,
            "R1_exists": r1_file in qc_files,
            "R2_exists": r2_file in qc_files
        })

# Create a DataFrame from the qc_data list
qc_df = pd.DataFrame(qc_data)

# Save the DataFrame to a samples_qc.txt file
qc_df.to_csv('samples_qc.txt', sep='\t', index=False)

# Notify if any files were missing
if missing_samples:
    print("Some samples are missing files:")
    for sample in missing_samples:
        print(f"Sample: {sample['sample']}")
        print(f"  R1 exists: {sample['R1_exists']}")
        print(f"  R2 exists: {sample['R2_exists']}")
        print()

print("samples_qc.txt has been created.")
```

#### Creating a collective samples_qc.txt for all projects

For the upcoming read recruitment, we no longer wanted to have the projects separated, so we concatenated the individual samples_qc.txt files we created above into one. The collective `samples_qc.txt` is also included in the [samples-txt folder](files/samples-txt.zip).
```bash
nano combine_and_check_samples_qc.py
```
script
```python
import os

# Define the paths to the files
file_paths = [
    "/fs/dss/groups/agecodatasci/RESOURCES/PUBLIC/METAGENOMES/OCEAN/samples_qc.txt",
    "/fs/dss/groups/agecodatasci/RESOURCES/PUBLIC/METAGENOMES/OCEAN/BATS/samples_qc.txt",
    "/fs/dss/groups/agecodatasci/RESOURCES/PUBLIC/METAGENOMES/OCEAN/BIOGEOTRACES/samples_qc.txt",
    "/fs/dss/groups/agecodatasci/RESOURCES/PUBLIC/METAGENOMES/OCEAN/BIOGOSHIP/samples_qc.txt",
    "/fs/dss/groups/agecodatasci/RESOURCES/PUBLIC/METAGENOMES/OCEAN/HOT/samples_qc_HOT1.txt",
    "/fs/dss/groups/agecodatasci/RESOURCES/PUBLIC/METAGENOMES/OCEAN/HOT/samples_qc_HOT3.txt",
    "/fs/dss/groups/agecodatasci/RESOURCES/PUBLIC/METAGENOMES/OCEAN/MALASPINA/samples_qc.txt",
    "/fs/dss/groups/agecodatasci/RESOURCES/PUBLIC/METAGENOMES/OCEAN/OSDay/samples_qc.txt"
]

# Define the output file path
output_file_path = "/fs/dss/groups/agecodatasci/PROJECTS/MicrObsFuture/samples_qc.txt"

# Initialize row counters
total_rows_before = 0
total_rows_after = 0

# Create or overwrite the output file
with open(output_file_path, 'w') as output_file:
    for i, file_path in enumerate(file_paths):
        with open(file_path, 'r') as input_file:
            lines = input_file.readlines()
            # Count the number of rows in the current file, excluding the header
            total_rows_before += len(lines) - 1 if i != 0 else len(lines)
            # Skip header line if it's not the first file
            if i != 0:
                lines = lines[1:]
            output_file.writelines(lines)

# Count the number of rows in the combined file
with open(output_file_path, 'r') as combined_file:
    total_rows_after = len(combined_file.readlines())

print(f"Total rows before combining: {total_rows_before}")
print(f"Total rows after combining: {total_rows_after}")

# Check if the counts match
if total_rows_before == total_rows_after:
    print("The number of rows matches.")
else:
    print("The number of rows does not match.")

# Check if the new file was created
if os.path.exists(output_file_path):
    print(f"The combined file was successfully created at {output_file_path}")
else:
    print("Failed to create the combined file.")
```
run
```
python3 combine_and_check_samples_qc.py
```

## Mapping metagenomic reads to SAR11 genomes (competitive read recruitment)

This section explains the steps we took to characterize the occurrence of each SAR11 isolate genome in metagenomes.

For that, we again made use of the anvi'o `metagenomics` snakemake workflow.

This time, that included:
- competitive read recruitment with `anvi_run_hmms`, `anvi_run_kegg_kofams`, `anvi_run_ncbi_cogs`, `anvi_run_scg_taxonomy`
- convert SAM -> BAM
- profile BAM file to get anvi'o profiles
- merge into a single anvi'o profile

Since we wanted to map our reads to the reference genomes, we needed not only had to specify the `samples_qc.txt` file, but also the `fasta.txt` file. Thus, we created one with the paths to our `all_fasta.fa` file containing all SAR11 reference genomes. 
```bash
name    path
reference_genomes       ../refGenomesSAR11/fastaDrep/dereplicated_genomes/all_fasta.fa
```

We also needed a new `config file` for the metagenomics workflow. To get that, we ran the following:
```bash
anvi-run-workflow -w metagenomics \
                  --get-default-config default-metagenomics-config.json
```

Importantly, we  
- set `"fasta_txt": "fasta.txt"` (this contains the path to the fasta file containing all reference genomes)
- set `"samples_txt": "samples_qc.txt"` (this contains the paths to all QCed metagenomes)
- put the `reference-mode` to `true` (this way, and by having all reference genomes in a single fasta file, we ensure that the read recruitment is competitive)
- turned on 
   - [`anvi_run_hmm`](https://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-run-hmms) - searches for HMMs against a [contigs-db](https://anvio.org/help/main/artifacts/contigs-db) and stores that information into the contigs-db’s [hmm-hits](https://anvio.org/help/main/artifacts/hmm-hits). 
   - [`anvi_run_kegg_kofams`](https://anvio.org/help/main/programs/anvi-run-kegg-kofams/) - annotates a [contigs-db](https://anvio.org/help/main/artifacts/contigs-db) with HMM hits from KOfam, a database of KEGG Orthologs (KOs). 
   -  [`anvi_run_ncbi_cogs`](https://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-run-ncbi-cogs) - annotate genes in your [contigs-db](https://anvio.org/help/main/artifacts/contigs-db) with [functions](https://anvio.org/help/main/artifacts/functions) from the NCBI’s [Clusters of Orthologus Groups](https://www.ncbi.nlm.nih.gov/COG/)
   - [`anvi_run_scg_taxonomy`](https://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-import-taxonomy)  - finds your single-copy core genes and assigns them taxonomy by searching them against GTDB

Other parts of the `config.json` file did not have to be specifically turned on, but ran automatically:
- `anvi-gen-contigs-database` - profile all contigs for SAR11 isolate genomes, and generate an anvi’o contigs database that stores for each contig the DNA sequence, GC-content, tetranucleotide frequency, and open reading frames
- `bowtie-build` - map short reads from the project metagenomes onto the scaffolds contained in fasta file containing the reference genomes using Bowtie2
- `samtools_view` - store the recruited reads as BAM files using samtools
- `anvi-profile` - process the BAM files and generate anvi’o PROFILE databases that contain the coverage and detection statistics of each SAR11 scaffold in a given metagenome
- `anvi-merge` - generate a merged anvi’o profile database from the individual PROFILE databases

<div class="extra-info" markdown="1">

<span class="extra-info-header">Citations</span>

This workflow relies on much more than just anvi'o or snakemake, but many other peoples' resources and work. Here are the citations to acknowledge that
- `anvi-gen-contigs-database`: Anvi'o will use 'prodigal' by Hyatt et al. 2010 to identify open reading frames in your data. Hyatt, D., Chen, GL., LoCascio, P.F. et al. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics 11, 119 (2010). https://doi.org/10.1186/1471-2105-11-119
- Bowtie2: Langmead, B., Salzberg, S. Fast gapped-read alignment with Bowtie 2. Nat Methods 9, 357–359 (2012). https://doi.org/10.1038/nmeth.1923
- samtools: Heng Li, Bob Handsaker, Alec Wysoker, Tim Fennell, Jue Ruan, Nils Homer, Gabor Marth, Goncalo Abecasis, Richard Durbin, 1000 Genome Project Data Processing Subgroup, The Sequence Alignment/Map format and SAMtools, Bioinformatics, Volume 25, Issue 16, August 2009, Pages 2078–2079, https://doi.org/10.1093/bioinformatics/btp352
- `anvi_run_kegg_kofams`: Anvi'o will annotate your database with the KEGG KOfam database, as described in Aramaki et al. 2020: Takuya Aramaki, Romain Blanc-Mathieu, Hisashi Endo, Koichi Ohkubo, Minoru Kanehisa, Susumu Goto, Hiroyuki Ogata, KofamKOALA: KEGG Ortholog assignment based on profile HMM and adaptive score threshold, Bioinformatics, Volume 36, Issue 7, April 2020, Pages 2251–2252, https://doi.org/10.1093/bioinformatics/btz859
- `anvi_run_ncbi_cogs`: Anvi'o will set up the COG20 version of NCBI COGs from Galperin et al. 2021: Michael Y Galperin, Yuri I Wolf, Kira S Makarova, Roberto Vera Alvarez, David Landsman, Eugene V Koonin, COG database update: focus on microbial diversity, model organisms, and widespread pathogens, Nucleic Acids Research, Volume 49, Issue D1, 8 January 2021, Pages D274–D281, https://doi.org/10.1093/nar/gkaa1018

</div>

---

<details markdown="1"><summary>Click to show/hide default-metagenomics-config.json</summary>


```
{
    "fasta_txt": "fasta.txt",
    "anvi_gen_contigs_database": {
        "--project-name": "{group}",
        "--description": "",
        "--skip-gene-calling": "",
        "--ignore-internal-stop-codons": "",
        "--skip-mindful-splitting": "",
        "--contigs-fasta": "",
        "--split-length": "",
        "--kmer-size": "",
        "--skip-predict-frame": "",
        "--prodigal-translation-table": "",
        "threads": ""
    },
    "centrifuge": {
        "threads": 2,
        "run": "",
        "db": ""
    },
    "anvi_run_hmms": {
        "run": true,
        "threads": 10,
        "--also-scan-trnas": true,
        "--installed-hmm-profile": "",
        "--hmm-profile-dir": "",
        "--add-to-functions-table": ""
    },
    "anvi_run_kegg_kofams": {
        "run": true,
        "threads": 20,
        "--kegg-data-dir": "",
        "--hmmer-program": "",
        "--keep-all-hits": "",
        "--log-bitscores": "",
        "--just-do-it": ""
    },
    "anvi_run_ncbi_cogs": {
        "run": true,
        "threads": 20,
        "--cog-data-dir": "",
        "--temporary-dir-path": "",
        "--search-with": ""
    },
    "anvi_run_scg_taxonomy": {
        "run": true,
        "threads": 6,
        "--scgs-taxonomy-data-dir": ""
    },
    "anvi_run_trna_scan": {
        "run": false,
        "threads": 6,
        "--trna-cutoff-score": ""
    },
    "anvi_script_reformat_fasta": {
        "run": false,
        "--prefix": "{group}",
        "--simplify-names": true,
        "--keep-ids": "",
        "--exclude-ids": "",
        "--min-len": "",
        "--seq-type": "",
        "threads": ""
    },
    "emapper": {
        "--database": "bact",
        "--usemem": true,
        "--override": true,
        "path_to_emapper_dir": "",
        "threads": ""
    },
    "anvi_script_run_eggnog_mapper": {
        "--use-version": "0.12.6",
        "run": "",
        "--cog-data-dir": "",
        "--drop-previous-annotations": "",
        "threads": ""
    },
    "samples_txt": "samples_qc.txt",
    "metaspades": {
        "additional_params": "--only-assembler",
        "threads": 7,
        "run": "",
        "use_scaffolds": ""
    },
    "megahit": {
        "--min-contig-len": 1000,
        "--memory": 0.4,
        "threads": 7,
        "run": "",
        "--min-count": "",
        "--k-min": "",
        "--k-max": "",
        "--k-step": "",
        "--k-list": "",
        "--no-mercy": "",
        "--no-bubble": "",
        "--merge-level": "",
        "--prune-level": "",
        "--prune-depth": "",
        "--low-local-ratio": "",
        "--max-tip-len": "",
        "--no-local": "",
        "--kmin-1pass": "",
        "--presets": "",
        "--mem-flag": "",
        "--use-gpu": "",
        "--gpu-mem": "",
        "--keep-tmp-files": "",
        "--tmp-dir": "",
        "--continue": "",
        "--verbose": ""
    },
    "idba_ud": {
        "--min_contig": 1000,
        "threads": 7,
        "run": "",
        "--mink": "",
        "--maxk": "",
        "--step": "",
        "--inner_mink": "",
        "--inner_step": "",
        "--prefix": "",
        "--min_count": "",
        "--min_support": "",
        "--seed_kmer": "",
        "--similar": "",
        "--max_mismatch": "",
        "--min_pairs": "",
        "--no_bubble": "",
        "--no_local": "",
        "--no_coverage": "",
        "--no_correct": "",
        "--pre_correction": "",
        "use_scaffolds": ""
    },
    "iu_filter_quality_minoche": {
        "run": false,
        "--ignore-deflines": true,
        "--visualize-quality-curves": "",
        "--limit-num-pairs": "",
        "--print-qual-scores": "",
        "--store-read-fate": "",
        "threads": ""
    },
    "gzip_fastqs": {
        "run": false,
        "threads": ""
    },
    "bowtie": {
        "additional_params": "--no-unal",
        "threads": 8
    },
    "samtools_view": {
        "additional_params": "-F 4",
        "threads": ""
    },
    "anvi_profile": {
        "threads": 10,
        "--sample-name": "{sample}",
        "--overwrite-output-destinations": true,
        "--report-variability-full": "",
        "--skip-SNV-profiling": "",
        "--profile-SCVs": "",
        "--description": "",
        "--skip-hierarchical-clustering": "",
        "--distance": "",
        "--linkage": "",
        "--min-contig-length": "",
        "--min-mean-coverage": "",
        "--min-coverage-for-variability": "",
        "--cluster-contigs": "",
        "--contigs-of-interest": "",
        "--queue-size": "",
        "--write-buffer-size-per-thread": "",
        "--fetch-filter": "",
        "--min-percent-identity": "",
        "--max-contig-length": ""
    },
    "anvi_merge": {
        "--sample-name": "{group}",
        "--overwrite-output-destinations": true,
        "--description": "",
        "--skip-hierarchical-clustering": "",
        "--enforce-hierarchical-clustering": "",
        "--distance": "",
        "--linkage": "",
        "threads": ""
    },
    "import_percent_of_reads_mapped": {
        "run": true,
        "threads": ""
    },
    "krakenuniq": {
        "threads": 3,
        "--gzip-compressed": true,
        "additional_params": "",
        "run": "",
        "--db": ""
    },
    "remove_short_reads_based_on_references": {
        "delimiter-for-iu-remove-ids-from-fastq": " ",
        "dont_remove_just_map": "",
        "references_for_removal_txt": "",
        "threads": ""
    },
    "anvi_cluster_contigs": {
        "--collection-name": "{driver}",
        "run": "",
        "--driver": "",
        "--just-do-it": "",
        "--additional-params-concoct": "",
        "--additional-params-metabat2": "",
        "--additional-params-maxbin2": "",
        "--additional-params-dastool": "",
        "--additional-params-binsanity": "",
        "threads": ""
    },
    "gen_external_genome_file": {
        "threads": ""
    },
    "export_gene_calls_for_centrifuge": {
        "threads": ""
    },
    "anvi_import_taxonomy_for_genes": {
        "threads": ""
    },
    "annotate_contigs_database": {
        "threads": ""
    },
    "anvi_get_sequences_for_gene_calls": {
        "threads": ""
    },
    "gunzip_fasta": {
        "threads": ""
    },
    "reformat_external_gene_calls_table": {
        "threads": ""
    },
    "reformat_external_functions": {
        "threads": ""
    },
    "import_external_functions": {
        "threads": ""
    },
    "anvi_run_pfams": {
        "run": "",
        "--pfam-data-dir": "",
        "threads": ""
    },
    "iu_gen_configs": {
        "--r1-prefix": "",
        "--r2-prefix": "",
        "threads": ""
    },
    "gen_qc_report": {
        "threads": ""
    },
    "merge_fastqs_for_co_assembly": {
        "threads": ""
    },
    "merge_fastas_for_co_assembly": {
        "threads": ""
    },
    "bowtie_build": {
        "additional_params": "",
        "threads": ""
    },
    "anvi_init_bam": {
        "threads": ""
    },
    "krakenuniq_mpa_report": {
        "threads": ""
    },
    "import_krakenuniq_taxonomy": {
        "--min-abundance": "",
        "threads": ""
    },
    "anvi_summarize": {
        "additional_params": "",
        "run": "",
        "threads": ""
    },
    "anvi_split": {
        "additional_params": "",
        "run": "",
        "threads": ""
    },
    "references_mode": true,
    "all_against_all": "",
    "kraken_txt": "",
    "collections_txt": "",
    "output_dirs": {
        "FASTA_DIR": "02_FASTA",
        "CONTIGS_DIR": "03_CONTIGS",
        "QC_DIR": "01_QC",
        "MAPPING_DIR": "04_MAPPING",
        "PROFILE_DIR": "05_ANVIO_PROFILE",
        "MERGE_DIR": "06_MERGED",
        "TAXONOMY_DIR": "07_TAXONOMY",
        "SUMMARY_DIR": "08_SUMMARY",
        "SPLIT_PROFILES_DIR": "09_SPLIT_PROFILES",
        "LOGS_DIR": "00_LOGS"
    },
    "max_threads": "",
    "config_version": "3",
    "workflow_name": "metagenomics"
}
```
</details>

---

We, again, did a dry run
```bash
anvi-run-workflow -w metagenomics -c default-metagenomics-config.json -A  -n -q
```
You were probably curious about what the dry run stats look like, so I included this one to appease your curiosity. It reflected that we are working with 1850 metagenomes but only have one fasta file containing the reference genomes. So all good.
```bash
Job stats:
job                                  count
---------------------------------  -------
annotate_contigs_database                1
anvi_gen_contigs_database                1
anvi_init_bam                         1850
anvi_merge                               1
anvi_profile                          1850
anvi_run_hmms                            1
anvi_run_kegg_kofams                     1
anvi_run_ncbi_cogs                       1
anvi_run_scg_taxonomy                    1
bowtie                                1850
bowtie_build                             1
import_percent_of_reads_mapped        1850
metagenomics_workflow_target_rule        1
samtools_view                         1850
total                                 9259
```

As all looked well, we submitted the job again using clusterize:
```bash
clusterize -j metagenomics_workflow \
                 -o metagenomics_workflow.log \
                 -n 1 \
                 "anvi-run-workflow -w metagenomics \
    -c default-metagenomics-config.json   -A \
   --cluster 'clusterize -o {log} -n {threads} -j {rule} ' \
   --resources nodes=400 \
   --jobs 400 \
   --rerun-incomplete \
   --keep-going"
```



## Filtering samples based on coverage and detection

In this section, we describe all the steps taken to filter the samples based on coverage and detection. Our cutoff was be to only keep samples, in which at least one of the metagenomes was detected with at least 0.5 detection and 10x coverage.

The filtering included multiple substeps:
- Generating a genomic collection for `anvi-split`
- using `anvi-split` to create individual, self-contained anvi’o projects for each reference genome and its recruited reads
- using `anvi-profile-blitz` to profile the BAM files to get contig-level coverage and detection stats
- filtering the samples based on coverage and detection of SAR11 reference genomes

### Generating a genomic collection

At this point, anvi’o still did not know how to link scaffolds to each isolate genome. In anvi’o, this kind of knowledge is maintained through ‘collections’. In order to link scaffolds to genomes of origin, we used the program [`anvi-import-collection`](https://anvio.org/help/main/programs/anvi-import-collection/) to create an anvi’o collection in our merged profile database. This program, however, needed a [`collection-txt`](https://anvio.org/help/main/artifacts/collection-txt/) artifact. 

We made that by:

For this, we first created a `contigs.txt` (ran this in the `readRecruitment` directory)
```bash
grep ">" ../refGenomesSAR11/fastaDrep/dereplicated_genomes/all_fasta.fa | sed "s/>//g" > contigs.txt
```
This removed the carrot ">" from the contig headers in our fasta file full of reference genomes by substituting it with nothing globally (not just once but anywhere it sees a ">" in the file) and parsed the output into a file called `contigs.txt`. 

These are the names by which the contigs were stored within anvi'o. 

---

<details markdown="1"><summary>Click to show/hide the head of the file</summary>


```bash
FZCC0015_000000000001
HIMB058_000000000001
HIMB058_000000000002
HIMB058_000000000003
HIMB058_000000000004
HIMB058_000000000005
HIMB058_000000000006
HIMB058_000000000007
HIMB058_000000000008
HIMB058_000000000009
```

</details>

---

Next, we created a file called `bins.txt` which contained only the identifier for the genome.
```bash
grep ">" ../refGenomesSAR11/fastaDrep/dereplicated_genomes/all_fasta.fa | sed "s/>//g" | cut -d "_" -f 1 > bins.txt
```
This does the same as the command above PLUS it removes the suffix of "\_000000000001" from the reference genome identifier: it cuts at the delimiter (`-d`) "_" and keeps the first field (`-f 1`), so the part in front of the delimiter)

---

<details markdown="1"><summary>Click to show/hide the head of the file</summary>


```bash
FZCC0015
HIMB058
HIMB058
HIMB058
HIMB058
HIMB058
HIMB058
HIMB058
HIMB058
HIMB058
HIMB058
```

</details>

---

Lastly, we combined the two outputs into a single file and fed that into a new file called `collection.txt`, which is what we needed (the two steps before were just for us to get here). 
```bash
paste contigs.txt bins.txt > collection.txt
```
[`collection.txt`](https://anvio.org/help/main/artifacts/collection-txt/) is a two-column TAB-delimited file without a header that describes a [collection](https://anvio.org/help/main/artifacts/collection) by associating items with [bin](https://anvio.org/help/main/artifacts/bin) names.

---

<details markdown="1"><summary>Click to see the first 104 lines of the file</summary>


```bash
FZCC0015_000000000001   FZCC0015
HIMB058_000000000001    HIMB058
HIMB058_000000000002    HIMB058
HIMB058_000000000003    HIMB058
HIMB058_000000000004    HIMB058
HIMB058_000000000005    HIMB058
HIMB058_000000000006    HIMB058
HIMB058_000000000007    HIMB058
HIMB058_000000000008    HIMB058
HIMB058_000000000009    HIMB058
HIMB058_000000000010    HIMB058
HIMB058_000000000011    HIMB058
HIMB058_000000000012    HIMB058
HIMB058_000000000013    HIMB058
HIMB058_000000000014    HIMB058
HIMB058_000000000015    HIMB058
HIMB058_000000000016    HIMB058
HIMB058_000000000017    HIMB058
HIMB058_000000000018    HIMB058
HIMB058_000000000019    HIMB058
HIMB058_000000000020    HIMB058
HIMB058_000000000021    HIMB058
HIMB058_000000000022    HIMB058
HIMB058_000000000023    HIMB058
HIMB058_000000000024    HIMB058
HIMB058_000000000025    HIMB058
HIMB058_000000000026    HIMB058
HIMB058_000000000027    HIMB058
HIMB058_000000000028    HIMB058
HIMB058_000000000029    HIMB058
HIMB058_000000000030    HIMB058
HIMB058_000000000031    HIMB058
HIMB058_000000000032    HIMB058
HIMB058_000000000033    HIMB058
HIMB058_000000000034    HIMB058
HIMB058_000000000035    HIMB058
HIMB058_000000000036    HIMB058
HIMB058_000000000037    HIMB058
HIMB058_000000000038    HIMB058
HIMB058_000000000039    HIMB058
HIMB058_000000000040    HIMB058
HIMB058_000000000041    HIMB058
HIMB058_000000000042    HIMB058
HIMB058_000000000043    HIMB058
HIMB058_000000000044    HIMB058
HIMB058_000000000045    HIMB058
HIMB058_000000000046    HIMB058
HIMB058_000000000047    HIMB058
HIMB058_000000000048    HIMB058
HIMB058_000000000049    HIMB058
HIMB058_000000000050    HIMB058
HIMB058_000000000051    HIMB058
HIMB058_000000000052    HIMB058
HIMB058_000000000053    HIMB058
HIMB058_000000000054    HIMB058
HIMB058_000000000055    HIMB058
HIMB058_000000000056    HIMB058
HIMB058_000000000057    HIMB058
HIMB058_000000000058    HIMB058
HIMB058_000000000059    HIMB058
HIMB114_000000000001    HIMB114
HIMB123_000000000001    HIMB123
HIMB123_000000000002    HIMB123
HIMB123_000000000003    HIMB123
HIMB123_000000000004    HIMB123
HIMB123_000000000005    HIMB123
HIMB123_000000000006    HIMB123
HIMB123_000000000007    HIMB123
HIMB123_000000000008    HIMB123
HIMB123_000000000009    HIMB123
HIMB123_000000000010    HIMB123
HIMB123_000000000011    HIMB123
HIMB123_000000000012    HIMB123
HIMB123_000000000013    HIMB123
HIMB123_000000000014    HIMB123
HIMB123_000000000015    HIMB123
HIMB123_000000000016    HIMB123
HIMB123_000000000017    HIMB123
HIMB123_000000000018    HIMB123
HIMB123_000000000019    HIMB123
HIMB123_000000000020    HIMB123
HIMB123_000000000021    HIMB123
HIMB123_000000000022    HIMB123
HIMB123_000000000023    HIMB123
HIMB123_000000000024    HIMB123
HIMB123_000000000025    HIMB123
HIMB123_000000000026    HIMB123
HIMB123_000000000027    HIMB123
HIMB123_000000000028    HIMB123
HIMB123_000000000029    HIMB123
HIMB123_000000000030    HIMB123
HIMB123_000000000031    HIMB123
HIMB123_000000000032    HIMB123
HIMB123_000000000033    HIMB123
HIMB123_000000000034    HIMB123
HIMB123_000000000035    HIMB123
HIMB123_000000000036    HIMB123
HIMB123_000000000037    HIMB123
HIMB123_000000000038    HIMB123
HIMB123_000000000039    HIMB123
HIMB123_000000000040    HIMB123
HIMB123_000000000041    HIMB123
HIMB1402_000000000001   HIMB1402
HIMB1402_000000000002   HIMB1402
```

</details>

---

We then used the program [`anvi-import-collection`](https://anvio.org/help/main/programs/anvi-import-collection/) to import this collection into the anvi’o profile database by naming this collection `SAR11COLLECTION`:

```bash
anvi-import-collection collection.txt -p 06_MERGED/reference_genomes/PROFILE.db -c 03_CONTIGS/reference_genomes-contigs.db -C SAR11COLLECTION --contigs-mode
```

We used the `--contigs-mode` because our `collection.txt` describes contigs rather than splits.

**Note:** The collection is not stored in a separate file but in the `PROFILE.db`.

### Creating individual, self-contained anvi’o projects for each reference genome and its recruited reads

We used `anvi-split` to create individual, self-contained anvi’o projects for each reference genome and its recruited reads.

We, again, used clusterize to submit the job:
```bash
clusterize -j SAR11_split_job \
           -o SAR11_split_job.log \
           -n 8 \
           -M 96000 \
           --nodelist mpcs045 \
           "anvi-split -p 06_MERGED/reference_genomes/PROFILE.db \
                       -c 03_CONTIGS/reference_genomes-contigs.db \
                       -C SAR11COLLECTION \
                       -o SAR11SPLIT"
```

The output was put into a directory called `SAR11SPLIT`, which then contained individual directories with `contig.db` and `profile.db` for each reference genome and recruited reads. 

### Determining which samples to continue working with


#### Getting contig-level coverage and detection stats
We used [`anvi-profile-blitz`](https://anvio.org/help/main/programs/anvi-profile-blitz/) to find out which samples met our filtering criteria if: my genomes should be covered 10x and detected 0.5). `anvi-profile-blitz` allows the fast profiling of BAM files to get contig- or gene-level coverage and detection stats.

We gave it the `BAM file`s that were created in the metagenomics workflow and were stored in the `04_MAPPING` directory, as well as the `CONTIGS.db`, and specified what the output should look like.

This is the base command
```bash
anvi-profile-blitz *.bam \
                   -c contigs-db \
                   -o OUTPUT.txt
```

However, since we did it for each `split`, so 51 times (once per reference genome), we adapted the base command a bit to do it in a loop.

```bash
nano run_anvi_profile_blitz.sh
```
content
```bash
#!/bin/bash

# Directory containing the subdirectories
BASE_DIR="SAR11SPLIT"
# Directory for the output files
OUTPUT_DIR="blitzOUTPUT"

# Create the output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Loop through each subdirectory in BASE_DIR (checks if the current item is a directory)
for DIR in $BASE_DIR/*; do
  if [ -d "$DIR" ]; then
    # Extract the directory name from the full path (to later name the output file)
    DIR_NAME=$(basename $DIR)
    # Construct the command
    anvi-profile-blitz 04_MAPPING/reference_genomes/*.bam \
                       -c $DIR/CONTIGS.db \
                       -o $OUTPUT_DIR/blitzOUTPUT_$DIR_NAME.txt
  fi
done
```
make it executable
```bash
chmod +x run_anvi_profile_blitz.sh
```
submit job
```bash
clusterize "./run_anvi_profile_blitz.sh" -n 1 -j anvi_profile_blitz_job
```

#### Combining BLITZ outputs into one file
We combined the BLITZ outputs into one dataframe. You can get the output [here](files/combined_blitzOUTPUT.txt.zip).
```bash
cd blitzOUTPUT/
nano combineOutputs.py
```
content
```python
import pandas as pd
import glob
import os

# Directory containing the output files
output_dir = "."

# List to hold individual data frames
dfs = []

# Read each file and append to the list
files = glob.glob(os.path.join(output_dir, "*.txt"))
for file_path in files:
    df = pd.read_csv(file_path, sep="\t")
    reference_genome = os.path.basename(file_path).replace("blitzOUTPUT_", "").replace(".txt", "")
    df['reference_genome'] = reference_genome
    
    # Extract project prefix from sample names (assuming prefix is before the first underscore)
    df['project'] = df['sample'].str.split('_').str[0]
    
    dfs.append(df)

# Combine all data frames
combined_blitzOUTPUT = pd.concat(dfs, ignore_index=True)

# Save the combined data frame as a .txt file
combined_blitzOUTPUT.to_csv('./combined_blitzOUTPUT.txt', sep='\t', index=False)

```
run
```bash
python combineOutputs.py
```


#### Getting weighted averages of coverage and detection
Even though some reference genomes had multiple contigs, we wanted to know how much of each reference genome was covered in a given sample / how much the entire reference genome was detected (so across contigs): so the mean_cov and detection information for the collective reference genome. However, we could not just take the information of all contigs from the same reference genome across a given sample and take the average but had to take the length of each contig compared to the length of the entire reference genome (length of its contigs combined) into consideration. 

To calculate the weighted averages of both `mean_cov` and `detection`, we wrote the following script.
```bash
nano calculate_weighted_coverage.py
```
content
```python
import pandas as pd

# Load the data into a pandas DataFrame
df = pd.read_csv('combined_blitzOUTPUT.txt', delim_whitespace=True)

# Calculate total length for each reference genome within each sample
df['total_length'] = df.groupby(['sample', 'reference_genome'])['length'].transform('sum')

# Calculate weighted mean coverage and detection
df['weighted_mean_cov'] = df['mean_cov'] * df['length'] / df['total_length']
df['weighted_detection'] = df['detection'] * df['length'] / df['total_length']

# Group by sample, reference_genome, and project to sum the weighted values
result = df.groupby(['sample', 'reference_genome', 'project']).agg({
    'weighted_mean_cov': 'sum',
    'weighted_detection': 'sum',
    'total_length': 'first'  # All values in the group should be the same, so just take the first one
}).reset_index()

# Save the results to a new tab-separated .txt file
result.to_csv('weighted_results_with_length.txt', sep='\t', index=False)
```
run
```bash
calculate_weighted_coverage.py
```

The output file is available [here](files/weighted_results_with_length.txt.zip) (zipped).


#### Getting overview stats of how many samples per project would stay with different detection and coverage cut-offs

How many samples per project would make it if we only took the samples for which at least one reference genome has at least 10x coverage and at least 0.5 detection (both need to be true for a sample to pass)? 

This is what we wanted to know for different coverage and detection combos. To find out we wrote the following script.
```bash
nano filter_samples_all_combinations.py
```
content
```python
import pandas as pd

# Load the data from the specified path
file_path = './weighted_results_with_length.txt'
combined_df = pd.read_csv(file_path, sep='\t')

# Calculate the number of samples per project before applying any filters
initial_stats = combined_df.groupby('project')['sample'].nunique().reset_index()
initial_stats.columns = ['project', 'num_samples_initial']

# List of coverage and detection thresholds
coverage_thresholds = list(range(1, 11))  # 1x to 10x
detection_thresholds = [0.25, 0.3, 0.4, 0.5]

# Initialize a DataFrame to hold all statistics
all_stats = initial_stats.copy()

# Function to generate a filter criteria string
def generate_filter_criteria(cov, det):
    return f'{cov}x_{int(det*100)}'

# Iterate through all combinations of coverage and detection thresholds
for cov in coverage_thresholds:
    for det in detection_thresholds:
        # Apply the filtering criteria
        filtered_df = combined_df[(combined_df['weighted_mean_cov'] >= cov) & (combined_df['weighted_detection'] >= det)]
        # Calculate the number of samples per project for the current criteria
        stats_filtered = filtered_df.groupby('project')['sample'].nunique().reset_index()
        filter_criteria = generate_filter_criteria(cov, det)
        stats_filtered.columns = ['project', f'num_samples_{filter_criteria}']
        # Merge the statistics with the all_stats DataFrame
        all_stats = pd.merge(all_stats, stats_filtered, on='project', how='outer')

# Save the combined statistics to a single file
output_combined_stats_path = './combined_stats_weighted_all_filters.txt'
all_stats.to_csv(output_combined_stats_path, sep='\t', index=False)

# Print the combined statistics
print("Combined statistics for all filtering criteria:")
print(all_stats)
```
run
```bash
python filter_samples_all_combinations.py
```

This gave us the `combined_stats_weighted_all_filters.txt` file, which is available [here](files/combined_stats_weighted_all_filters.txt) and gives the number of samples ber project given different detection and coverage cut-offs.

```bash
project num_samples_initial     num_samples_1x_25       num_samples_1x_30       num_samples_1x_40       num_samples_1x_50       num_samples_2x_25       num_samples_2x_30       num_samples_2x_40       num_samples_2x_50       num_samples_3x_25       num_samples_3x_30       num_samples_3x_40       num_samples_3x_50       num_samples_4x_25       num_samples_4x_30       num_samples_4x_40       num_samples_4x_50       num_samples_5x_25       num_samples_5x_30       num_samples_5x_40       num_samples_5x_50       num_samples_6x_25       num_samples_6x_30       num_samples_6x_40       num_samples_6x_50       num_samples_7x_25       num_samples_7x_30       num_samples_7x_40       num_samples_7x_50       num_samples_8x_25       num_samples_8x_30       num_samples_8x_40       num_samples_8x_50       num_samples_9x_25       num_samples_9x_30       num_samples_9x_40       num_samples_9x_50       num_samples_10x_25      num_samples_10x_30      num_samples_10x_40      num_samples_10x_50
BATS    40      40      40      40      40      40      40      40      40      38      38      38      38      36      36      36      36      35      35      35      35      31      31      31      31      29      29      29      29      21      21      21      21      17      17      17      17      16      16      16      16
BGS     969     941     941     938     907     828     828     828     828     691     691     691     691     579     579     579     579     485     485     485     485     408     408     408     408     335     335     335     335     286     286     286     286     236     236     236     236     201     201     201     201
BGT     323     313     313     313     310     309     309     308     306     291     291     291     291     271     271     271     271     231     231     231     231     193     193     193     193     170     170     170     170     149     149     149     149     127     127     127     127     104     104     104     104
HOT1    28      28      28      28      28      28      28      28      28      26      26      26      26      25      25      25      25      23      23      23      23      22      22      22      22      22      22      22      22      22      22      22      22      19      19      19      19      17      17      17      17
HOT3    230     182     182     182     182     160     160     160     160     151     151     151     150     140     140     140     140     131     131     131     131     120     120     120     120     104     104     104     104     92      92      92      92      86      86      86      86      74      74      74      74
MAL     16      16      16      16      16      16      16      16      16      16      16      16      16      16      16      16      16      16      16      16      16      16      16      16      16      16      16      16      16      16      16      16      16      16      16      16      16      16      16      16      16
OSD     149     82      82      68      44      48      48      48      43      29      29      29      28      22      22      22      22      14      14      14      14      9       9       9       9       7       7       7       7       6       6       6       6       5       5       5       5       4       4       4       4
TARA    95      95      94      94      94      95      94      94      94      94      94      94      94      94      94      94      94      94      94      94      94      94      94      94      94      94      94      94      94      94      94      94      94      94      94      94      94      94      94      94      94
```


#### Filtering based on 0.5 detection and 10x coverage

Knowing the above, we set our filtering cutoff to 10x coverage and 0.5 detection. We wrote the following script to get the number of samples per project that pass that filtering  threshold, as well as their coverage and detection values.
```bash
nano filter_samples_detNcov.py
```
content
```python
import pandas as pd

# Load the data from the specified path
file_path = './weighted_results_with_length.txt'
combined_df = pd.read_csv(file_path, sep='\t')

# Calculate the number of samples per project before applying any filters
initial_stats = combined_df.groupby('project')['sample'].nunique().reset_index()
initial_stats.columns = ['project', 'num_samples_initial']

# Filter for samples with at least one reference genome having at least 0.5 weighted detection
filter_05_detection = combined_df[combined_df['weighted_detection'] >= 0.5]

# Filter for samples with at least one reference genome having at least 10x weighted mean coverage
filter_10x = combined_df[combined_df['weighted_mean_cov'] >= 10]

# Combine filters to get samples that passed both criteria
combined_filters = filter_05_detection.merge(filter_10x, on=['sample', 'total_length', 'weighted_detection', 'weighted_mean_cov', 'reference_genome', 'project'])

# Generate statistics for 0.5 weighted detection filter
stats_05_detection = filter_05_detection.groupby('project')['sample'].nunique().reset_index()
stats_05_detection.columns = ['project', 'num_samples_05_detection']

# Generate statistics for 10x weighted mean coverage filter
stats_10x = filter_10x.groupby('project')['sample'].nunique().reset_index()
stats_10x.columns = ['project', 'num_samples_10x']

# Print the statistics
print("Initial statistics for samples per project:")
print(initial_stats)

print("\nStatistics for samples with at least 0.5 weighted detection:")
print(stats_05_detection)

print("\nStatistics for samples with at least 10x weighted mean coverage:")
print(stats_10x)

# Merge all statistics into a single DataFrame
all_stats = initial_stats.merge(stats_05_detection, on='project', how='outer').merge(stats_10x, on='project', how='outer')

# Save all statistics to a single .txt file
output_stats_path = './filtered05N10x_stats.txt'
all_stats.to_csv(output_stats_path, sep='\t', index=False)

# Save the combined filtered DataFrame to a .txt file
output_filtered_path = './filtered05N10x_combined_df.txt'
combined_filters.to_csv(output_filtered_path, sep='\t', index=False)
```
run
```bash
python filter_samples_detNcov.py
```

This resulted in two files called `filtered05N10x_stats.txt` (containing the number of samples of each project that passed the filter criteria, available [here](files/filtered05N10x_stats.txt)) and `filtered05N10x_combined_df.txt` (the sample_IDs, contigs, weighted_mean_cov, weighted_detection of the samples that passed the filtering, available [here](files/filtered05N10x_combined_df.txt)).

### Generating anvi'o projects for samples that passed the filtering

Based on the information we had post-detection0.5-coverage10x-filtering, we knew which samples to continue working with. To continue working with them, however, we needed to generate a new anvi'o `PROFILE.db` with those samples only.

We created a .txt file that only lists the names of the samples that passed (and only lists each once).
```bash
# Extract the first column (sample names), sort them, and keep only unique entries
awk 'NR>1 {print $1}' ./blitzOUTPUT/filtered05N10x_combined_df.txt | sort | uniq > ./blitzOUTPUT/unique_samples.txt

```

#### Profiling the mapping results with anvi’o
To profile the mapping results with anvi'o, we generated individual `PROFILE.db`s for each `.bam` file of the samples we want. The `.bam` files were stored in `./04_MAPPING/reference_genomes/`.

To generate these `profile.db`s we used anvi'o's [`anvi-profile`](https://anvio.org/help/7.1/programs/anvi-profile/) program. In its simplest form, it looks like this:
```bash
anvi-profile -i SAMPLE-01.bam -c contigs.db
```

However, we adapted it to our needs as such: Submit multiple `anvi-profile` jobs to a cluster using the `clusterize` tool. The script reads a list of sample names from a file, processes the corresponding BAM files, and submits the jobs in a controlled manner to avoid overwhelming the cluster scheduler. It includes functionality to stagger job submissions and limit the number of jobs running simultaneously.

We kept the number of threads = 16, as described here: [https://merenlab.org/2017/03/07/the-new-anvio-profiler/](https://merenlab.org/2017/03/07/the-new-anvio-profiler/).

```bash
nano run_anvi_profile_cluster.sh
```
content
```bash
#!/bin/bash

# Maximum number of jobs to submit at once
max_jobs=10

# Iterate over each line in the unique_samples.txt file
while IFS=$'\t' read -r sample
do
    # Define the path to the BAM file
    bam_file="./04_MAPPING/reference_genomes/${sample}.bam"

    # Define the output directory
    output_dir="./filteredPROFILEdb/$sample"

    # Define the log file for each sample
    log_file="./filteredPROFILEdb/${sample}_anvi_profile.log"

    # Check if the BAM file exists
    if [ -f "$bam_file" ]; then
        # Define the command to run anvi-profile with 16 threads
        cmd="anvi-profile -c 03_CONTIGS/reference_genomes-contigs.db \
                          -i $bam_file \
                          --num-threads 16 \
                          -o $output_dir"

        # Submit the job using clusterize
        clusterize -j "anvi_profile_$sample" \
                   -o "$log_file" \
                   -n 16 \
                   --mail-user raissa.meyer@awi.de \
                   "$cmd"

        # Check the number of jobs in the queue, and wait if too many are running
        while [ $(squeue -u $USER | wc -l) -ge $max_jobs ]; do
            sleep 60
        done

    else
        echo "BAM file for sample $sample not found, skipping..."
    fi

done < ./blitzOUTPUT/unique_samples.txt
```
give it permission
```bash
chmod +x run_anvi_profile_cluster.sh
```
and run
```bash
./run_anvi_profile_cluster.sh
```

{:.notice}
Following this, there are multiple avenues one can take. We will describe the one that led to the figure seen in the manuscript in detail. The other avenues we took on the way to get there were 1) generating anvi'o projects such that all samples that had any reference genome at the required ≥10x cov and ≥0.5 det were included and 2) generating anvi'o projects such that each reference genome project includes only those samples where the specific reference genome for that project was found at ≥10x cov and ≥0.5 det. They will be described a bit more briefly first (mainly, the scripts will all be collapsed but can be extended if you are interested).

#### Generating anvi'o projects for all samples that passed filtering (inclusive) 

Here, each reference genome project included all samples that had **any of the 51 reference genomes at ≥10x coverage and ≥0.5 detection, even if the sample met these thresholds for a different reference genome than the one the project is based on**.

For this, we first generated a merged anvi’o profile db that combines all sample-level `profile-db`s (of our selected samples) into a single `profile-db` using [`anvi-merge`](https://anvio.org/help/7.1/programs/anvi-merge/). Basically, this takes the alignment data from each sample (each contained in its own sample-level `profile-db`) and combines them into a single database that anvi’o can look through more easily.

---

<details markdown="1"><summary>Click to see how we did that</summary>

The basic command is:
```bash
anvi-merge -c cool_contigs.db \
            Single_profile_db_1 Single_profile_db_2 \
            -o cool_contigs_merge
```

We adjusted it to our directory structure and naming and ran it in `filteredPROFILEdb/`
```bash
clusterize -j merge_profile_db \
-o merge_profile_db.log \
-n 1 \
"anvi-merge */PROFILE.db -o SAR11-MERGED -c ../03_CONTIGS/reference_genomes-contigs.db"
```

</details>

---

We then created self-contained anvi'o projects for our 51 reference genomes and associated metagenomes using `anvi-split` on this merged `PROFILE.db`, to which we first added a collection again. This gave us an anvi'o project per reference genome.

---

<details markdown="1"><summary>Click to see how we did that</summary>
    
First, we had to re-create a COLLECTION that will be stored in the PROFILE.db. For that, we reused the `collection.txt` that we made before, but had to re-create the COLLECTION in itself: specifying the new PROFILE.db.

Here we go:
```bash
anvi-import-collection collection.txt \
   -p filteredPROFILEdb/SAR11-MERGED/PROFILE.db \
   -c 03_CONTIGS/reference_genomes-contigs.db \
   -C SAR11COLLECTION \
   --contigs-mode
```

Then the actual creation of the individual, self-contained anvi'o projects for each reference genome and its recruited reads can start. 

The output went into a directory called `SAR11_postFilter`.

```bash
clusterize -j SAR11_split_post_filter_job \
           -o SAR11_split_post_filter_job.log \
           -n 8 \
           -M 96000 \
           --nodelist mpcs052 \
           "anvi-split -p filteredPROFILEdb/SAR11-MERGED/PROFILE.db \
                       -c 03_CONTIGS/reference_genomes-contigs.db \
                       -C SAR11COLLECTION \
                       -o SAR11SPLIT_postFilter"
```

</details>

---

#### Generating anvi'o projects for all samples that passed filtering (exclusive)

Here, each reference genome project included only the samples **where the specific reference genome for that project was found at ≥10x coverage and ≥0.5 detection**.

Instead of creating a merged `PROFILE.db` with ALL samples and then splitting it into the different reference genomes later, we created individual merged `PROFILE.db`s, one per reference genome. We just had to know exactly which samples should be added for which reference genomes. We could get that information from the `weighted_results_with_legth.txt` that we got as part of the filtering above. We still have to perform a split at the end, but just to get the proper `CONTIGS.db`


---

<details markdown="1"><summary>Click to see how we did that</summary>

To merge the sample-level `PROFILE.db`s into reference-genome-level-merged `PROFILE.db`s, we wrote this script:
```
#!/bin/bash

# Loop through all unique reference genomes
for genome in $(cut -f2 ../blitzOUTPUT/weighted_results_with_length.txt | tail -n +2 | sort | uniq); do
    # Step 1: Create the directory if it doesn't exist
    mkdir -p "mergeSpecifics/${genome}-profiles"

    # Step 2: Extract samples where the conditions (0.5 det and 10x cov) are met
    awk -v ref_genome="$genome" '$2 == ref_genome && $4 >= 10 && $5 >= 0.5 {print $1}' ../blitzOUTPUT/weighted_results_with_length.txt | sort | uniq > "mergeSpecifics/${genome}-profiles/matching_samples_${genome}.txt"

    # Step 3: Remove land samples if they exist
    for sample in SAMEA2619802 SAMEA3275469 SAMEA3275484 SAMEA3275486 SAMEA3275509 SAMEA3275519 SAMEA3275520 SAMEA3275533 SAMEA3275581 SAMEA3275586; do
        sed -i "/_${sample}_/d" "mergeSpecifics/${genome}-profiles/matching_samples_${genome}.txt"
    done

    # Step 4: Collect all PROFILE.db files in a variable
    # Restrict the search to the current directory only
    profile_dbs=$(while read -r sample; do
        find . -maxdepth 1 -type d -name "$sample" -exec find {} -name "PROFILE.db" \;
    done < "mergeSpecifics/${genome}-profiles/matching_samples_${genome}.txt" | tr '\n' ' ')

    # Step 5: Run the anvi-merge command with contigs-db and collected PROFILE.db files
    clusterize -j merge_profile_db \
        -o mergeSpecifics/${genome}_merge_profile_db.log \
        -n 1 \
        "anvi-merge $profile_dbs -o mergeSpecifics/${genome}-profiles/${genome}-MERGED -c ../03_CONTIGS/reference_genomes-contigs.db"

done
```

checked for errors
```
$ grep -i "error" *.log
HIMB1412_merge_profile_db_stdyhADzTp.log:anvi-merge: error: the following arguments are required: SINGLE_PROFILE(S)
HIMB1420_merge_profile_db_sdAfWZpUKt.log:anvi-merge: error: the following arguments are required: SINGLE_PROFILE(S)
LSUCC0530_merge_profile_db_nvSystbihu.log:anvi-merge: error: the following arguments are required: SINGLE_PROFILE(S)
```
These reference genomes did not work because none of the samples passed the filter criteria of ≥10x coverage and ≥ 0.5 detection. I removed their directories so as to not get confused later.
```
rm -r  HIMB1412-profiles/  HIMB1420-profiles/ LSUCC0530-profiles/
```

We built the collection and do splits with the following script

```
#!/bin/bash

# Loop through all unique reference genomes
for genome in $(cut -f2 ../blitzOUTPUT/weighted_results_with_length.txt | tail -n +2 | sort | uniq); do
    
    # Step 1: Create collection file and import collection
    grep "$genome" ../collection.txt > "mergeSpecifics/${genome}-profiles/collection_${genome}.txt"
    
    anvi-import-collection "mergeSpecifics/${genome}-profiles/collection_${genome}.txt" \
        -p "mergeSpecifics/${genome}-profiles/${genome}-MERGED/PROFILE.db" \
        -c ../03_CONTIGS/reference_genomes-contigs.db \
        -C "${genome}COLLECTION" \
        --contigs-mode

    # Step 2: Run anvi-split based on the collection
    clusterize -j SAR11_split_post_filter_${genome}_job \
        -o SAR11_split_post_filter_${genome}_job.log \
        -n 8 \
        -M 96000 \
        "anvi-split -p mergeSpecifics/${genome}-profiles/${genome}-MERGED/PROFILE.db \
                    -c ../03_CONTIGS/reference_genomes-contigs.db \
                    -C ${genome}COLLECTION \
                    -o mergeSpecifics/${genome}-profiles/${genome}_postFilter"

done
```
We checked for errors again, and as expected, the directories I removed in the previous step came up. But that is no problem.
```
$ grep -i "error" *.log
SAR11_split_post_filter_HIMB1412_job_hZBIlRsgLA.log:File/Path Error: No such file: 'mergeSpecifics/HIMB1412-profiles/HIMB1412-MERGED/PROFILE.db' :/
SAR11_split_post_filter_HIMB1420_job_oONSdvODyq.log:File/Path Error: No such file: 'mergeSpecifics/HIMB1420-profiles/HIMB1420-MERGED/PROFILE.db' :/
SAR11_split_post_filter_LSUCC0530_job_qSABXiLKqv.log:File/Path Error: No such file: 'mergeSpecifics/LSUCC0530-profiles/LSUCC0530-MERGED/PROFILE.db' :/
```

</details>

---

#### Generating anvi'o projects for all samples that passed filtering for reference genomes HTCC1002 and HIMB2204 (as done for the manuscript figure) 

{:.notice}
This and the following steps will actually get you to the figure in the manuscript. So let's go!

While above, we created anvi'o projects to inspect each reference genome across the different samples, for the final visualization used in the manuscript, we focused on the `HTCC1002` and `HIMB2204` reference genomes.

Thus, we wanted anvi'o projects that only included the samples relevant to them. For that, we first had to identify which those samples were. To do so, we used the following script:

```bash
#!/bin/bash

# Define the genomes of interest
genomes=("HIMB2204" "HTCC1002")

# Step 1: Create the directory for combined genomes if it doesn't exist
mkdir -p "mergeSpecifics/HIMB2204_HTCC1002-profiles"

# Step 2: Extract samples where the conditions are met for both HIMB2204 and HTCC1002
awk '($2 == "HIMB2204" || $2 == "HTCC1002") && $4 >= 10 && $5 >= 0.5 {print $1}' ../blitzOUTPUT/weighted_results_with_length.txt | sort | uniq > "mergeSpecifics/HIMB2204_HTCC1002-profiles/matching_samples_HIMB2204_HTCC1002.txt"

# Step 3: Remove land samples if they exist
for sample in SAMEA2619802 SAMEA3275469 SAMEA3275484 SAMEA3275486 SAMEA3275509 SAMEA3275519 SAMEA3275520 SAMEA3275533 SAMEA3275581 SAMEA3275586; do
    sed -i "/_${sample}_/d" "mergeSpecifics/HIMB2204_HTCC1002-profiles/matching_samples_HIMB2204_HTCC1002.txt"
done

# Step 4: Collect all PROFILE.db files in a variable
profile_dbs=$(while read -r sample; do
    find . -maxdepth 1 -type d -name "$sample" -exec find {} -name "PROFILE.db" \;
done < "mergeSpecifics/HIMB2204_HTCC1002-profiles/matching_samples_HIMB2204_HTCC1002.txt" | tr '\n' ' ')

# Step 5: Run the anvi-merge command with contigs-db and collected PROFILE.db files
if [ -n "$profile_dbs" ]; then
    clusterize -j merge_profile_db \
        -o mergeSpecifics/HIMB2204_HTCC1002_merge_profile_db.log \
        -n 1 \
        --nodelist mpcs043 \
        "anvi-merge $profile_dbs -o mergeSpecifics/HIMB2204_HTCC1002-profiles/HIMB2204_HTCC1002-MERGED -c ../03_CONTIGS/reference_genomes-contigs.db"
else
    echo "No PROFILE.db files found for HIMB2204 and HTCC1002"
fi
```

Then, we created a collection with HTCC1002 and HIMB2204 samples and added it to their collective `PROFILE.db`. 

For this, we started from the `collection.txt` we have been using throughout this workflow, but selected only the HIMB2204 and HTCC1002 entries:
```
grep -E "HIMB2204|HTCC1002" collection.txt > collection_HIMB2204_HTCC1002.txt
```
This resulted in the following file `collection_HIMB2204_HTCC1002.txt`. 
```
HIMB2204_000000000001   HIMB2204
HIMB2204_000000000002   HIMB2204
HTCC1002_000000000001   HTCC1002
HTCC1002_000000000002   HTCC1002
HTCC1002_000000000003   HTCC1002
HTCC1002_000000000004   HTCC1002
```
This is what we then added to the `PROFILE-db` as a new collection.
```
anvi-import-collection collection_HIMB2204_HTCC1002.txt \
    -p filteredPROFILEdb/mergeSpecifics/HIMB2204_HTCC1002-profiles/HIMB2204_HTCC1002-MERGED/PROFILE.db \
    -c 03_CONTIGS/reference_genomes-contigs.db \
    -C HIMB2204_HTCC1002_COLLECTION \
    --contigs-mode
```

Once this was done, we used the `anvi-split` program again to get separate anvi'o projects for HTCC1002 and HIMB2204.
```
clusterize -j SAR11_split_post_filter_HIMB2204nHTCC1002_job \
           -o SAR11_split_post_filter_HIMB2204nHTCC1002_job.log \
           -n 8 \
           -M 96000 \
           --nodelist mpcs052 \
           -t 04:00:00 \
           "anvi-split -p filteredPROFILEdb/mergeSpecifics/HIMB2204_HTCC1002-profiles/HIMB2204_HTCC1002-MERGED/PROFILE.db \
                       -c 03_CONTIGS/reference_genomes-contigs.db \
                       -C HIMB2204_HTCC1002_COLLECTION \
                       -o HIMB2204nHTCC1002_postFilter"
```
And checked for errors in the log file.
```
grep -i "error" SAR11_split_post_filter_HIMB2204nHTCC1002_job_VibnNEZDgI.log
```

## Visualising anvi'o phylograms

In this section, we will explain how we used `anvi-interactive` to visualize our metagenomes (Panels C and D of Figure 1 of the manuscript).

For `anvi-interactive` to give us what we want, we needed
- to decide which reference genomes to focus on (HIMB2203 and HTCC1002)
- a collection and bin to feed to `--gene-mode`
- to prepare and import the metadata we want to order our layers (metagenomes) by in the interactive interface

### Deciding which SAR11 reference genomes to visualise

To decide which SAR11 reference genomes to focus on, we checked which are found across most of the projects.

We wrote a little script to output which projects each reference genome is found in:

```bash
nano countSamplesEachRefGenome.py
```
content
```python
import pandas as pd

# Load the data from the txt file into a DataFrame
df = pd.read_csv('filtered05N10x_combined_df.txt', sep='\t')

# Group by 'reference_genome' and aggregate the unique projects into a list
genome_project_details = df.groupby('reference_genome').agg({
    'project': lambda x: list(x.unique())
}).reset_index()

# Count the number of unique projects per reference genome
genome_project_details['num_projects'] = genome_project_details['project'].apply(len)

# Sort by the number of projects in descending order
sorted_genomes = genome_project_details.sort_values(by='num_projects', ascending=False)

# Determine the total number of unique projects
total_projects = df['project'].nunique()

# Check which genomes are found in all projects
genomes_in_all_projects = sorted_genomes[sorted_genomes['num_projects'] == total_projects]

# Display the results
print("Genomes found in all projects:")
print(genomes_in_all_projects)

print("\nNumber of projects each genome is found in and the project names (sorted):")
print(sorted_genomes)
```
output
```bash
Genomes found in all projects:
Empty DataFrame
Columns: [reference_genome, project, num_projects]
Index: []

Number of projects each genome is found in and the project names (sorted):
   reference_genome                                  project  num_projects
0          FZCC0015  [BATS, BGS, BGT, HOT1, HOT3, MAL, TARA]             7
13         HIMB1483  [BATS, BGS, BGT, HOT1, HOT3, MAL, TARA]             7
46             RS39  [BATS, BGS, BGT, HOT1, HOT3, MAL, TARA]             7
37           HIMB83  [BATS, BGS, BGT, HOT1, HOT3, MAL, TARA]             7
34         HIMB2305  [BATS, BGS, BGT, HOT1, HOT3, MAL, TARA]             7
33         HIMB2250  [BATS, BGS, BGT, HOT1, HOT3, MAL, TARA]             7
32         HIMB2226  [BATS, BGS, BGT, HOT1, HOT3, MAL, TARA]             7
31         HIMB2204  [BATS, BGS, BGT, HOT1, HOT3, MAL, TARA]             7
30         HIMB2201  [BATS, BGS, BGT, HOT1, HOT3, MAL, TARA]             7
29         HIMB2200  [BATS, BGS, BGT, HOT1, HOT3, MAL, TARA]             7
28         HIMB2187  [BATS, BGS, BGT, HOT1, HOT3, MAL, TARA]             7
27         HIMB2104  [BATS, BGS, BGT, HOT1, HOT3, MAL, TARA]             7
47             RS40  [BATS, BGS, BGT, HOT1, HOT3, MAL, TARA]             7
12         HIMB1456  [BATS, BGS, BGT, HOT1, HOT3, MAL, TARA]             7
7          HIMB1413  [BATS, BGS, BGT, HOT1, HOT3, MAL, TARA]             7
10         HIMB1437  [BATS, BGS, BGT, HOT1, HOT3, MAL, TARA]             7
11         HIMB1444  [BATS, BGS, BGT, HOT1, HOT3, MAL, TARA]             7
6          HIMB1409  [BATS, BGS, BGT, HOT1, HOT3, MAL, TARA]             7
5          HIMB1402  [BATS, BGS, BGT, HOT1, HOT3, MAL, TARA]             7
3           HIMB123  [BATS, BGS, BGT, HOT1, HOT3, MAL, TARA]             7
9          HIMB1436  [BATS, BGS, BGT, HOT1, HOT3, MAL, TARA]             7
1           HIMB058        [BATS, BGS, BGT, HOT3, MAL, TARA]             6
39         HTCC7211              [BATS, BGS, BGT, MAL, TARA]             5
40         HTCC7214              [BATS, BGS, BGT, MAL, TARA]             5
45              NP1               [BGS, BGT, MAL, OSD, TARA]             5
38         HTCC1002               [BGS, BGT, MAL, OSD, TARA]             5
41         HTCC8051                    [BGS, BGT, MAL, TARA]             4
42         HTCC9022                    [BGS, BGT, MAL, TARA]             4
43         HTCC9565                    [BGS, BGT, MAL, TARA]             4
4           HIMB140                         [BGT, MAL, TARA]             3
24         HIMB1685                         [BGT, MAL, TARA]             3
14         HIMB1485                         [BGT, MAL, TARA]             3
17         HIMB1509                         [BGT, MAL, TARA]             3
25         HIMB1695                         [BGT, MAL, TARA]             3
23         HIMB1573                         [BGT, MAL, TARA]             3
16         HIMB1495                         [BGT, MAL, TARA]             3
21         HIMB1559                         [BGT, MAL, TARA]             3
26         HIMB1709                         [BGT, MAL, TARA]             3
15         HIMB1488                              [BGT, TARA]             2
44         IMCC9063                              [BGS, TARA]             2
20         HIMB1556                              [BGT, TARA]             2
18         HIMB1517                                   [TARA]             1
19         HIMB1520                                   [TARA]             1
22         HIMB1564                                   [TARA]             1
36            HIMB5                                   [TARA]             1
35            HIMB4                                   [TARA]             1
2           HIMB114                                   [TARA]             1
8          HIMB1427                                   [TARA]             1
```

{:.notice}
This showed that OSD was only included in NP1 and HTCC1002. In the manuscript, we went with HTCC1002, as this was found in both OSD samples, while NP1 was only found in one OSD sample.
HIMB2204 on the other hand is one of the reference genomes that was found in metagenomes from many different (7 of the 8) projects.


We also visualised the occurrence (above 10x cov and 0.5 detection) of reference genomes across samples in R

```R
## visualise anvi-profile-blitz output

library(ggplot2)
library(dplyr)

# Load data
blitz_covNdet <- read.table("/Users/rameyer/Documents/_P3/P3dataAnalysis/P3_visualise/outputBLITZ/filtered05N10x_combined_df.txt", header = TRUE, sep = "\t")
head(blitz_covNdet)

# Convert reference_genome to a factor
blitz_covNdet$reference_genome <- as.factor(blitz_covNdet$reference_genome)

# Define the custom color palette for projects
custom_palette <- c(
  "BATS" = "#2B5D9C", "BGS" = "#B13A8C", "BGT" = "#40B2B2", "HOT1" = "#875D9B",
  "HOT3" = "#f15c40", "MAL" = "#fdb64e", "OSD" = "#fce84e", "TARA" = "#469F77"
)

# Create a dataframe for project labels
project_labels <- blitz_covNdet %>%
  select(sample, project) %>%
  distinct() %>%
  arrange(sample)

# Assign colors from the custom palette to projects
project_colors <- project_labels %>%
  distinct(project) %>%
  mutate(color = custom_palette[seq_len(n())])

# Merge project colors back to project labels
project_labels <- project_labels %>%
  left_join(project_colors, by = "project")

# Ensure samples are factors with levels matching the order in project_labels
blitz_covNdet$sample <- factor(blitz_covNdet$sample, levels = project_labels$sample)

# Create a named vector of colors for the samples
sample_colors <- setNames(project_labels$color, project_labels$sample)

# Create the ggplot object
ggplot(blitz_covNdet, aes(x = sample, y = reference_genome, size = weighted_mean_cov, color = weighted_detection)) +
  geom_point(alpha = 0.7) +
  scale_color_viridis_c() +
  scale_size_continuous(range = c(1, 10)) +
  labs(title = "Reference Genomes Presence, Mean Coverage, and Detection in Samples",
       x = "Sample",
       y = "Reference Genome") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5, color = sample_colors))
```

[![blitz](images/blitzOutputR.pdf)](images/){:.center-img .width-90}


### Getting collection and bin for `--gene-mode`

First, we started the anvio-dev conda environment on our local machine
```bash
conda list env
conda activate anvio-dev
```

Then we looked at the [anvi-script-add-default-collection](https://anvio.org/help/main/programs/anvi-script-add-default-collection/) program, which can help us get a collection and bin.

The basic command is:

```bash
anvi-script-add-default-collection -c [contigs-db](https://anvio.org/help/main/artifacts/contigs-db) \
                                   -p [profile-db](https://anvio.org/help/main/artifacts/profile-db)
```
The program will add a new collection into the profile database named `DEFAULT`, which will contain a single bin that describes all items in the database named `EVERYTHING`. 


For HIMB2204 and HTCC1002, this is what we ran:

Add default collection to HTCC1002
```
anvi-script-add-default-collection -c HIMB2204nHTCC1002_postFilter/HIMB2204/CONTIGS.db -p HIMB2204nHTCC1002_postFilter/HIMB2204/PROFILE.db 
```
 
Add default collection to HTCC1002
```
anvi-script-add-default-collection -c HIMB2204nHTCC1002_postFilter/HTCC1002/CONTIGS.db -p HIMB2204nHTCC1002_postFilter/HTCC1002/PROFILE.db 
```

### Getting genes databases

To generate a gene database, anvi'o offers to create one when `anvi-interactive` is started in `--gene-mode`, or alternatively with the program `anvi-gen-gene-level-stats-databases`. We used the latter.

The basic command is
```bash
anvi-gen-gene-level-stats-databases -c contigs-db \
                                    -p profile-db \
                                    -C collection
```

Adapted to our data and aim, it is

HIMB2204
```
anvi-gen-gene-level-stats-databases -c HIMB2204nHTCC1002_postFilter/HIMB2204/CONTIGS.db \
                                    -p HIMB2204nHTCC1002_postFilter/HIMB2204/PROFILE.db  \
                                    -C DEFAULT
```

HTCC1002
```
anvi-gen-gene-level-stats-databases -c HIMB2204nHTCC1002_postFilter/HTCC1002/CONTIGS.db \
                                    -p HIMB2204nHTCC1002_postFilter/HTCC1002/PROFILE.db  \
                                    -C DEFAULT
```

The genes database will automatically be stored in the directory in which the PROFILE.db and CONTIGS.db are also in, but in its own subdirectory called `GENES/`, which contains the `DEFAULT-EVERYTHING.db` GENES database.

---

add metadata to them (in `P3_visualise`)

```
import pandas as pd

# Load the weighted results file and metadata file
weighted_results = pd.read_csv('../digital_sup/filter_covNdet/weighted_results_with_length.txt', sep='\t')
metadata = pd.read_csv('metadata/metagenomes_filtered_anvio.txt', sep='\t')

# Step 1: Filter the weighted results dataframe for both HIMB2204 and HTCC1002
filtered_results = weighted_results[
    ((weighted_results['reference_genome'] == 'HIMB2204') | 
     (weighted_results['reference_genome'] == 'HTCC1002')) &
    (weighted_results['weighted_mean_cov'] >= 10) &
    (weighted_results['weighted_detection'] >= 0.5)
]

# Step 2: Extract the unique sample names for HIMB2204 and HTCC1002
filtered_samples = filtered_results['sample'].unique()

# Step 3: Filter the metadata dataframe using the filtered sample names
filtered_metadata = metadata[metadata['sample'].isin(filtered_samples)]

# Step 4: Save the filtered metadata to a new file
filtered_metadata.to_csv('HIMB2204_HTCC1002_metagenomes_filtered_anvio.txt', sep='\t', index=False)

```
check
```
wc -l HIMB2204_HTCC1002_metagenomes_filtered_anvio.txt 
     184 HIMB2204_HTCC1002_metagenomes_filtered_anvio.txt
```
good (1 line more than we have samples between them: header line)

### Import the metadata we want to order our layers (metagenomes) by

#### Prepare

To visualise the metagenomes capturing each reference genome with their metadata in mind (here, following the latitudinal gradient), we can use `anvi-import-misc-data`.

[`anvi-import-misc-data`](https://anvio.org/help/main/programs/anvi-import-misc-data/) is a program to populate additional data or order tables in pan or profile databases for items and layers, OR additional data in contigs databases for nucleotides and amino acids. We want to use it to populate additional data in profile databases for layers.

This blog post gives some more information on this program: [https://merenlab.org/2017/12/11/additional-data-tables/#layers-additional-data-table](https://merenlab.org/2017/12/11/additional-data-tables/#layers-additional-data-table).


Okay, so we need a table to give this program. A table that should look something like this:
```bash
samples numerical_01    numerical_02    categorical     stacked_bar!X   stacked_bar!Y   stacked_bar!Z
c1      100     5       A       1       2       3
c2      200     4       B       2       3       1
c3      300     3       B       3       1       2
```

The first column of this file needs to be identical to our layer (sample) names, and every column describes a property of a given layer.

The basic command is structured as follows: 
- `layers_additional_data.txt` is the table we see above, containing any metadata on the layers we want to add
- `-p profile.db` is the PROFILE.db containing all the information for how our metagenomes short reads map to the reference genome contigs (in our case the PROFILE.db made with only the samples that passed the filtering)
- `--target-data-table layers` notes that what we are providing our extra metadata to is the `layers` object


```bash
anvi-import-misc-data layers_additional_data.txt \
                         -p profile.db \
                         --target-data-table layers
```
and if you want to delete it later 
```bash
anvi-delete-misc-data -p profile.db \
                      --target-data-table layers 
                      --just-do-it
```
So now we needed to make sure that our metadata is following the format that this program expects. That means using the sample names we used as part of this workflow, but also to subset the main metadata table we made in [public-marine-omics-metadata](https://github.com/merenlab/public-marine-omics-metadata/tree/main) to only keep info on the samples that passed our filtering based on coverage and detection here.

Let's do it!

Our metadata file (metagenomes.txt) still had metadata for more samples than the ones that passed the cov and det filtering. 

We first subsetted `metagenomes.txt` to **only include those samples that passed the filtering** and were thus listed in the `unique_samples.txt` we created earlier.

Further, the `metagenomes.txt` file currently has multiple rows per sample if there are multiple runs associated with one sample. We will make it such that there is **only one row per sample** (so one row per layer that we want to associate this with in anvi'o). Of course, ensuring that any differing values across different runs from the same sample are concatenated, and identical values retained as they are.

```bash
nano filterMetagenomesTxt.py
```
content
```python
import pandas as pd

# Step 1: Read the unique samples file
with open('../../../digital_sup/metadataForAnvio/unique_samples.txt', 'r') as f:
    unique_samples = f.read().splitlines()

# Step 2: Read the metagenomes file
metagenomes_df = pd.read_csv('../data/metagenomes.txt', sep='\t')

# Step 3: Create a dictionary to map `biosample` to the full sample name
sample_mapping = {sample.split('_')[1]: sample for sample in unique_samples}

# Step 4: Filter the metagenomes data where `biosample` matches the key in `sample_mapping`
filtered_df = metagenomes_df[metagenomes_df['biosample'].isin(sample_mapping.keys())].copy()

# Step 5: Add a new column `sample` with the full sample name
filtered_df['sample'] = filtered_df['biosample'].map(sample_mapping)

# Step 6: Group by the 'sample' column and aggregate each column by concatenating unique values
aggregated_df = filtered_df.groupby('sample').agg(lambda x: '|'.join(sorted(set(x.dropna().astype(str)))))

# Step 7: Reset the index to make 'sample' a column again
aggregated_df.reset_index(inplace=True)

# Step 8: Reorder columns to make `sample` the first column
cols = ['sample'] + [col for col in aggregated_df.columns if col != 'sample']
aggregated_df = aggregated_df[cols]

# Step 9: Write the resulting DataFrame to a new file
aggregated_df.to_csv('../data/metagenomes_filtered.txt', sep='\t', index=False)
```
run
```
python3 filterMetagenomesTxt.py
```

We then **selected only the columns we wanted to bring into anvi'o**.

```bash
nano prepAnvioMetadata.py
```
content
```python
import pandas as pd

# Load the dataset
file_path = '../data/metagenomes_filtered.txt'  # replace with your file path
df = pd.read_csv(file_path, sep='\t')

# Select only the desired columns and create a copy
columns_to_keep = ['sample', 'run', 'latitude', 'longitude', 'bases', 'spots', 'distance_from_equator', 'project', 'depth', 'temperature_degC', 'year', 'model', 'layer', 'season', 'sub_region_seavox', 'region_seavox', 'provdescr_longhurst', 'environment', 'env_biome', 'env_feature', 'env_material', 'pressure_dbar', 'salinity_pss', 'oxygen_umolKg', 'phosphate_umolKg', 'silicate_umolKg', 'nitrate_umolKg', 'nitrite_umolKg', 'dic_umolKg', 'doc_umolKg', 'chlorophyll', 'chla_ngL', 'chlb_ngL', 'chlc_ngL']
filtered_df = df[columns_to_keep].copy()

# Convert the 'year' column to float
filtered_df['year'] = filtered_df['year'].astype(float)

# Save the filtered DataFrame to a new file
filtered_df.to_csv('../data/metagenomes_filtered_anvio.txt', sep='\t', index=False)
```

The output file `metagenomes_filtered_anvio.txt` still has metadata for all samples that somehow passed the filtering, but - for now - we **only wanted the metadata associated with the samples that are relevant to HIMB2204 and HTCC1002**. For that, we wrote another script.

```python
import pandas as pd

# Load the weighted results file and metadata file
weighted_results = pd.read_csv('../digital_sup/filter_covNdet/weighted_results_with_length.txt', sep='\t')
metadata = pd.read_csv('metadata/metagenomes_filtered_anvio.txt', sep='\t')

# Step 1: Filter the weighted results dataframe for both HIMB2204 and HTCC1002
filtered_results = weighted_results[
    ((weighted_results['reference_genome'] == 'HIMB2204') | 
     (weighted_results['reference_genome'] == 'HTCC1002')) &
    (weighted_results['weighted_mean_cov'] >= 10) &
    (weighted_results['weighted_detection'] >= 0.5)
]

# Step 2: Extract the unique sample names for HIMB2204 and HTCC1002
filtered_samples = filtered_results['sample'].unique()

# Step 3: Filter the metadata dataframe using the filtered sample names
filtered_metadata = metadata[metadata['sample'].isin(filtered_samples)]

# Step 4: Save the filtered metadata to a new file
filtered_metadata.to_csv('HIMB2204_HTCC1002_metagenomes_filtered_anvio.txt', sep='\t', index=False)
```

In addition to the metadata already included in the metadata file, we wanted to assign the samples to the different **climate zones**. To do so, we ran:
```python
import pandas as pd

# Step 1: Load the data
file_path = 'HIMB2204_HTCC1002_metagenomes_filtered_anvio.txt'
data = pd.read_csv(file_path, sep='\t')  # Assuming tab-separated file

# Step 2: Function to classify marine climate zone based on latitude
def classify_marine_zone(lat, lon):
    if abs(lat) <= 10:
        return "A Equatorial"
    elif 10 < abs(lat) <= 23.5:
        return "B Tropical"
    elif 23.5 < abs(lat) <= 40:
        return "C Subtropical"
    elif 40 < abs(lat) <= 60:
        return "D Temperate"
    else:
        return "E Polar"

# Step 3: Apply the function to each row
data['Marine_Climate_Zone'] = data.apply(lambda row: classify_marine_zone(row['latitude'], row['longitude']), axis=1)

# Step 4: Save the updated data to a new file
output_file = 'HIMB2204_HTCC1002_metagenomes_with_marine_climate_zones.txt'
data.to_csv(output_file, sep='\t', index=False)

print(f"File saved as {output_file}")
```


#### Do it

Since we are using `--gene-mode`, we need to import the metadata into the GENES database.

This is the metadata file `metadata/HIMB2204_HTCC1002_metagenomes_filtered_anvio.txt` we prepared above.

add metadata to genes db HIMB2204
```
anvi-import-misc-data ../HIMB2204_HTCC1002_metagenomes_with_marine_climate_zones.txt -p HIMB2204nHTCC1002_postFilter/HIMB2204/GENES/DEFAULT-EVERYTHING.db --target-data-table layers
```

add metadata to genes db HTCC1002
```
anvi-import-misc-data ../HIMB2204_HTCC1002_metagenomes_with_marine_climate_zones.txt -p HIMB2204nHTCC1002_postFilter/HTCC1002/GENES/DEFAULT-EVERYTHING.db --target-data-table layers
```

Since we wanted to also show it on the top of the phylogram, we also needed to import it into the PROFILE.db.

add metadata to profile db HIMB2204
```
anvi-import-misc-data ../HIMB2204_HTCC1002_metagenomes_with_marine_climate_zones.txt -p HIMB2204nHTCC1002_postFilter/HIMB2204/PROFILE.db --target-data-table layers
```

add metadata to genes db HTCC1002
```
anvi-import-misc-data ../HIMB2204_HTCC1002_metagenomes_with_marine_climate_zones.txt -p HIMB2204nHTCC1002_postFilter/HTCC1002/PROFILE.db --target-data-table layers
```

### Look at it

To look at it, we finally used the long-awaited `anvi-interactive` command. With all the prep we did above, it was all set to show the genes in the metagenomes and allow us to sort layers (metagenomes) based on the metadata keys we imported as well as adjust what will be visualized (coverage, detection, ...) how the genes will be ordered in the metagenomes (synteny, detection, ...) and so on.


Looked at it: HIMB2204
```
anvi-interactive -c HIMB2204nHTCC1002_postFilter/HIMB2204/CONTIGS.db -p HIMB2204nHTCC1002_postFilter/HIMB2204/PROFILE.db -C DEFAULT -b EVERYTHING --gene-mode
```

Anvi'o allows us to define a state and export and import that state. 

To do so, we will adjust any parameter we want to adjust and then click on `Save` under the mention of `State` in the left-hand panel in the interactive interface. That prompts us to select a name: `default_HIMB2204nHTCC1002`. This will save the state in the GENES.db

We can also use the same state with other projects. To do so, we can use `anvi-export-state` on the GENES.db in which we defined the state and `anvi-import-state` on the GENES.db that we want to import the state to. 

We exported the state defined in HIMB2204 and then exported it with:
```
anvi-export-state -s default_HIMB2204nHTCC1002_climate -p HIMB2204nHTCC1002_postFilter/HIMB2204/GENES/DEFAULT-EVERYTHING.db -o default_HIMB2204nHTCC1002_climate.json
```

Then imported it into HTCC1002
```
anvi-import-state -p HIMB2204nHTCC1002_postFilter/HTCC1002/GENES/DEFAULT-EVERYTHING.db -s default_HIMB2204nHTCC1002_climate.json -n default_HIMB2204nHTCC1002
```

Looked at it: HTCC1002
```
anvi-interactive -c HIMB2204nHTCC1002_postFilter/HTCC1002/CONTIGS.db -p HIMB2204nHTCC1002_postFilter/HTCC1002/PROFILE.db -C DEFAULT -b EVERYTHING --gene-mode
```

{:.notice}
The state `default_HIMB2204nHTCC1002_climate` was used to visualize both phylograms in Panel C. You can get the state file [here](files/default_HIMB2204nHTCC1002_climate.json).

{:.notice}
A different state, called `default_HIMB2204nHTCC1002_geneDet` was used to visualize the phylogram zoom-in in Panel D. You can get the state file [here](files/default_HIMB2204nHTCC1002_geneDet.json).



## Visualising world maps

For the world maps seen in Panels A and B of Figure 1, we mainly relied on R.

### Importing all the data we need into R

First, we imported all relevant data into R and made sure the columns were understood correctly. This resulted in the script becoming quite long, so we will hide it here, but feel free to expand it if you are curious.

The files we imported with were
- metagenomes_filtered_anvio.txt (metadata)
- weighted_results_with_length.txt (coverage and detection data)

We then combined both into one dataframe.

<details markdown="1"><summary>Click to see how we imported the data into R</summary>

```
# import metadata used for anvio

# 1. metadata
# 2. Coverage and detection data

####################################
## prep
####################################

# Set working directory
setwd("/Users/rameyer/Documents/_P3/P3dataAnalysis/P3_metadata/public-marine-omics-metadata/")

# Load necessary libraries
library(readr)
library(dplyr)

####################################
## go 1. metadata
####################################


# import metadata table for metagenomes that passed anvi'o filtering of
# having at least one reference genome with 0.5 detection and 10x coverage

# make sure all are understood properly

# Read in the data
metagenomes_filtered_anvio <- read_delim(
  "data/metagenomes_filtered_anvio.txt", 
  delim = "\t", 
  escape_double = FALSE, 
  col_types = cols(
    latitude = col_number(), 
    longitude = col_number(), 
    distance_from_equator = col_number(), 
    depth = col_number(), 
    temperature_degC = col_number(), 
    year = col_number(),
    spots = col_character()
  ), 
  trim_ws = TRUE
)

# Use mutate() to adjust the factors
metagenomes_filtered_anvio <- metagenomes_filtered_anvio %>%
  mutate(
    project = factor(project, levels = c("BATS", "HOT3", "HOT1", "MAL", "BGS", "BGT", "TARA", "OSD")),
    model = factor(model, levels = c("NextSeq 500", "NextSeq 550", "Illumina HiSeq 2000", "Illumina HiSeq 4000", "Illumina NovaSeq 6000", "Illumina MiSeq", "Illumina Genome Analyzer IIx")),
    layer = factor(layer, levels = c("surface water", "deep chlorophyll maximum")),
    season = factor(season, levels = c("Spring", "Fall", "Summer", "Winter")),
    sub_region_seavox = factor(sub_region_seavox, levels = c(
      "NORTHWEST ATLANTIC OCEAN (40W)", "LABRADOR SEA", "CELTIC SEA", 
      "NORTHEAST ATLANTIC OCEAN (40W)", "SOUTHWEST ATLANTIC OCEAN (20W)", 
      "INDIAN OCEAN", "ARABIAN SEA", "BAY OF BENGAL", 
      "NORTHEAST PACIFIC OCEAN (180W)", "SOUTHEAST PACIFIC OCEAN (140W)", 
      "SOUTHERN OCEAN", "SOUTHEAST ATLANTIC OCEAN (20W)", 
      "SOUTHWEST PACIFIC OCEAN (140W)", "CORAL SEA", "TASMAN SEA", 
      "BASS STRAIT", "ICELAND SEA", "LILLEBAELT", "ADRIATIC SEA", 
      "MEDITERRANEAN SEA, EASTERN BASIN", "STRAIT OF SICILY", "IONIAN SEA", 
      "RED SEA", "LAKSHADWEEP SEA", "MOZAMBIQUE CHANNEL", "DRAKE PASSAGE", 
      "WEDDELL SEA", "CARIBBEAN SEA", "GULF OF MEXICO")),
    region_seavox = factor(region_seavox, levels = c(
      "ATLANTIC OCEAN", "INDIAN OCEAN", "PACIFIC OCEAN", 
      "SOUTHERN OCEAN", "ARCTIC OCEAN", "BALTIC SEA", 
      "MEDITERRANEAN REGION")),
    provdescr_longhurst = factor(provdescr_longhurst, levels = c(
      "Westerlies - N. Atlantic Subtropical Gyral Province (West) (STGW)",
      "Polar - Atlantic Arctic Province", 
      "Westerlies - N. Atlantic Drift Province (WWDR)",
      "Coastal - NW Atlantic Shelves Province", 
      "Westerlies - Gulf Stream Province", 
      "Coastal - NE Atlantic Shelves Province",
      "Trades - N. Atlantic Tropical Gyral Province (TRPG)", 
      "Trades - Western Tropical Atlantic Province",
      "Trades - South Atlantic Gyral Province (SATG)", 
      "Westerlies - S. Subtropical Convergence Province",
      "Westerlies - Subantarctic Province", 
      "Trades - Indian S. Subtropical Gyre Province",
      "Trades - Indian Monsoon Gyres Province", 
      "Coastal - NW Arabian Upwelling Province", 
      "Trades - N. Pacific Tropical Gyre Province",
      "Trades - N. Pacific Equatorial Countercurrent Province", 
      "Trades - Pacific Equatorial Divergence Province",
      "Westerlies - S. Pacific Subtropical Gyre Province", 
      "Polar - Antarctic Province", 
      "Polar - Austral Polar Province",
      "Coastal - SW Atlantic Shelves Province", 
      "Trades - Archipelagic Deep Basins Province", 
      "Coastal - East Australian Coastal Province",
      "Coastal - Benguela Current Coastal Province", 
      "Coastal - E. Africa Coastal Province", 
      "Westerlies - N. Atlantic Subtropical Gyral Province (East) (STGE)", 
      "Westerlies - Mediterranean Sea, Black Sea Province", 
      "Coastal - Canary Coastal Province (EACB)", 
      "Coastal - Red Sea, Persian Gulf Province", 
      "Coastal - Chile-Peru Current Coastal Province", 
      "Coastal - California Upwelling Coastal Province", 
      "Coastal - Central American Coastal Province", 
      "Coastal - Guianas Coastal Province", 
      "Trades - Caribbean Province")),
    environment = factor(environment, levels = c("Westerlies", "Polar", "Coastal", "Trades"))
  )

#View(metagenomes_filtered_anvio)
head(metagenomes_filtered_anvio)



####################################
## go 2. coverage and detection data
####################################

######
# https://github.com/merenlab/world-map-r/blob/master/generate-PDFs.R expects 
# relative abundance. I don't have that, so will use coverage or detection 
# instead. Let's get those into R (from the blitz outputs)
######

# Step 1: Read the file into R
weighted_results <- read.delim("/Users/rameyer/Documents/_P3/P3dataAnalysis/digital_sup/filter_covNdet/weighted_results_with_length.txt", sep="\t")

# Step 2: Reshape the data for weighted_mean_cov
# We want to create one column for each reference_genome and fill it with weighted_mean_cov values
weighted_cov_reshaped <- weighted_results %>%
  select(sample, reference_genome, weighted_mean_cov) %>%  # Select relevant columns
  pivot_wider(names_from = reference_genome, values_from = weighted_mean_cov, names_prefix = "cov_")  # Pivot to wide format

# Step 3: Reshape the data for weighted_detection
weighted_detection_reshaped <- weighted_results %>%
  select(sample, reference_genome, weighted_detection) %>%  # Select relevant columns
  pivot_wider(names_from = reference_genome, values_from = weighted_detection, names_prefix = "det_")  # Pivot to wide format

# Step 4: Merge the reshaped data with metagenomes_filtered_anvio
# First, merge the cov_ columns
metagenomes_filtered_anvio_cov <- metagenomes_filtered_anvio %>%
  left_join(weighted_cov_reshaped, by = "sample")

# Then, merge the det_ columns
metagenomes_filtered_anvio_covNdev <- metagenomes_filtered_anvio_cov %>%
  left_join(weighted_detection_reshaped, by = "sample")

# Step 5: Check the result
head(metagenomes_filtered_anvio_covNdev)
colnames(metagenomes_filtered_anvio_covNdev)

View(metagenomes_filtered_anvio_covNdev)
```

</details>


### Normalising coverage data
Then we went on to normalize the coverage data. The different projects had very different sequencing depth. This normalization was to avoid the inter-project differences visually overshadowing the ecological differences.

The normalization simply was to divide the coverage by the number of reads in the sample (given as `spots` in the metadata). Projects with deeper sequencing had a much higher number of reads. In cases where multiple runs from the same sample were available, we first averaged the number of reads across runs. To make sure the resulting coverage numbers were not too tiny, we divided the average number of reads per sample by 1 million before using it to normalize the coverage data.

```
# normalize coverage data
# Load the dplyr package for data manipulation

library(dplyr)

normalise_df <- metagenomes_filtered_anvio_covNdev

# 1. Calculate `average_spots` if not already done
calculate_average_spots <- function(spots) {
  spot_values <- as.numeric(unlist(strsplit(spots, "\\|")))
  return(mean(spot_values))
}

normalise_df$average_spots <- sapply(normalise_df$spots, calculate_average_spots)

# Between 1 and 2: Sanity check
### Get range by project 
# Group by project and calculate the range of average_spots
range_by_project <- normalise_df %>%
  group_by(project) %>%
  summarize(min_spots = min(average_spots, na.rm = TRUE), 
            max_spots = max(average_spots, na.rm = TRUE))

# View the range summary
print(range_by_project)

# 2. Scale `average_spots` by dividing by 1 million
normalise_df$scaled_average_spots <- normalise_df$average_spots / 1e6

# 3. Normalize the `cov_[REFERENCE_GENOME]` columns using the scaled `average_spots`
cov_columns <- grep("^cov_", colnames(normalise_df), value = TRUE)

# Create new normalized columns `norm_cov_[REFERENCE_GENOME]`
for (cov_column in cov_columns) {
  new_col_name <- paste0("norm_", cov_column)  # Create new column name
  normalise_df[[new_col_name]] <- normalise_df[[cov_column]] / normalise_df$scaled_average_spots
}

# View the updated dataframe
head(normalise_df)
colnames(normalise_df)

### Get range by project 

# Group by project and calculate the range of norm_cov_HTCC7211
range_by_project <- normalise_df %>%
  group_by(project) %>%
  summarize(min_spots = min(cov_HTCC7211, na.rm = TRUE), 
            max_spots = max(cov_HTCC7211, na.rm = TRUE))

print(range_by_project)

# Group by project and calculate the range of norm_cov_HTCC7211
range_by_project <- normalise_df %>%
  group_by(project) %>%
  summarize(min_spots = min(norm_cov_HTCC7211, na.rm = TRUE), 
            max_spots = max(norm_cov_HTCC7211, na.rm = TRUE))

print(range_by_project)
```

### Creating the global world maps

We then went on to create the world maps for HTCC1002 and HIMB2204.

What is important to note here is that we had any sample that passed the coverage and detection filtering for any reference genome in the background, while the samples from either HTCC1002 or HIMB2204 were highlighted. We also scaled the size of points between the two.
```
# make world maps for HTCC1002 and HIMB2204 where the size of the points 
# (representing normalized coverage) is scaled the same between the two genomes

# Set working directory (customize this as needed)
setwd("/Users/rameyer/Documents/_P3/P3dataAnalysis/P3_visualise/")

##############################
# User Input Section
##############################

# Specify the cutoff values for coverage and detection
coverage_cutoff <- 10
detection_cutoff <- 0.5

# Specify whether to use the highlight column (TRUE or FALSE)
use_highlight <- TRUE  # Set to FALSE if you don't want to use highlight functionality

# Option to save the plots as PDFs
save_pdf <- TRUE  # Set to FALSE if you don't want to save the plots as PDFs

# Specify the PDF file dimensions
PDF_WIDTH <- 12
PDF_HEIGHT <- 5.5

dataframe <- normalise_df

##############################
# Data Preparation: Focus on HTCC1002 and HIMB2204
##############################

# Load necessary libraries
library(ggplot2)
library(scales)
library(maps)

# Define the specific reference genomes to loop over
reference_genomes <- c("HTCC1002", "HIMB2204")

# Calculate global min and max normalized coverage across both genomes for consistent scaling
all_norm_coverage <- c()

for (reference_genome in reference_genomes) {
  norm_coverage_col <- paste0("norm_cov_", reference_genome)
  all_norm_coverage <- c(all_norm_coverage, dataframe[[norm_coverage_col]])
}

global_min_norm_cov <- min(all_norm_coverage, na.rm = TRUE)
global_max_norm_cov <- max(all_norm_coverage, na.rm = TRUE)

# Loop through each reference genome
for (reference_genome in reference_genomes) {
  
  # Define the name of the coverage, detection, and normalized coverage columns for the current reference genome
  coverage_col <- paste0("cov_", reference_genome)
  detection_col <- paste0("det_", reference_genome)
  norm_coverage_col <- paste0("norm_cov_", reference_genome)  # New column for normalized coverage
  
  # Create a dynamic name for the new dataframe
  new_df_name <- paste0("metagenomes_", reference_genome)
  
  # Assign the metagenomes_filtered_anvio_covNdev dataframe to the new dynamically named dataframe
  assign(new_df_name, dataframe)
  
  # Now the dynamically named dataframe is accessible via its name
  df <- get(new_df_name)
  
  ############################################
  # Optional: Add Highlight Column if `use_highlight` is TRUE
  if (use_highlight) {
    # Highlight samples only if both coverage and detection values are above or equal to the cutoffs
    df$highlight <- ifelse(
      df[[coverage_col]] >= coverage_cutoff & df[[detection_col]] >= detection_cutoff, TRUE, FALSE
    )
  }
  
  ############################################
  # Create the World Map Plot
  ############################################
  
  # Basic world map data
  world_map <- map_data("world")
  
  # Define custom color palettes (light, dark, and extralight) for different projects
  custom_palette_extralight <- c(
    "BATS" = "#D8EBF8", "BGS" = "#F6C2D5", "BGT" = "#D1ECEC", "HOT1" = "#E3D2E8",
    "HOT3" = "#fdd2c1", "MAL" = "#ffeacd", "OSD" = "#FBEEC2", "TARA" = "#DCF1E0"
  )
  
  custom_palette_light <- c(
    "BATS" = "#B0D9F2", "BGS" = "#F19CBB", "BGT" = "#A8D8D8", "HOT1" = "#C9B2D6",
    "HOT3" = "#f9a98e", "MAL" = "#ffd8a1", "OSD" = "#F9E29C", "TARA" = "#B9E3C6"
  )
  
  custom_palette_dark <- c(
    "BATS" = "#2B5D9C", "BGS" = "#B13A8C", "BGT" = "#40B2B2", "HOT1" = "#875D9B",
    "HOT3" = "#f15c40", "MAL" = "#fdb64e", "OSD" = "#D8A600", "TARA" = "#469F77"
  )
  
  # Ensure 'project' is a factor and the palettes have colors for each project
  df$project <- factor(df$project, levels = names(custom_palette_light))
  
  # Assign colors based on whether we use highlight or not
  if (use_highlight) {
    # If using highlight, use different colors for highlighted and non-highlighted samples
    df$outline_color <- ifelse(
      df$highlight == TRUE,
      custom_palette_dark[df$project],  # Darker outline for highlighted samples
      custom_palette_extralight[df$project]  # Lighter outline for others
    )
    
    # Assign fill colors based on whether the sample is highlighted
    df$fill_color <- ifelse(
      df$highlight == TRUE,
      custom_palette_light[df$project],  # Light fill for highlighted samples
      custom_palette_extralight[df$project]  # Extralight fill for non-highlighted samples
    )
  } else {
    # If not using highlight, use the light palette for both fill and outline
    df$fill_color <- custom_palette_light[df$project]
    df$outline_color <- custom_palette_light[df$project]
  }
  
  # Split data into highlighted and non-highlighted groups if highlight is used
  if (use_highlight) {
    df_no_highlight <- df[df$highlight == FALSE, ]
    df_highlight <- df[df$highlight == TRUE, ]
  }
  
  # Set outline size for all points
  stroke_size <- 1.5  # Adjust this for thicker outlines
  
  # Create the base plot with the world map
  p <- ggplot() +
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
                 fill = "gray90", color = "white") +
    coord_fixed(1.3)  # Maintains aspect ratio
  
  # Add non-highlighted points if `use_highlight` is TRUE, otherwise add all points
  if (use_highlight) {
    # First layer: plot non-highlighted samples with extralight fill colors
    p <- p + geom_point(data = df_no_highlight,
                        aes(x = longitude, y = latitude, 
                            fill = fill_color,   # Use extralight fill for non-highlighted samples
                            color = outline_color, 
                            size = !!sym(norm_coverage_col)),
                        shape = 21, stroke = stroke_size, alpha = 0.8)
    
    # Second layer: plot highlighted samples with light fill colors (these will appear on top)
    p <- p + geom_point(data = df_highlight,
                        aes(x = longitude, y = latitude, 
                            fill = fill_color,   # Use light fill for highlighted samples
                            color = outline_color, 
                            size = !!sym(norm_coverage_col)),
                        shape = 21, stroke = stroke_size, alpha = 0.8)
  } else {
    # Plot all samples in a single layer if no highlight is used
    p <- p + geom_point(data = df,
                        aes(x = longitude, y = latitude, 
                            fill = fill_color,  # Default to light palette if not using highlight
                            color = outline_color, 
                            size = !!sym(norm_coverage_col)),
                        shape = 21, stroke = stroke_size, alpha = 0.8)
  }
  
  # Add color scale and other aesthetics
  p <- p + scale_fill_identity() +  # Use the fill colors from the fill_color column
    scale_color_identity() +         # Use the color from the outline_color column
    scale_size_continuous(range = c(4, 10), limits = c(global_min_norm_cov, global_max_norm_cov), guide = guide_legend(title = "Normalized Coverage")) +  # Global normalized coverage for size scaling
    theme_minimal() +
    labs(
      title = paste("Project Locations with", reference_genome, "Normalized Coverage"),
      x = "Longitude", y = "Latitude",
      fill = "Project", color = "Highlight", size = "Normalized Coverage per 1 mio reads"
    ) +
    theme(
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.position = "right"
    )
  
  # Show plot
  print(p)
  
  ############################################
  # Optionally Save the Plot as PDF
  ############################################
  
  if (save_pdf) {
    # Create the filename dynamically for each reference genome
    pdf_filename <- paste0("maps_normalised_lighter_spec/Project_Locations_with_", reference_genome, "_Normalized_Coverage.pdf")
    
    # Save the plot as a PDF
    ggsave(filename = pdf_filename, plot = p, width = PDF_WIDTH, height = PDF_HEIGHT)
    cat(paste0("Plot saved as PDF: ", pdf_filename, "\n"))
  }
}

```

### Creating regional world maps

For the regional biogeographical patterns, we focused on HTCC1002 and HTCC7221. To create the four different plots, the selection of `reference genome` (either `reference_genome <- "HTCC1002"` or `reference_genome <- "HTCC7211"`) and the zoom limits need to be manually adjusted (all options given in the script, just commented out).


```
# Zoom in and map for specific reference genome (e.g., HTCC7211)


# Load necessary libraries
library(ggplot2)
library(scales)
library(maps)

# Set working directory (customize as needed)
setwd("/Users/rameyer/Documents/_P3/P3dataAnalysis/P3_visualise/")

# Define the reference genome for which you want to zoom in
reference_genome <- "HTCC1002"

# Define the specific reference genomes you want to consider for the point sizing
reference_genomes <- c("HTCC1002", "HTCC7211")

# Calculate global min and max normalized coverage across both genomes for consistent scaling
all_norm_coverage <- c()

for (reference_genome in reference_genomes) {
  norm_coverage_col <- paste0("norm_cov_", reference_genome)
  all_norm_coverage <- c(all_norm_coverage, dataframe[[norm_coverage_col]])
}

global_min_norm_cov <- min(all_norm_coverage, na.rm = TRUE)
global_max_norm_cov <- max(all_norm_coverage, na.rm = TRUE)

# Specify the cutoff values for coverage and detection
coverage_cutoff <- 10
detection_cutoff <- 0.5

# Specify whether to use the highlight column (TRUE or FALSE)
use_highlight <- TRUE  # Set to FALSE if you don't want to use highlight functionality

# Option to save the zoomed plots as PDFs
save_pdf <- TRUE  # Set to FALSE if you don't want to save the zoomed plots as PDFs

# Specify the PDF file dimensions
PDF_WIDTH <- 12
PDF_HEIGHT <- 5.5

#Specify zoom area for the map (longitude and latitude limits)
zoom_limits <- list(
  xlim = c(-70, -20),  # Set longitude limits for NORTH ATLANTIC
  ylim = c(25, 70)      # Set latitude limits for NORTH ATLANTIC
)

# # Specify zoom area for the map (longitude and latitude limits)
# zoom_limits <- list(
#   xlim = c(-104, -40),  # Set longitude limits for south america
#   ylim = c(-74, -28)      # Set latitude limits for south america
# )

# # Specify zoom area for the map (longitude and latitude limits)
# zoom_limits <- list(
#   xlim = c(-104, 15),  # Set longitude limits for south america and south africa
#   ylim = c(-73, -29)      # Set latitude limits for south america and south africa 
# )

# Extract the relevant columns for the chosen reference genome
coverage_col <- paste0("cov_", reference_genome)
detection_col <- paste0("det_", reference_genome)
norm_coverage_col <- paste0("norm_cov_", reference_genome)

# Assign the metagenomes_filtered_anvio_covNdev dataframe to the new dynamically named dataframe
df <- dataframe  # Assuming dataframe is already loaded

# Optional: Add Highlight Column if `use_highlight` is TRUE
if (use_highlight) {
  df$highlight <- ifelse(
    df[[coverage_col]] >= coverage_cutoff & df[[detection_col]] >= detection_cutoff, TRUE, FALSE
  )
}

# Basic world map data
world_map <- map_data("world")

# Define custom color palettes (light, dark, and extralight) for different projects
custom_palette_extralight <- c(
  "BATS" = "#D8EBF8", "BGS" = "#ffe3edff", "BGT" = "#D1ECEC", "HOT1" = "#E3D2E8",
  "HOT3" = "#fdd2c1", "MAL" = "#ffeacd", "OSD" = "#FBEEC2", "TARA" = "#DCF1E0"
)

custom_palette_light <- c(
  "BATS" = "#B0D9F2", "BGS" = "#F19CBB", "BGT" = "#A8D8D8", "HOT1" = "#C9B2D6",
  "HOT3" = "#f9a98e", "MAL" = "#ffd8a1", "OSD" = "#F9E29C", "TARA" = "#B9E3C6"
)

custom_palette_dark <- c(
  "BATS" = "#2B5D9C", "BGS" = "#B13A8C", "BGT" = "#40B2B2", "HOT1" = "#875D9B",
  "HOT3" = "#f15c40", "MAL" = "#fdb64e", "OSD" = "#D8A600", "TARA" = "#469F77"
)

# Ensure 'project' is a factor and the palettes have colors for each project
df$project <- factor(df$project, levels = names(custom_palette_light))

# Assign colors based on whether we use highlight or not
if (use_highlight) {
  df$outline_color <- ifelse(
    df$highlight == TRUE,
    custom_palette_dark[df$project],  # Darker outline for highlighted samples
    custom_palette_extralight[df$project]  # Lighter outline for non-highlighted samples
  )
  df$fill_color <- ifelse(
    df$highlight == TRUE,
    custom_palette_light[df$project],  # Light fill for highlighted samples
    custom_palette_extralight[df$project]  # Extralight fill for non-highlighted samples
  )
} else {
  # If not using highlight, use the light palette for both fill and outline
  df$fill_color <- custom_palette_light[df$project]
  df$outline_color <- custom_palette_light[df$project]
}

# # Determine the minimum and maximum normalized coverage values for sizing
# min_norm_cov <- min(df[[norm_coverage_col]], na.rm = TRUE)
# max_norm_cov <- max(df[[norm_coverage_col]], na.rm = TRUE)

# Split data into highlighted and non-highlighted groups if highlight is used
if (use_highlight) {
  df_no_highlight <- df[df$highlight == FALSE, ]
  df_highlight <- df[df$highlight == TRUE, ]
}

# Set outline size for all points
stroke_size <- 1.5  # Adjust this for thicker outlines

# Create the base plot with the world map
p <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = "gray90", color = "white") +
  coord_quickmap(xlim = zoom_limits$xlim, ylim = zoom_limits$ylim)  # Zoom into the specified area and maintain aspect ratio

# Add non-highlighted points if `use_highlight` is TRUE, otherwise add all points
if (use_highlight) {
  p <- p + geom_point(data = df_no_highlight,
                      aes(x = longitude, y = latitude, fill = fill_color, color = outline_color, size = !!sym(norm_coverage_col)),
                      shape = 21, stroke = stroke_size, alpha = 0.8)
  
  p <- p + geom_point(data = df_highlight,
                      aes(x = longitude, y = latitude, fill = fill_color, color = outline_color, size = !!sym(norm_coverage_col)),
                      shape = 21, stroke = stroke_size, alpha = 0.8)
} else {
  p <- p + geom_point(data = df,
                      aes(x = longitude, y = latitude, fill = fill_color, color = outline_color, size = !!sym(norm_coverage_col)),
                      shape = 21, stroke = stroke_size, alpha = 0.8)
}

# Add color scale and other aesthetics
p <- p + scale_fill_identity() +  # Use the fill colors from the fill_color column
  scale_color_identity() +         # Use the outline colors from the outline_color column
  scale_size_continuous(range = c(6, 12), limits = c(global_min_norm_cov, global_max_norm_cov), guide = guide_legend(title = "Normalized Coverage")) +
  theme_minimal() +
  labs(
    title = paste("Zoomed Project Locations with", reference_genome, "Normalized Coverage"),
    x = "Longitude", y = "Latitude",
    fill = "Project", color = "Highlight", size = "Normalized Coverage"
  ) +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.position = "right"
  )

# Show plot
print(p)

# Optionally Save the Zoomed Plot as PDF
if (save_pdf) {
  # Create the filename dynamically for the zoomed map
  pdf_filename <- paste0("maps_normalised_lighter_zoom/Zoomed_Project_Locations_with_", reference_genome, "_Normalized_Coverage_north.pdf")
  
  # Save the plot as a PDF
  ggsave(filename = pdf_filename, plot = p, width = PDF_WIDTH, height = PDF_HEIGHT)
  cat(paste0("Zoomed plot saved as PDF: ", pdf_filename, "\n"))
}
```

## Checking the ANI between the visualized reference genomes

To understand how similar or dissimilar the reference genomes we show in Figure 1 are from one another, we - as the last step - checked their ANI using FastANI (https://github.com/ParBLiSS/FastANI) version 1.33.

To do so, we specified the location of the fasta files of each reference genomes, as well as the output file. 

First for `HTCC1002` and `HIMB2204`
```
fastANI -q fastaDrep/dereplicated_genomes/HTCC1002.fa -r fastaDrep/dereplicated_genomes/HIMB2204.fa --visualize -o output_fastANI_HTCC1002nhIMB2204.txt
```

**Output:** fastaDrep/dereplicated_genomes/HTCC1002.fa      fastaDrep/dereplicated_genomes/HIMB2204.fa       78.0557 199     441


And second for `HTCC1002` and `HTCC7211`
```bash
fastANI -q fastaDrep/dereplicated_genomes/HTCC1002.fa -r fastaDrep/dereplicated_genomes/HTCC7211.fa --visualize -o output_fastANI_HTCC1002nHTCC7211.txt
```
**Output:** fastaDrep/dereplicated_genomes/HTCC1002.fa      fastaDrep/dereplicated_genomes/HTCC7211.fa       79.4378 241     441

To visualise them, we utilized the R script FastANI provides:

First for `HTCC1002` and `HIMB2204`
```bash
Rscript visualizeFastANI.R ../../P3_refGenomes/FASTA_collection_used/HTCC1002.fa ../../P3_refGenomes/FASTA_collection_used/HIMB2204.fa output_fastANI_HTCC1002nhIMB2204.txt.visual
```

[![dRep]({{images}}/output_fastANI_HTCC1002nhIMB2204.txt.visual.pdf)]({{images}}/output_fastANI_HTCC1002nhIMB2204.txt.visual.pdf){:.center-img .width-70}


And second for `HTCC1002` and `HTCC7211`
```bash
Rscript visualizeFastANI.R ../../P3_refGenomes/FASTA_collection_used/HTCC1002.fa ../../P3_refGenomes/FASTA_collection_used/HTCC7211.fa output_fastANI_HTCC1002nHTCC7211.txt.visual
```

[![dRep]({{images}}/output_fastANI_HTCC1002nHTCC7211.txt.visual.pdf)]({{images}}/output_fastANI_HTCC1002nHTCC7211.txt.visual.pdf){:.center-img .width-70}












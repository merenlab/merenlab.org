---
layout: page
title: A reproducible workflow for Meyer & Eren manuscript 2024
modified: 2024-08-30
excerpt: "A bioinformatics workflow for our metagenomics data integration across observatories"
comments: true
authors: [raissa]
---

<div class="extra-info" markdown="1">

<span class="extra-info-header">Summary</span>

**The purpose of this page** is to provide access to reproducible data products and analyses for the manuscript "**Marine Microbial Observatories for the Future: From Samples to Data to Legacy using integrated omics strategies**" by Meyer & Eren


In this study, we 
-  integrate metagenomics metadata and data from observatories (Hawaii Ocean Time Series and Bermuda Atlantic Time-series Study), sampling expeditions (Bio-GO-SHIP, bioGEOTRACES, Malaspina, and Tara Oceans), and citizen science initiatives (Ocean Sampling Day)
-  generate an anvi'o contigs database describing 51 SAR11 isolate genomes (reference genomes)
-  recruite reads from the above listed projects' metagenomes to the SAR11 reference genomes, and profil the recruitment results
-  investigate the patterns in genes across metagenomes recruited to the individual SAR11 reference genomes

Sections in this document will detail all the steps of downloading and processing SAR11 genomes and metagenomes, mapping metagenomic reads onto the SAR11 genomes, as well as analysing and visualising the outcomes.

For the curation of metadata, please consult the [public-marine-omics-metadata](https://github.com/merenlab/public-marine-omics-metadata/tree/main) GitHub repository, where all steps of metadata gathering and standardisation are described in detail.

{:.notice}
If you have any questions, notice an issue, and/or are unable to find an important piece of information here, please feel free to leave a comment down below, send an e-mail to [us](/people/), or get in touch with us through Discord:

{% include _join-anvio-discord.html %}

Information for the projects included in this analysis (information only for depth up to 100 m)
Project | Acronym | Accession | Years covered | # Samples 
| - | - | - | - | - |
[Bermuda Atlantic Time-series Study](https://bats.bios.asu.edu/about/) | BATS | PRJNA385855 | 2003 - 2004, 2009 | 40 
[Bio-GO-SHIP](https://biogoship.org) | BGS | PRJNA656268 | 2011, 2013, 2014, 2016 - 2018 | 969 
[bioGEOTRACES](https://www.nature.com/articles/sdata2018176) | BGT | PRJNA385854 | 2010, 2011 | 323 
[Hawaii Ocean Time-series](http://hahana.soest.hawaii.edu/hot/hot_jgofs.html) | HOT1 \| HOT3 |  PRJNA385855 \| PRJNA352737 | 2003, 2004 \| 2014 - 2017 | 28 \| 230
[Malaspina](https://www.nature.com/articles/s41597-024-02974-1) | MAL | PRJEB52452 | 2011 | 16
[Ocean Sampling Day 2014](https://doi.org/10.1186/s13742-015-0066-5) | OSD | PRJEB8682 | 2014 | 149
[Tara Oceans](https://fondationtaraocean.org/en/home/) | TARA | PRJEB1787 | 2009 - 2012 | 93


</div>




## Metadata
For the curation of metadata, please consult the [public-marine-omics-metadata](https://github.com/merenlab/public-marine-omics-metadata/tree/main) GitHub repository, where all steps of metadata gathering and standardisation are described in detail.

Remember, team: Data without metadata is like a recipe without ingredients listed - you know what you'd get out but not how you got there. Our data is only worth as much as our metadata can support it.

## SAR11 cultivar genomes
This section explains how to prepare the set of SAR11 isolate genomes, ending up with 51 reference genomes. 

The isolate genomes available at the time of this analysis are included in `SAR11_June2024_bycontig.fa`. To see more information and the source of each isolate genome, expand the section below.

---

<details markdown="1"><summary>Click here to show/hide an overview of the reference genomes used</summary>

<details>
<summary>Click here to see an overview of the reference genomes used</summary>

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
| HIMB123                             | Ia.3.V    | Ia.3.V         | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1420                            | Ia.3.II   | Ia.3.II        | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB1427                            | Ia.3.II   | Ia.3.II        | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB1520                            | Ia.3.II   | Ia.3.II        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1564                            | Ia.3.II   | Ia.3.II        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1412                            | Ia.3.II   | Ia.3.II_1412   | Type | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB1505                            | Ia.3.III  | Ia.3.III       | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1527                            | Ia.3.III  | Ia.3.III       | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1549                            | Ia.3.III  | Ia.3.III       | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1409                            | Ia.3.V    | Ia.3.V         | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB1413                            | Ia.3.V    | Ia.3.V         | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB1444                            | Ia.3.V    | Ia.3.V         | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB1456                            | Ia.3.V    | Ia.3.V         | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB2201                            | Ia.3.V    | Ia.3.V         | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB2250                            | Ia.3.V    | Ia.3.V         | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB1430                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB1485                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1488                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1490                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1491                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1493                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1494                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1495                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1506                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1507                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1509                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1513                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1518                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1521                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1526                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1542                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1552                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1556                            | Ia.3.VI   | Ia.3.VI        | Type | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1559                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1573                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1577                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1587                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1593                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1597                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1611                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1623                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1631                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1636                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1641                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1662                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1685                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1695                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1701                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1702                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1709                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1710                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1715                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1723                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1746                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1748                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1758                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1765                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1770                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1782                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB2211                            | Ia.3.VI   | Ia.3.VI        | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB1483                            | Ia.4.1483 | Ia.4.1483      | Type | Freel et al., in prep        | coastal, Oahu, Hawaii, USA  |
| HIMB1402                            | Ia.4.II   | Ia.4.II        | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB1437                            | Ia.4.II   | Ia.4.II        | Type | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB1863                            | Ia.4.II   | Ia.4.II        | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB2104                            | Ia.4.II   | Ia.4.II        | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB2200                            | Ia.4.II   | Ia.4.II        | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB2204                            | Ia.4.II   | Ia.4.II        | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB2215                            | Ia.4.II   | Ia.4.II        | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB1436                            | Ia.4.RS39 | Ia.4.RS39      | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB2187                            | Ia.4.RS39 | Ia.4.RS39      | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB2226                            | Ia.4.RS39 | Ia.4.RS39      | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB2305                            | Ib.1      | Ib.1           | No   | Freel et al., in prep        | offshore, Oahu, Hawaii, USA |
| HIMB114                             | IIIa      | IIIa           | Type | Freel et al., in prep        |                             |
| HIMB1517                            | IIIa      | IIIa           | No   | Freel et al., in prep        |                             |
| HIMB1565                            | IIIa      | IIIa           | No   | Freel et al., in prep        |                             |
| IMCC9063                            | IIIa      | IIIa           | No   | Oh et al., 2011              |                             |
| LSUCC0530                           | IIIb      | IIIb           | Type | Henson et al., 2018          |                             |
| HIMB58                              | IIa       | IIa            | Type | Brandon M.S. thesis (2006)   |                             |

</details>

--- 

For the next steps, we need to separate `SAR11_June2024_bycontig.fa` into individual .fa files - one per reference genome.

For that, use the following script, running it in the same directory as you have the `SAR11_June2024_bycontig.fa` file.
```
nano separateFasta.py
```
content
```
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
            if line.startswith('>'):  # Check if the line is a header
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
```
python3 separateFasta.py
```

You should now hav 99 individual fasta files, ready to be evaluated.

### Use checkM to evaluate quality of isolate genomes

We will use CheckM to evaluate the completeness and contamination of the isolate genomes. Since they are isolate genomes, we expect them to be of quite high completeness and low contamination.

> [!NOTE]
There is also an option to include CheckM in the dRep step that is following, however, that does not seem to work for everyone, so we will do it separately.

Citation
> CheckM citation:
> 
> - Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW. 2015. [CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes](http://genome.cshlp.org/content/25/7/1043.short). Genome Research, 25: 1043–1055.
> 
> CheckM relies on several other software packages:
> 
> - [pplacer](http://matsen.fhcrc.org/pplacer/): Matsen FA, Kodner RB, Armbrust EV. 2010. pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree. BMC Bioinformatics 11: doi:10.1186/1471-2105-11-538.
> - [prodigal](http://prodigal.ornl.gov/): Hyatt D, Locascio PF, Hauser LJ, Uberbacher EC. 2012. Gene and translation initiation site prediction in metagenomic sequences. Bioinformatics 28: 2223–2230.
> - [HMMER](http://hmmer.org/): http://hmmer.org/



run checkM on all files in the directory `fastaOriginal/` that end on .fa and add the ouput into a directoy called `check output/`. 
```
clusterize -j checkM -o checkm.log -n 1 "checkm lineage_wf -t 40 -x fa ./fastaOriginal ./checkMoutput/ -f out_checkM.tab --tab_table"
```

---

<details>
<summary>Click here to see the checkM output table</summary>
 
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

Besides HIMB123 (completness: 86.73), all isolate genomes have a completeness of ≥96.00. However, even 86.73 is sufficient for our purposes.

The highest contamination value is HIMB1427 (4.29%), followed by HIMB1863 (3.32%), HIMB1517 (2,38%), HIMB4 (2.30%), and HIMB1748, HIMB1488, and HIMB123 (all 1.90%), HIMB1505 (1.43%) HIMB1509 and HIMB1506, HIMB1485 (all 1.42 %). All others are below 1% contamination.

We are keeping all isolate genomes, since they are all above 80% completness and below 5% contamination.

### Dereplicating isolate genomes

To dereplicate the isolate genomes, we will be using `dRep`, a python program which performs rapid pair-wise comparison of genome sets.

> dRep citation:
> Olm, M., Brown, C., Brooks, B. et al. dRep: a tool for fast and accurate genomic comparisons that enables improved genome recovery from metagenomes through de-replication. ISME J 11, 2864–2868 (2017). https://doi.org/10.1038/ismej.2017.126

To use dRep, let us first write a bash script telling it to dereplicate at 95% ANI (preclustering at 90% ANI) and then submit it with `clusterize`.

```
nano drep_clusterize.sh
```
content
```
#!/bin/bash 

# Load necessary modules (if applicable)
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
```
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

**51 genomes passed**

---

<details>
<summary>Click to see primary clustering dendrogram</summary>

blue and purple stars: representatives after my dRep

![image](https://github.com/user-attachments/assets/6c33933e-576a-4ffa-a3cc-b00fa7cfaf7b)

</details>

---

### Concatonate FASTAs of dereplicated reference genomes and simplify deflines

We are concatonating all FASTA files into one because we will, further down, be performing competitive read recruitment. And for it to be competitive, all reference genomes need to be in a singular fasta file.

In order to use anvio (which does not like special characters in the deflines of FASTA files) and in order to be able to tell which reference genome is which as we proceed with this analysis, we are simplifying the deflines using anvi'o's `anvi-script-reformat-fasta` program with the flags `--prefix` (for the name of the reference genome) and `--simplify-names` (removing any special characters and simplifying the names), as well as the `--report-file` which will show the mapping from input fasta file name to output fasta file name. In the last step we are concatenating all fasta files with the reformatted deflines into a single `all_fasta.fa`.

do this in the dir in which the fasta files are 
```
# assuming each .fa file is named according to genome name

for g in *.fa; do
  name=${g%.fa}
  anvi-script-reformat-fasta --prefix $name --simplify-names --report-file ${name}-reformat-report.txt -o ${name}_reformatted.fa $g
done

# afterwards, concatenate
cat *_reformatted.fa > all_fasta.fa

```

To check if it worked, we grep for the number of carrots across input files VS in the all_fasta.fa file. The values should match exactly.
```
grep -c ">" *reformatted.fa | awk -F: '{s+=$2} END {print s}'
grep -c ">" all_fasta.fa
```




## Metagenomes

This section explains how to download and quality filter short metagenomic reads from the Hawaii Ocean Time-Series ([Biller et al., 2018](https://doi.org/10.1038/sdata.2018.176), [Mende et al., 2017](https://doi.org/10.1038/s41564-017-0008-3)), the Bermuda Atlantic Time-series Study ([Biller et al. 2018](https://doi.org/10.1038/sdata.2018.176)) OceaTARA Oceans project ([Sunagawa et al., 2015](https://doi.org/10.1126/science.1261359)), the Malaspina Expedition ([Sánchez et al., 2024](https://doi.org/10.1038/s41597-024-02974-1)), the Bio-GO-SHIP project ([Larkin et al., 2021](https://doi.org/10.1038/s41597-021-00889-9)), the bioGEOTRACES project ([Biller et al. 2018](https://doi.org/10.1038/sdata.2018.176)), and the Ocean Sampling Day project ([Kopf et al., 2015](https://doi.org/10.1186/s13742-015-0066-5)).

> [!NOTE]
> All metagenomes we analyzed are publicly available through the European Nucleotide Archive (ENA) repository and NCBI. The project accession numbers are given above in the Summary table.

### Downloading the metagenomes
To download the metagenomes, we will use anvi'o's `[sra_download`](https://anvio.org/help/main/workflows/sra-download/). For this, we need a`SRA_accession_list.txt` for each project and a `download_config.json` config file.

The SRA_accession_list.txt artifact should look like this: no headers, only a list of run accession numbers.
```
$ cat SRA_accession_list.txt
ERR6450080
ERR6450081
SRR5965623
```

#### Prep the download

To create it, we will utilize what we did on metadata curation in [public-marine-omics-metadata](https://github.com/merenlab/public-marine-omics-metadata/tree/main), part of which was the creation of the `SRA_accession_list.tx` artifacts for each project.

For the different projects, they are called: 
```
SRA_accession_PRJEB1787_TARA.txt
SRA_accession_PRJEB8682_OSD.txt
SRA_accession_PRJEB52452_MAL.txt
SRA_accession_PRJNA385854_BGT.txt
SRA_accession_PRJNA385855_BATS.txt
SRA_accession_PRJNA385855_PRJNA352737_HOT_combined.txt
SRA_accession_PRJNA656268_BGS.txt
```

We created a directory for each of the observatories, with a `00_WORKFLOW_FILES` direcotry, to which these files were added.

The `download_config.json` needed for the the `sra_download` we get by running
```
anvi-run-workflow -w sra_download --get-default-config download_config.json
```
In each of the `donwload_config.json` files, we exchanged the `"SRA_accession_list": SRA_accession_list,txt",` with `"SRA_accession_list": "00_WORKFLOW_FILES/SRA_accession_[PROJECT_ACCESSION_NO]_[PROJECT_ACRONYM].txt",` (of course adding the respective project accession number and project acronym).
```
{
    "SRA_accession_list": "00_WORKFLOW_FILES/SRA_accession_[PROJECT_ACCESSION_NO]_[PROJECT_ACRONYM].txt",
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

We can now use this `download_config.json` in combination with anvi'o's `anvi-run-workflow` and the workflow `sra_download` and send it!

do a dry run
```
anvi-run-workflow -w sra_download -c 00_WORKFLOW_FILES/download_config.json -A -n -q
```
and the real deal
```
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

Check if there are any errors in the log files
```
grep -i ERROR sra_download_*.log
```

Check if there is a file for each accession number in the .txt

```
nano check_fastq_files.sh
```
new script
```
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
```
chmod +x check_fastq_files.sh
```
run script
```
./check_fastq_files.sh
```



### QC

We used part of the anvi'o [`metagenomics` workflow](https://anvio.org/help/main/workflows/metagenomics/) to remove noise from raw reads prior to mapping. Namely, the `illumina-utils` program ([Eren et al., 2013](https://doi.org/10.1371/journal.pone.0066643)).

For that, we need a samples-txt file following this structure (group is added bonus if one needs it):
```
column -t samples.txt
sample     group  r1                                           r2
sample_01  G01    three_samples_example/sample-01-R1.fastq.gz  three_samples_example/sample-01-R2.fastq.gz
sample_02  G02    three_samples_example/sample-02-R1.fastq.gz  three_samples_example/sample-02-R2.fastq.gz
sample_03  G02    three_samples_example/sample-03-R1.fastq.gz  three_samples_example/sample-03-R2.fastq.gz
```
It needs to note the sample name and the locations of raw paired-end reads for each sample.

#### Prep QC

Let us create one of those files in two steps

1. Create a more extensive mapping of biosample to run to project to fancy_sample_name

To create said fancy_sample_name, we created a a detailed mapping (which also has the Project accession, Real Sample accession and run accession, as well as the renamed sample names [added based on this schema: prefix: [PROJECT_ACRONYM], then the value sample_accession, then the value in depth, then the value in collection_date, all separated by underscores])

> [!NOTE]
Note, that we are again using files created as part of the [public-marine-omics-metadata](https://github.com/merenlab/public-marine-omics-metadata/tree/main) GitHub repository.

This is what this looked like in practice:

```
nano sampleMapping.py
```
script
```
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

2. Create the sample-txt artifact from that and the files in our 02_FASTA directories

Now, we create the samples-txt artifact. For that we have to account for there being two BioSample_to_SRA_accession_[PROJECT_ACRONYM].csv files in the HOT directory. Further, anvi'o does not appreciate the use of '-' in sample names, so those will be substituted with '_'.

Create and run this in the directory above your project directories.
```
nano makeSamples-txt.py
```
script
```
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

At this step, we will also count how many samples and runs were downloaded
```
nano countSamplesNruns.py
```
script
```
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

#### Do QC
Now, we are ready to get the config file we need:

```
anvi-run-workflow -w metagenomics --get-default-config QC_config.json 
```

We will
- replace `"samples_txt": "samples.txt"` with `"samples_txt": "samples_raw.txt"` 
- turn off everything besides 
   - `idba_ud` (needs to be on to trick snakemake, but we're not acually running it, because when submitting the job we tell it to stop after gzipping), 
   - `iu_filter_quality_minoche` (and set threads to 4), and 
   - `gzip_fastqs`
   - `anvi_script_reformat_fasta` (needs to be on for anvi'o not to complain but is not used here)

---

<details>
<summary>Click here to see the content of the contigs file</summary>

```
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

We do a dry run to look at the stats.

```
anvi-run-workflow -w metagenomics -c QC_config.json -A --until gzip_fastqs -n -q
```

Assumuing the stats look good, we submit the job using clusterize.

> [!NOTE]
> Note that we are using the flag `--unil gzip-fastqs` so that we are not running the entire `metagenomics` workflow but only the quality control portion of that.
```
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
Once it is done, we grep for errors in the log file.
```
grep -i ERROR QC_workflow_*.log
```

#### Make samples_qc.txt

Since we are splitting the metagenomics workflow in two, we are now creating samples_qc.txt, which we will use as the `samples-txt` artifact in the next steps.


```
nano makeQCsample-txt.p`
```
script
```
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

#### create collective samples_qc.txt for all projects

For the upcomming read recruitment, we no longer want to have the projects separated, so we are concatenating the individual samples_qc.txt files we just created into one.
```
nano combine_and_check_samples_qc.py
```
script
```
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

## Mapping metagenomic reads to SAR11 genomes (competitive read recruitment)

This section explains various steps to characterize the occurrence of each SAR11 isolate genome in metagenomes.

For that, we will make use of the anvi'o `metagenomics` snakemake workflow again.

This time, that will include:
-  competitive read recruitment with anvi_run_hmms, anvi_run_kegg_kofams, anvi_run_ncbi_cogs, anvi_run_scg_taxonomy
- convert SAM -> BAM
- profile BAM file to get anvi'o profiles
- merge into single anvi'o profile

Since we now want to map our reads to the reference genomes, we need to now only specify the samples_qc.txt file, but also a `fasta.txt` file:

Create a fasta.txt with the parth to our `all_fasta.fa``file containing all SAR11 reference genomes.
```
name    path
reference_genomes       ../refGenomesSAR11/fastaDrep/dereplicated_genomes/all_fasta.fa
```

get config file
```
anvi-run-workflow -w metagenomics \
                  --get-default-config default-metagenomics-config.json
```

Importantly, we will 
- set `"fasta_txt": "fasta.txt"` (this contains the path to the fasta file containing all reference genomes)
- set `"samples_txt": "samples_qc.txt"` (this contains the paths to all QCed metagenomes)
- put the `reference-mode` to `true` (this way, and by having all reference genomes in a single fasta file, we ensure that the read recruitment is competitive)
- turn on 
   - [`anvi_run_hmm`](https://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-run-hmms) - searches for HMMs against a [contigs-db](https://anvio.org/help/main/artifacts/contigs-db) and stores that information into the contigs-db’s [hmm-hits](https://anvio.org/help/main/artifacts/hmm-hits). 
   - [`anvi_run_kegg_kofams`](https://anvio.org/help/main/programs/anvi-run-kegg-kofams/) - annotates a [contigs-db](https://anvio.org/help/main/artifacts/contigs-db) with HMM hits from KOfam, a database of KEGG Orthologs (KOs). 
   -  [`anvi_run_ncbi_cogs`](https://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-run-ncbi-cogs) - annotate genes in your [contigs-db](https://anvio.org/help/main/artifacts/contigs-db) with [functions](https://anvio.org/help/main/artifacts/functions) from the NCBI’s [Clusters of Orthologus Groups](https://www.ncbi.nlm.nih.gov/COG/)
   - [`anvi_run_scg_taxonomy`](https://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-import-taxonomy)  - finds your single-copy core genes and assigns them taxonomy by searching them against GTDB

Other parts of the config.json file do not have to be specifically turned on, but will run automatically:
- `anvi-gen-contigs-database` - profile all contigs for SAR11 isolate genomes, and generate an anvi’o contigs database that stores for each contig the DNA sequence, GC-content, tetranucleotide frequency, and open reading frames
- `bowtie-build` - map short reads from the project metagenomes onto the scaffolds contained in fasta file containing the reference genomes using Bowtie2
- `samtools_view` - store the recruited reads as BAM files using samtools
- `anvi-profile` - process the BAM files and to generate anvi’o PROFILE databases that contain the coverage and detection statistics of each SAR11 scaffold in a given metagenome
- `anvi-merge` - generate a merged anvi’o profile database from the individual PROFILE databases

This workflow relies on much more than just anvi'o or snakemake, but many other peoples' resources and work. Here are the citations to acknowledge that
- `anvi-gen-contigs-database`: Anvi'o will use 'prodigal' by Hyatt et al. 2010 (doi:10.1186/1471-2105-11-119) to identify open reading frames in your data. 
- Bowtie2: Langmead and Salzberg, 2012
- samtools: Li et al., 2009
- `anvi_run_kegg_kofams`: Anvi'o will annotate your database with the KEGG KOfam database, as described in Aramaki et al (doi:10.1093/bioinformatics/btz859)
- `anvi_run_ncbi_cogs`: Anvi'o will set up the COG20 version of NCBI COGs from Galperin et al. 2021 (https://doi.org/10.1093/nar/gkaa1018)


---

<details>
<summary>Click to see default-metagenomics-config.json</summary>

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

We, again, do a dry run
```
anvi-run-workflow -w metagenomics -c default-metagenomics-config.json -A  -n -q
```
The stats we are getting out reflect that we are working with 1850 metagenomes but only one fasta file containing the reference genomes.
```
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

As all is looking well, we submitted the job again using clusterize:
```
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



## Filter samples based on coverage and detection

In this section, we describe all the steps taken to filter the samples based on coverage and detection. Our cutoff will be to only keep samples, in which at least one of the metagenomes was detected with at least 0.5 detection and 10x coverage.

It will include multiple substeps
- Generating a genomic collection for `anvi-split`
- use `anvi-split` to create individual, self-contained anvi’o projects for each reference genome and its recruited reads
- use `anvi-profile-blitz` to profile the BAM files to get contig-level coverage and detection stats
- filter the samples based on coverage and detection of SAR11 reference genomes

### Generating a genomic collection

At this point anvi’o still doesn’t know how to link scaffolds to each isolate genome. In anvi’o, this kind of knowledge is maintained through ‘collections’. In order to link scaffolds to genomes of origin, we will use the program `anvi-import-collection` to create an anvi’o collection in our merged profile database. This program, however, needs a contigs.txt artifact. 

We will make that by:

For this, we first create a `contigs.txt` (run this in the readRecruitment dir)
```
grep ">" ../refGenomesSAR11/fastaDrep/dereplicated_genomes/all_fasta.fa | sed "s/>//g" > contigs.txt
```
This removes the carrot ">" from the contig headers in our fasta file full of reference genomes by substituting it with nothing globally (not just once but anywhere it sees a ">" in the file) and parses the output into a file called `contigs.txt`. 

These are the names by which the contigs are stored within anvio. 

---

<details>
<summary>Click to see the head of the file</summary>

```
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

Next, we create a file called `bins.txt` which contains only the identifier for the genome.
```
grep ">" ../refGenomesSAR11/fastaDrep/dereplicated_genomes/all_fasta.fa | sed "s/>//g" | cut -d "_" -f 1 > bins.txt
```
This does the same as the command above PLUS it removes the suffix of "_000000000001" from the reference genome identifier: it cuts at the delimiter (`-d`) "_" and keeps the first field (`-f 1`), so the part in front of the delimiter)

---

<details>
<summary>Click to see the head of the file</summary>

```
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

Lastly, we combine the two into a single file and feed that into a new file called `collection.txt`, which is what we need. 
```
paste contigs.txt bins.txt > collection.txt
```
[`collection.txt`](https://anvio.org/help/main/artifacts/collection-txt/) is a two-column TAB-delimited file without a header that describes a [collection](https://anvio.org/help/main/artifacts/collection) by associating items with [bin](https://anvio.org/help/main/artifacts/bin) names.

---
<details>
<summary>Click to see the first 104 lines of the file</summary>


```
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

```
anvi-import-collection collection.txt -p 06_MERGED/reference_genomes/PROFILE.db -c 03_CONTIGS/reference_genomes-contigs.db -C SAR11COLLECTION --contigs-mode
```

We are using `--contigs-mode` because our`collection.txt` describes contigs rather than split.

The collection is not stored in a separate file, but in the PROFILE.db.

### Creating individual, self-contained anvi’o projects for each reference genome and its recruited reads

We will be using `anvi-split` to create individual, self-contained anvi’o projects for each reference genome and its recruited reads.

We are, again, using clusterize to submit the job:
```
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

The output will be put into a directory called `SAR11SPLIT`, which will then contain individual directories with contig.db and profile.db for each reference genome and recruited reads. 

### Determine which samples to continue working with


#### get contig-level coverage and detection stats
We will be using [`anvi-profile-blitz`](https://anvio.org/help/main/programs/anvi-profile-blitz/) to find out which samples are useful (e.g. cutoff of: my genomes should be covered 10x). 

`anvi-profile-blitz` allows the fast profiling of BAM files to get contig- or gene-level coverage and detection stats.

We will give it the `BAM file` that were created in the 04_MAPPING, the `CONTIGS.db` and specify what the output should look like.

This is the base command
```
  
anvi-profile-blitz *.bam \
                   -c contigs-db \
                   -o OUTPUT.txt
```

However, since we are doing it for each `split`, so 51 time (once per reference genome), we are adapting it a bit.



#### combine outputs into one file
Combine the BLITZ outputs into one dataframe
```
cd blitzOUTPUT/
nano combineOutputs.py
```
content
```
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
```
python combineOutputs.py
```


#### Get weighted averages 
Even though some reference genomes have multiple contigs, I want to know how much of each reference genomed was covered in a given sample / how much the entire reference genome was detected (so across contigs): ao the mean_cov and detection information for the collective reference genome. However, we cannot just take the information of all contigs from the same reference genome across a given sample and take the average, but have to take the length of each contig in respective to the length of the entire reference genome (length of its contigs combined). 

```
nano calculate_weighted_coverage.py
```
content
```
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
```
calculate_weighted_coverage.py
```



#### get overview stats of how many samples per project would stay with different detection and coverage cut-offs

How many samples per project make it if we only take the samples for which at least one reference genome has at least 10x coverage and at least 0.5 detection (both need to be true for a sample to pass)? And do that for different coverage and detection combos.

```
nano filter_samples_all_combinations.py
```
content
```
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
```
python filter_samples_all_combinations.py
```

This gives you `combined_stats_weighted_all_filters.txt`

```
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


#### Filter based on 0.5 detection and 10x coverage

```
nano filter_samples_detNcov.py
```
content
```
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
```
python filter_samples_detNcov.py
```

This results in two files called `filtered05N10x_stats.txt` (containing the number of samples of each project that passed the filter criteria) and `filtered05N10x_combined_df.txt` (the sample_IDs, contigs, weighted_mean_cov, weighted_detection of the samples that passed the filtering)

```
sample  reference_genome        project weighted_mean_cov       weighted_detection      total_length
BATS_SAMN07137064_10_2003_04_22 HIMB1436        BATS    17.363279079539108      0.8795767950849176      1208525
BATS_SAMN07137064_10_2003_04_22 HIMB2187        BATS    13.489588372836309      0.8361808067690548      1188940
BATS_SAMN07137064_10_2003_04_22 HIMB2226        BATS    12.678559899501952      0.8206248871785373      1193257
BATS_SAMN07137064_10_2003_04_22 HTCC7211        BATS    21.92   0.9167  1456888
BATS_SAMN07137064_10_2003_04_22 HTCC7214        BATS    11.12   0.7453  1375060
BATS_SAMN07137064_10_2003_04_22 RS39    BATS    17.65   0.8745  1200090
BATS_SAMN07137074_60_2004_03_23 FZCC0015        BATS    13.3    0.8425  1364101
BATS_SAMN07137074_60_2004_03_23 HIMB123 BATS    10.71272055510127       0.8193143926122783      1271840
BATS_SAMN07137074_60_2004_03_23 HIMB1456        BATS    11.999417281233908      0.8549003880194077      1353489
```

### Generate anvi'o PROFILE.db with selected samples

Based on the information we have post-detection0.5-coverage10x-filtering, we know which samples to continue working with. In order to continue working with them, however, we need to generate a new anvi'o PROFILE.db with those samples only.

We will create a .txt file which only lists the sample names (and only lists each once).
```
# Extract the first column (sample names), sort them, and keep only unique entries
awk 'NR>1 {print $1}' ./blitzOUTPUT/filtered05N10x_combined_df.txt | sort | uniq > ./blitzOUTPUT/unique_samples.txt

```

#### Profile the mapping results with anvi’o
First, that requires generating individual PROFILE.dbs for each .bam file of the samples we want. The .bam files are stored in `./04_MAPPING/reference_genomes/`.

To generate these profile.dbs we are using anvi'o's [`anvi-profile`](https://anvio.org/help/7.1/programs/anvi-profile/) program. In its simplest form, that looks like:
```
anvi-profile -i SAMPLE-01.bam -c contigs.db
```

However, we will adapt it to our needs as such:

Submit multiple anvi-profile jobs to a cluster using the clusterize tool. The script reads a list of sample names from a file, processes the corresponding BAM files, and submits the jobs in a controlled manner to avoid overwhelming the cluster scheduler. It includes functionality to stagger job submissions and limit the number of jobs running simultaneously.

We are keeping the number of threads = 16, as described here: https://merenlab.org/2017/03/07/the-new-anvio-profiler/.


```
nano run_anvi_profile_cluster.sh
```
content
```
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
```
chmod +x run_anvi_profile_cluster.sh
```
and run
```
./run_anvi_profile_cluster.sh
```


#### Generate a merged anvi’o profile

Use [`anvi-merge`](https://anvio.org/help/7.1/programs/anvi-merge/) to convert multiple single-profile-dbs (of our selected samples) into a single profile-db (also called a merged profile database). Basically, this takes the alignment data from each sample (each contained in its own single-profile-db) and combines them into a single database that anvi’o can look through more easily.

The basic command is:
```
anvi-merge -c cool_contigs.db \
            Single_profile_db_1 Single_profile_db_2 \
            -o cool_contigs_merge
```

We adjusted it to our directory structure and naming and ran it in `filteredPROFILEdb/`
```
clusterize -j merge_profile_db \
-o merge_profile_db.log \
-n 1 \
"anvi-merge */PROFILE.db -o SAR11-MERGED -c ../03_CONTIGS/reference_genomes-contigs.db"

```


## Create self-contained anvi'o projects for my reference genomes and associated metagenomes, again

We are again using `anvi-split`, but this time, we are using it on the collective PROFILE.db for the selected samples (for which at least one reference genome is above 10x coverage and 0.5 detection) only.

### create COLLECTION, again
First, we have to re-create a COLLECTION that will be stored in the PROFILE.db. For that, we can reuse the collection.txt that we used before, but have to re-create the COLLECTION in itself: specifying the new PROFILE.db.

Here we go:
```
anvi-import-collection collection.txt \
   -p filteredPROFILEdb/SAR11-MERGED/PROFILE.db \
   -c 03_CONTIGS/reference_genomes-contigs.db \
   -C SAR11COLLECTION \
   --contigs-mode
```

### do the splits, again

The output will go into a directory called `SAR11_postFilter`.

```
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



## Visualise 

In this section we will explain how we used `anvi-interactive` to visualise our metagenomes.

For `anvi-interactive` to give us what we want, we need
- to decide which reference genomes to focus on
- a collection and bin to feed to `--gene-mode`
- to prepare and import the metadata we want to order our layers (metagenomes) by in the interactive interface

### decide which SAR11 reference genomes to viusalise

To decide which SAR11 reference genomes to focus on, we will see which are found across most of the projects.

We wrote a little script to output which projects each reference genome is found in:

```
nano countSamplesEachRefGenome.py
```
content
```
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
```
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

> [!NOTE]
OSD only included in NP1 and HTCC1002...

We will focus on those metagenomes that are found in 7 of the 8 projects.


We can also visualise the occurrence (above 10x cov and 0.5 detection) of reference genomes across samples in R

```
## visualise anvi-profile-blitz output

library(ggplot2)
library(dplyr)


blitz_covNdet <- read.table("/Users/rameyer/Documents/_P3/P3dataAnalysis/P3_visualise/outputBLITZ/filtered05N10x_combined_df.txt", header = TRUE, sep = "\t")
head(blitz_covNdet)



# Convert reference_genome to a factor
blitz_covNdet$reference_genome <- as.factor(blitz_covNdet$reference_genome)


# Create a dataframe for project labels
project_labels <- blitz_covNdet %>%
  select(sample, project) %>%
  distinct() %>%
  arrange(sample)

# Create a color palette for projects
project_colors <- project_labels %>%
  distinct(project) %>%
  mutate(color = rainbow(n()))

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
![Uploading file..._ly99umh4u]()






### get collection and bin for `--gene-mode`


First, we start the anvio-dev conda environement
```
conda list env
conda activate anvio-dev
```

Then we look at the [anvi-script-add-default-collection](https://anvio.org/help/main/programs/anvi-script-add-default-collection/) program, which can help us get a collection and bin.

The basic command is:

```
anvi-script-add-default-collection -c [contigs-db](https://anvio.org/help/main/artifacts/contigs-db) \
                                   -p [profile-db](https://anvio.org/help/main/artifacts/profile-db)
```
the program will add a new collection into the profile database named `DEFAULT`, which will contain a single bin that describes all items in the database named `EVERYTHING`. 

I will need to do it by specifying every single PROFILE.db I have (one per ref genome).

```
anvi-script-add-default-collection -c SAR11SPLIT_postFilter/FZCC0015/CONTIGS.db \
                                   -p SAR11SPLIT_postFilter/FZCC0015/PROFILE.db
```

To avoid doing this for every single reference genome manually, we will do it in a loop

```
for dir in SAR11SPLIT_postFilter/*/; do
    contigs_db="${dir}CONTIGS.db"
    profile_db="${dir}PROFILE.db"
    anvi-script-add-default-collection -c "$contigs_db" -p "$profile_db"
done

```

### get gene database

To generate a gene database, anvi'o offers creating one when `anvi-interactive` is started in `--gene-mode`, or alternatively with the program `anvi-gen-gene-level-stats-databases`. We will use the latter here because we need it for many reference genomes.

The basic command is
```
anvi-gen-gene-level-stats-databases -c contigs-db \
                                    -p profile-db \
                                    -C collection
```

Adapted to our data and aim, it is

```
for dir in SAR11SPLIT_postFilter/*/; do
    contigs_db="${dir}CONTIGS.db"
    profile_db="${dir}PROFILE.db"
    anvi-gen-gene-level-stats-databases -c "$contigs_db" -p "$profile_db" -C DEFAULT
done

```

The genes database will automatically be stored in the directory in which the PROFILE.db and CONTIGS.db are also in, but in its own subdirectory called `GENES/`, which contains the `DEFAULT-EVERYTHING.db` GENES database.


### import the metadata we want to order our layers (metagenomes) by

#### Prepare

To visualise the metagenomes capturing each reference genome with their metadata in mind (here, following the latitudinal gradient), we can use `anvi-import-misc-data`.

[`anvi-import-misc-data`](https://anvio.org/help/main/programs/anvi-import-misc-data/) is a program to populate additional data or order tables in pan or profile databases for items and layers, OR additional data in contigs databases for nucleotides and amino acids. We want to use it to propulate additinal data in profile databases for layers.

This blog post gives some more information on this program: https://merenlab.org/2017/12/11/additional-data-tables/#layers-additional-data-table.


Okay, so we need a table to give this program. A table that should look something like this:
```
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


```
anvi-import-misc-data layers_additional_data.txt \
                         -p profile.db \
                         --target-data-table layers
```
and if you want to delete it later 
```
anvi-delete-misc-data -p profile.db \
                      --target-data-table layers 
                      --just-do-it
```
So now we need to make sure that our metadata is following the format that this program expects. That means to use the sample names we used as part of this workflow, but also to subset the main metadata table we made in [public-marine-omics-metadata](https://github.com/merenlab/public-marine-omics-metadata/tree/main) to only keep info on the samples that passed our filtering based on coverage and detection here.

Let's do it!

Our metadata file (metagenomes.txt) currently still has metadata for more samples than the ones that passed the cov and det filtering. 

We will first subset `metagenomes.txt` to only include those samples that passed the filtering and are thus listed in the `unique_samples.txt` we created earlier.

Further, the `metagenomes.txt` file currently has multiple rows per sample if there are multiple runs associated with one sample. We will make it such that there is only one row per sample (so one row per layer that we want to associate this with in anvi'o). Of course, ensuring that any differing values across different runs from the same sample are concatenated, and identical values retained as they are.

```
nano filterMetagenomesTxt.py
```
content
```
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
aggregated_df = filtered_df.groupby('sample').agg(lambda x: ','.join(sorted(set(x.dropna().astype(str)))))

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

Now only select the columns we want to bring into anvi'o

```
nano prepAnvioMetadata.py
```
content
```
import pandas as pd

# Load the dataset
file_path = '../data/metagenomes_filtered.txt'  # replace with your file path
df = pd.read_csv(file_path, sep='\t')

# Select only the desired columns
columns_to_keep = ['sample', 'latitude', 'longitude', 'project', 'depth', 'temperature_degC', 'year', 'model', 'layer', 'season', 'sub_region_seavox', 'region_seavox', 'provdescr_longhurst', 'environment', 'env_biome', 'env_feature', 'env_material']
filtered_df = df[columns_to_keep]

# Save the filtered DataFrame to a new file
filtered_df.to_csv('../data/metagenomes_filtered_anvio.txt', sep='\t', index=False)

```

#### Do it

Since we are using `--gene-mode`, we need to import the metadata into the GENES database.

This is the metadata file `metadata/metagenomes_filtered_anvio.txt` we prepared above.


We will loop through the directories containing the `GENES.db`s for our reference genomes.
```
for dir in SAR11SPLIT_postFilter/*/; do
    genes_db="${dir}GENES/DEFAULT-EVERYTHING.db"
    anvi-import-misc-data metadata/metagenomes_filtered_anvio.txt -p "$genes_db" --target-data-table layers
done
```


Since we want to also show it on the top of the circle (ah, indeed, my lady, I dare say I am most proficient in the distinguished tongue of `anvi'o`), we also need to import it into the PROFILE.db.

We will, again, loop through the directories containing the `PROFILE.db`s for our reference genomes.
```
for dir in SAR11SPLIT_postFilter/*/; do
    profile_db="${dir}PROFILE.db"
    anvi-import-misc-data metadata/metagenomes_filtered_anvio.txt -p "$profile_db" --target-data-table layers
done
```

```
anvi-import-misc-data metadata/metagenomes_filtered_anvio.txt \
                         -p SAR11SPLIT_postFilter/HIMB2204/PROFILE.db \
                         --target-data-table layers

```

### Look at it

To look at it, we are using the `anvi-interactive` command again. With all the prep we did above, it will now show the genes in the metagenomes and allow us to sort layers (metagenomes) based on the metadata keys we imported as well as adjust what will be visualised (coverage, detection, ...) how the genes will be ordered in the metagenomes (synteny, detection, ...) and so on.

We are showing it here with the FZCC0015 reference genome.
```
anvi-interactive -c SAR11SPLIT_postFilter/FZCC0015/CONTIGS.db                  -p SAR11SPLIT_postFilter/FZCC0015/PROFILE.db                  -C DEFAULT                  -b EVERYTHING                  --gene-mode
```

Anvi'o allows us to define a state and export and import that state. 

To do so, we will adjust any parameter we want to adjust and then click on `Save` under the mention of `State` in the left hand panel in the interactive interface. That prompts us to select a name: `visualise`. This will save the state in the GENES.db

If we want to viualise the same project again, we can load that state with the `--state-autoload` flag.
```
anvi-interactive -c SAR11SPLIT_postFilter/FZCC0015/CONTIGS.db \
   -p SAR11SPLIT_postFilter/FZCC0015/PROFILE.db \
   -C DEFAULT \
   -b EVERYTHING \
   --gene-mode \
   --state-autoload visualise
```

However, we can also use the same state with other projects. 

To do so, we will use `anvi-export-state` on the GENES.db in which we defined the state and `anvi-import-state` on the GENES.db that we want to import the state to. After that, we can specify it in the `anvi-interactive` command as above.

Export the state:
```
anvi-export-state -p SAR11SPLIT_postFilter/FZCC0015/GENES/DEFAULT-EVERYTHING.db \
   -s visualise \
   -o visualise.json
```

Import the state to a different project
```
anvi-import-state -p SAR11SPLIT_postFilter/HIMB140/GENES/DEFAULT-EVERYTHING.db \
   -n visualise \
   -s visualise.json 
```

And use `anvi-interactive` to look at it
```
anvi-interactive -c SAR11SPLIT_postFilter/HIMB140/CONTIGS.db \
  -p SAR11SPLIT_postFilter/HIMB140/PROFILE.db \
  -C DEFAULT \
  -b EVERYTHING \
  --gene-mode \
  --state-autoload visualise
```


























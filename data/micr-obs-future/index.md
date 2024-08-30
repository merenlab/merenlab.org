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
- 

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




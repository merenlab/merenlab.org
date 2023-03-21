---
layout: page
title: Ecology of a cryptic plasmid, pBI143
modified: 2021-10-21
excerpt: "by Fogarty et al, 2023"
comments: true
authors: [meren, emily]
---

**The purpose of this page** is to provide access to reproducible data products that underlie our key findings in the study "**A highly conserved and globally prevalent cryptic plasmid is among the most numerous mobile genetic elements in the human gut**" by [Emily Fogarty](https://twitter.com/emily_fogarty11) et al.

{:.notice}
If you have any questions and/or if you are unable to find an important piece of information here, please feel free to leave a comment down below, send an e-mail to [us]({{ site.url }}/people/), or get in touch with us through Discord:

{% include _join-anvio-discord.html %}

## Reproducible / Reusable Data Products

The following data items are compatible with [anvi'o](https://anvio.org) version `v7.1` or later. The anvi'o {% include ARTIFACT name="contigs-db" text="contigs databases" %} and {% include ARTIFACT name="profile-db" text="profile databases' " %} in them can be further analyzed using any program in the anvi'o ecosystem, or they can be used to report summary data in flat-text files to be imported into other analysis environments.

* [doi:10.6084/m9.figshare.22308949](https://doi.org/10.6084/m9.figshare.22308949): Metagenomic read recruitment results from 4,515 global gut metagenomes from healthy individuals using **pBI143 Version 1**.
* [doi:10.6084/m9.figshare.22308964](https://doi.org/10.6084/m9.figshare.22308964): Metagenomic read recruitment results from 4,515 global gut metagenomes from healthy individuals using **pBI143 Version 2**.
* [doi:10.6084/m9.figshare.22308967](https://doi.org/10.6084/m9.figshare.22308967): Metagenomic read recruitment results from 4,515 global gut metagenomes from healthy individuals using **pBI143 Version 3**.
* [doi:10.6084/m9.figshare.22308991](https://doi.org/10.6084/m9.figshare.22308991): Metagenomic read recruitment results from 1,419 gut metagenomes from individuals who were diagnosed with IBD using **pBI143 Version 1**.
* [doi:10.6084/m9.figshare.22309003](https://doi.org/10.6084/m9.figshare.22309003): Metagenomic read recruitment results from 1,419 gut metagenomes from individuals who were diagnosed with IBD using **pBI143 Version 2**.
* [doi:10.6084/m9.figshare.22309006](https://doi.org/10.6084/m9.figshare.22309006): Metagenomic read recruitment results from 1,419 gut metagenomes from individuals who were diagnosed with IBD using **pBI143 Version 3**.
* [doi:10.6084/m9.figshare.22309015](https://doi.org/10.6084/m9.figshare.22309015): Anvi'o single-copy core gene taxonomy for healthy gut metagenomes. Each file in this data pack represents the output of {% include PROGRAM name="anvi-estimate-scg-taxonomy" %} for a given ribosomal protein run on **3,036 gut metagenomes generated from healthy individuals** across the globe.
* [doi:10.6084/m9.figshare.22309033](https://doi.org/10.6084/m9.figshare.22309033): Anvi'o single-copy core gene taxonomy for healthy gut metagenomes. Each file in this data pack represents the output of {% include PROGRAM name="anvi-estimate-scg-taxonomy" %} for a given ribosomal protein run on on **1,419 gut metagenomes generated from individuals who were diagnosed with IBD**.
* [doi:10.6084/m9.figshare.22309105](https://doi.org/10.6084/m9.figshare.22309105): Metagenomic read recruitment results from **936 non-human gut metagenomes** using the **pBI143 version 1**. Metagenomes used for this analysis include open and coastal marine metagenomes from the Tara Oceans and Ocean Sampling Day projects, human oral and skin metagenomes from the Human Microbiota Project, primate gut metagenomes, dog metagenomes, and sewage metagenomes around the globe. Their accession numbers are available via the Supplementary Table 1 (alternative_environment tab) of the study.
* [doi:10.6084/m9.figshare.22309114](https://doi.org/10.6084/m9.figshare.22309114): Metagenomic read recruitment results from 4,515 global gut metagenomes from healthy individuals using **the crAssphage genome**.
* [doi:10.6084/m9.figshare.22298215](https://doi.org/10.6084/m9.figshare.22298215): Metagenomic read recruitment results from four mother-infant datasets from Finland, Italy, Sweden, and the United States. These data enabled us to study the population genetics of pBI143 in mothers and their infants to establish an understanding of the vertical transmission of the cryptic plasmid. A reproducible bioinformatics workflow for the analysis that rely upon these files is below.



## Sequences of the cryptic plasmid pBI143

In our study we identified three sequence-discreet version of pBI143 across human gut metagenomes. We used all three versions extensively for our analyses, and the summary figure below shows the distribution of the pBI143 and its versions across the globe:

{% include IMAGE path="images/pBI143_across_the_globe.png" width=80 caption="The draft Figure 1 from our study that depicts pBI143 prevalence and abundance in globally distributed human populations." %}


### pBI143 Version 1

Assembled from a metagenome internally labeled as `CHI0054` (NCBI accession: [SRX023606](https://www.ncbi.nlm.nih.gov/sra/SRX023606%5Baccn%5D)):

```
>pBI143_V1
AAAGCCGCTGGAGATGTCGGATGTCCGTCAATTAACTCCAGCGGCTTGTATTTTGGTAACTCTCCCGTATCGGTCGGCTA
CCGTTTGGTCAGCAACAAAGGTAGTACTTTATTGTCGCAAATCCAAATGTTTTCCCCTTTATTTGTGCATTTACAGAACA
TTTGCTAACTTCGCACCTATGCAGTATGGTGATAGCACACTAGCACCCTAAAGGGTGCGTGTGTAAGGGGACAAGCCCCC
TTAACCCCCTGTCAGCCGGCTGTCACCGGCTGTTATCGATTGTCACAAAAATGTTAACTTATGGGAGCAACAAGTATTCA
TGTACAAGCAGTGAAGCCGGGGAGTGAGATTCACAACTTTAGGGAAAAAGAGTTGGACTATGTTCGTCCCGAACTTAGTC
ATTTGAATGAAAGCTGGGTTGGAGATAGCATCTCCCATCGGCTGGAAAGTGCGAAACAAAGATATCTCGATACGGTTGGG
CAGAAGATGCAGGCTAAAGCCGCACCCATACGAGAGGGAGTAATAGTAATCAAACAAGAAACCACCATGCAGGAACTCCA
GCAGTTTGCCACGGTCTGCAAAGAACGTTTCGGTATCGAAGCATTTCAAATCCATATACACAAGGACGAAGGATACATGA
ACGCAAAGCAATGGACACCTAACCTTCATGCCCATGTAGTTTTCGATTGGACGCAGCCGAATGGGAAAAGTGTGCGTTTA
TCGCGTGATGACATGGCAGAACTCCAAACCATAGCATCTGAAACACTGGGCATGGAACGTGGCGTTTCTTCTGACCGCAA
ACACTTATCGGCTATGCAGTACAAGACCGAATGTGCGAAAGAACAGCTTCAGGAACTATCAAACGATATATCCAGTGCAT
TAGACAAACACAAGGACGTACAAAACCAGCTTCTTCAGCTTCAGAAAGAACTACGTTCCATTGAAACAAAGAAGAACGTC
CAAAAACTCATTTCTAAGGCTTCAGAGAAGTTTTACGGCTTGATTGGAAAAACAGTCAACGATAGGGAAAAGGATGCCTT
AAAAGCCAAAGTAAAGGCTTTAGAGGGGGAAAATGAACAGCTATCCGATAGACTGGGAAAGGCTATACTTGAAAAAGAGC
AAAACGGCACAAAGGCATTCAAAGCCGAGAATGACAAGGAGTATTATAGACAACAAATGGACAACGCAAGGACAACCAGT
AACCGATTAAGGACGGAGAACCAAGAACTAAAGACCGAAACCAAGGAATTGAAAAAAGAACTTGGTAAGATGAAAGACTT
GTTCAATTCCGAACAACTGGAGGCTTTAAGGCATCATTTCCCAAACATATCAAAAGCTATGGAAGAGGGGAAAGACCTAC
TCAAGCAAATCACCAGAAGCAGAGGTTTCGGTATGGGTATGTAACCCCCCCCCGCCCCCAAAGGGGGAGCGACCAAACGG
CAGCCTCTCTCAATGGAGTGTTACGTCGTTCGGATTCGGAAGAAATTTTGCAGGGTTCGGAAAACGGTTGTGTGCTTATT
CTAAAAAAAGTTTATGGTGTGGCGCTCGTCGGGACGGGTGGTCCCGCCCCTTGCCGCTGGGCGCTCCCCGACCGATGATT
TTTAAGTGACTGATTTTGTGCTGTTTTGGGGGTATATTAAGAATAGAAGAAATAGAATAAGTTAAGTACTTGATACACAA
TATAGGGCATTTTCCATATTGGAAATTCTCATTTTCCAATCTAGAAAATACTGATTTCCTAATATAATACTAAACAGGAA
AATAATATTTCCCTTAATATTGTTTTTATGGAAAATAATAGTTTACTTTGTGGAGAATAATATTTCCCAAAAACACATCA
AAATGGAAAATAAAAAAGCAGTTAAGTTAACCGATTTTCAAAAGAACGAAGAAAATCCTTTTATGAAACAAGCTATAGAA
GGTATTGAAAATCATGTTGTTAAAAAGTATAAGAGTAATAGTGGTGGCGATAAGAGAGCCGTAGTAGCTTTAGCCGACAC
TGAAACTGGAGAAGTGTTTAAGACTTCGTTTATCCGTCAAATAGAAGTAGATGAAGAACAATTCACTAAATTGTATCTTT
CTAACTTTGCTGCATTCTTTGACCTATCACAAGCAGCTATTCGGGTTTTTGGTTACTTTATGACCTGCATGAAACCCAAA
AATGATTTAATCATCTTCAATAGAAAAAAATGCCTAGAATATACCAAATACAAAACAGACAAAGCCGTTTATAAAGGACT
TGCAGAACTTGTAAAAGCTGAAATCATAGCCCGAGGACCAGCCGATAATCTTTGGTTTATTAATCCTCTGATAGTATTCA
ATGGTGACCGAGTGACATTTGCTAAAACATACGTTCGGAAAAAGACTTTAGCTGCCCAAAAGAAAGAAGAAGCAGAGAAA
CGACAATTATCACTTGGCTTTGATGAACAGTAACACTCCATTGAGTGAAGCTGCCGTTTGGTCGCTCCCCCTTTGGGGGC
GGGGGGGGGATAGATAAAGTTCCTCTATGTAAAGTTATAATGGGGGATGAAAGGCAAGGTTCGCTAACCTTACCCCGATT
GTTTATGCTTCCCAGCATTCAATCTATCCCCGTTTATGGGTTTCTCTTTCACCCTACAAAGATAACCTCATGGGGGAAAA
ATGTCCAAGAATATGGGGTAAAAACTATCAAGTCGGTAGAAAATAAGTATCTTTGAAGTACATTTTTAGTAGAGGTACTG
GTATGCCTAGACGAAGCAGATAGGCAAAAAT
```

### pBI143 Version 2

Assembled from a metagenome internally labeled as `ISR0084` (NCBI accession: [SRR341705](https://www.ncbi.nlm.nih.gov/sra/?term=SRR341705)):

```
>pBI143_V2
AAAGCCGCTGGGGATGTCGGATGTCCGTCAATTAACTCCAGCGGCTTGTATTTTGGTAACTCTCCCGTATCGGTCAGCTA
CCGTTTGGTCAGCAACAAAGGTAGTACTTTATTGTCGCAAATCCAAATGTTTTCCCCTTTATTTGTGCATTTACAAAACA
TTCACTAACTTCGCACCTATGCAGTATGGTGATAGCACACTAGCACCCTAAAGGGTACGTGTGTAAGGGGACAAGCCCCC
TTAACCCCCTGTCAGCCGGCTGTCACCGGCTGTTATCGATTGTCACAAAAATGTTAACTTATGGGAGCAACAAGTATTCA
TGTACAAGCAGTGAAGCCGGGGAGCGAGATTCACAACTTTAGGGAAAAAGAGTTGGACTATGTTCGTCCGGAACTTAGTC
ATTTGAATGAAAGCTGGGTTGGAGATAGCATCTCCCATCGGCTGGAGAGTGCGAAACAGAGATATTTCGATACGGTTGGG
CAGAAGATGCAGACTAAAGCAGCACCCATACGGGAGGGGGTAATAGTAATCAAGCAAGAAACCACCATGCAGGAACTCCA
GCAATTCGCTGCGGTCTGCAAGGAACGTTTCGGTATTGAAGCGTTTCAAATCCATATACACAAGGACGAAGGATACATGA
ACGCAAAGCAATGGACACCTAACCTTCACGCCCATGTAGTTTTCGATTGGACGCAGCCGAACGGGAAGAGTGTGCGTTTA
TCGCGTGATGACATGGCAGAACTCCAAACCATAGCATCTGAAGCACTGGGCATGGAACGTGGTGTTTCTTCTGATCGCAA
ACACTTATCGGCTATGCAGTACAAGACCGAATGTGCGAAAGAACAGCTTCAGGAACTATCCAACGATATATCCAGTGCGT
TAGACAAGCACAAGGACGTGCAAAACCAGCTTCTTCAGCTTCAGAAAGAACTACGTTCCATTGAAACAAAGAAAAACGTC
CAAAAACTCATTTCTAAGGCTTCAGAGAAGTTTTACGGCTTGATAGGGAAAACGGTTAACGATAGGGAAAAAGATGCCTT
AAAAGCCAAAATAAAGGCTTTAGAGGGGGAAAATGAACAGCTATCCGATAGACTGGGAAAGGCTATACTTGAAAAAGAGC
GAAACGGCACAAAGGTATTCAAAGCCGAGAACGACAAGGAGTATTATAGACAGCAAATGGACAATGCAAGGACAACCAGT
AACCGACTAAGAACGGAGAACCAAGAACTAAAGACCGAAACCAAGGAATTGAAAAAAGAACTTGGTAAGATGAAAGACTT
GTTCAATTCCGAACAACTGGAGGCTTTAAGGCATCATTTCCCAAACATATCAAAAGCTATGGAAGAGGGGAAAGACCTAC
TCAAGCAAATCACCAGAAGCAGAGGTTTCGGCATGGGTATGTAACCCCCCCCCGCCCCCAAAGGGGGAGCGACCAAACGG
CAGCTTCACTCAATGGAGTGTTACGCCGTTCGGGTTCGGAAGAAATTTTGCAGGGTTCGAAAAACGGTTGTGTGCTTATT
CGAAAAAAAGTTTATGGTGTGGGCGCTCGTCGGGACGGGTGGTCCCGCCCCTTGCCGCTGGGCGCTCCCCGACCGATGAT
TTTTTAAGTAGTTGATTATATGCCGTTTTTAGGGTATAATAAGAACTTAAGAAATAGAATAAGTTAAGTAGTTGATACTC
AATATAGGGTATTTTCCATATTGGAAAATCTGATTTTCCATATTAGAAAGTCACTCTTTCCTATTAGATGTTAAAAAGGA
AAATATATTTTGCTTTTTTATTGTATTAATGGAAATTTAAATATTCCTTTGTGGAAAATCAAATTATACTTTTTCACCGA
TATATAGAAAGTCAAATTTGCCTAAATAGAATGTTTATGGAAAATAAGAAAGTACTCAAATTAACCGATTTTCAAAAAAA
CGAAGAAAACCCCTTTATGAAACAAGCAATTGAAGGGATTGAAAATCACGTTGTTAAGAAATACAAAAGTAGTACGGGAA
GTGATAAAAAAGCTGTAGTTGCTGTTGCTGATACAGATACAGGAGAAGTATTCAAAACCTCGTTTATCCGACAAATAGAA
GTAGATGAAGAGCAATTTGCAAAACTATACCTTGCGAATTTTGCCGCTTTTTTCAATCTATCACAATCAGCTATTCGGGT
TTTCGGGTATATTCTAACTTGCATGAAGCCTAAGAATGATATGATAATTTTCGATAGGCGAAAATGCTTAGAATACACCC
AATACAAATCAGACAAAGCCATATACAAAGGATTGGCAGAGCTCGTACAAAGTGAAATCATAGCAAGAGGTCCAAACGAA
TATAATTGGTTCATTAATCCGTTGATTGTCTTTAATGGGGATAGGGTTTCATTTACAAAAACATACGTTCGGAAAAAGAC
TTTAGCTGCCCAAAAGAAAGAAGAAGCAGAGAAACGACAGTTATCACTTGGTTTTGATGAACTGTAACACTCCATTGAGT
GAAGCTGCCGTTTGGTCGCTCCCCCTTTGGGGGCGGGGGGGGATAGATAAAGTTCCTCTATGTAAAGTTATAATGGGGGA
TGAAAGGCAAGGTTCGCTAACCTTACCCCGATAGTTTATGCTTCCCAGCATTCAATCTATCCCCGTTTATGGGTTTCTCT
TTCACCCTACAAAGATAACCTCATGGGGGAAAAATGTCCAAGAATATGGGGTAAAAACTATCAAGTCGGTAGAAAATAAG
TATCTTTGAAGTACATTTTTAGTAGAGGTACTGGTATGCCTAGACGAAGCAGATAGGCAAAAAT
```

### pBI143 Version 3

Assembled from a metagenome internally labeled as `ISR0084` (NCBI accession: [ERR1136776](https://www.ncbi.nlm.nih.gov/sra/?term=ERR1136776)):

```
>pBI143_V3
AAAGCCGCTGGAGATGTCGGATATCCGTCAATTAACTCCAGCGGCTTGTATTTTGGTAACTCTCCCGTATCGGTCAGCTA
CCGTTTGGTCAGCAACAAAGGTAATGCTTTATCATCGCAATTCCAAATGTTTTTCCCTTTATTTGTGCATTTACAAAACA
TTCACTAACTTCGCACCTATGCAGTATGGTGATAGCGCACTAGCACCCTAAAGGGTGCGTGTGTAAGGGGACAAGCCCCC
TTAACCCCCTGTCAGCCGGCTGTCACCGGCTGTTATCGATTGTCACAAAATGTTAACTTATGGGAGCAACAAGTATTCAT
GTACAAGCAGTGAAGCCGGGGAGCGAGATTCACAACTTTAGGGAAAAAGAGTTGGACTATGTTCGTCCGGAACTTAGTCA
TTTGAATGAAAGCTGGGTTGGAGATAGCATCTCCCATCGGCTGGAGAGTGCGAAACAGAGATATTTCGATACGGTTGGGC
AGAAGATGCAGACTAAAGCAGCACCCATACGGGAGGGGGTAATAGTAATCAAGCAAGAAACCACCATGCAGGAACTCCAG
CAATTCGCTGCGGTCTGCAAGGAACGTTTCGGTATTGAAGCGTTTCAAATCCATATACACAAGGACGAAGGATACATGAA
CGCAAAGCAATGGACACCTAACCTTCACGCCCATGTAGTTTTCGATTGGACGCAGCCGAACGGGAAGAGTGTGCGCTTAT
CGCGTGATGACATGGCAGAACTCCAAACCATAGCATCTGAAGCACTGGGCATGGAACGTGGTGTTTCTTCTGACCGCAAA
CACTTATCGGCTATGCAGTACAAGACCGAATGTGCGAAAGAACAGCTTCAGGAACTATCCAACGATATATCCAGTGCGTT
AGACAAGCACAAGGACGTGCAAAACCAGCTTCTTCAGCTTCAGAAAGAACTACGTTCCATTGAAACAAAGAAAAACGTCC
AAAAACTCATTTCTAAGGCTTCAGAGAAGTTTTACGGCTTGATAGGGAAAACGGTTAACGATAGGGAAAAAGATGCCTTA
AAAGCCAAAATAAAGGCTTTAGAGGGGGAAAATGAACAGCTATCCGATAGACTGGGAAAGGCTATACTTGAAAAAGAGCG
AAACGGCACAAAGGTATTCAAAGCCGAGAACGACAAAGAGTATTATAGACAGCAAATGGACAATGCAAGGACAACCAGTA
ACCGACTAAGAACGGAGAACCAAGAACTAAAGACCGAAACCAAGGAATTGAAAAAAGAACTTGGTAAAATGAAAGACTTG
TTCAATTCCGAACAACTGGAGACTTTAAGACATCATTTCCCGAACATATCGAAAGCTATGGAAGAGGGGAAAGACCTACT
CAAGCAAATCACCAAAAGCAGGGGTTTCGGCATGGGTATGTAACCCCCCCCGCCCCCAAAGGGGGAGCGACCAAACGGCA
GCTTCACTCAATGGAGTGTTACGCCGTTCGGATTCAGAAACTTCACTCAATGGAGTGTTACGCCGTTCGGATTCAGAAAA
CATTTTGCAGGGTTCGAAAAACGGTTGTGTGCTTATCCGAAAAAAAGTTTATGGTGTGGGCGCTCGTCGGGACGGGTGGT
CCCGCCCCTTGCCGCTGGGCGCTCCCCGACCGATGATTTTTAAGTGGCTGATTTTGTGCTGTTTTGGGGGTGTATTAAGA
ACTATAAGAAACAACATAAGCTAAGTAATTGATAATCAATATAGGCTAGTTCACACTGTATGAACTGCCAGTTCAATCTA
GCTTGAACTCAGAGTTCAATAAATCATGTTAAAAATCGAACTGTGTGTTTATTTTTGTTGGACTGATAGTCTATATTTGT
AGGCGAAAAATCGAACTAATAGAATAATATACCATGGCAAAAAACGAACTAAAATTATCCGACTTCCAAAGGAATAAGGA
AAATCCATTCATGAAGCAAGCAATAGAAGATATTGAAAATCATGTTGTTAAAAAATATAAGAGCAGTACAGGAAGTGACA
AAAAAGCAGTTGTTGCCGTTGCTGATACAGATACAGGAGAAGTGTTCCGAACTTCCTTTATTCGACAGATAGAAGTAGAT
GAAGAACAATTTGCAAAACTCTATCTGAACAATTTTGCTGCATTTTTTGACCTTTCACAAGCCGCAATAAGGGTATTTGG
CTATATCATGACTTGTATGAAACCCAAGAATGATATGATAATGTTTATCCTAGAAGATTGTCTTGAATATACGAAATATA
CAAGCAAAGGTACAGTTTATCGAGGGCTTGCAGAACTGGTTAAAGCGGAGATTATAGCCAGAGGTATAAACGAGAACTTA
TGGTTCATCAATCCCCTTATAGTTTTCAATGGGGATAGAGTTTCTTTTACCAAAACATACGTTCGGAAAAAATCTCTATC
AACTAAAAAGAAAGCTGAAGAAGATGAACGCCAGCTATCTTTAGGTTTCAGTTCAGAGTAACACTCCATTGAGCGAAGCT
GCCGTTTGGTCGCTTCCCCCCTTGGGGGACGGGAGGGGGATAGGATAAAGTTCCTCTATGTAAAGTTATAATGGGGGATA
AAGGCAAGGTTCGCTAACCTTACCCCGATTGTCTATGCTTCCCAGCATCCAATCTACCCCCGTTTATGGGTTTCTCTTTC
ACTTTACAAAGATAACCTCATGGGGGAAAAATGTCCAAGAATATGGGGGAAAATCTATCAAATCAGTAGAAAATTAGTAT
CTTTGGGGGACACTTCTAATGGGGGTACTGGTATAGCCTAGACGAAGCGGATAGGCGAAAAT
```

## Analyzing mother-infant metagenomes to find evidence for plasmid transfer through single-nucleotide variants (SNVs)

The purpose of the following workflow is to use read recruitment results obtained from multiple mother-infant gut metagenomes using pBI143 Version 1 sequence to investigate whether unique SNVs occur primarily between family members.

If you are planning to reproduce this workflow, I would suggest you to download [this jupyter notebook](files/pBI143_SNVs.ipynb) file on your computer, and follow it from within your local jupyter environment. If you wish to do that, all you need to do is the following, assuming you are in an `anvio-dev` installed environment (the installation instructions for anvio-dev is [here](https://anvio.org/install/)):

* Download [this file](files/pBI143_SNVs.ipynb) on your computer.
* In your terminal go to the directory in which you have downloaded the file.
* In your terminal type `jupyter notebook`, and select `pBI143_SNVs.ipynb` to start.

Alternatively, you can read through the following steps and look at the output files the workflow reports.

### Acquiring the data pack

The primary input for this analysis is anvi'o project files for four mother-infant datasets from Finland, Italy, Sweden, and the United States. The generation of these project files is very straightforward with the program {% include PROGRAM name="anvi-run-workflow" %}, to which we essentially provided the plasmid sequence and a list of metagenomes of interest. The program {% include PROGRAM name="anvi-run-workflow" %} simply (1) recruited reads using the plasmid sequence from each metagenome it was given, (2) profiled each read recruitment result using the program {% include PROGRAM name="anvi-profile" %}, and finally (3) merged all single profile databases using the program {% include PROGRAM name="anvi-merge" %}. We then moved these final anvi'o project files into four directories (which we will download for reproducibility in a second), where each directory contained a single {% include ARTIFACT name="contigs-db" %} and a single merged {% include ARTIFACT name="profile-db" %} that will enable us to perform the analyses down below.

Let's start with downloading the data pack, which is at [doi:10.6084/m9.figshare.22298215](https://doi.org/10.6084/m9.figshare.22298215). In this jupyter notebook environment, I will use an anvi'o function to directly download it to my work directory:


```python
# from anvi'o utils library import necessary functions
from anvio.utils import download_file, gzip_decompress_file, tar_extract_file

# instruct jupyter notebook to download the datapack
download_file('https://figshare.com/ndownloader/files/39659905',
              output_file_path='MOTHER_INFANT_pBI143_POP_GEN.tar.gz')

# if we are here, the download is finished. now we will
# first decompress it (and get rid of the original file
# while at it).
gzip_decompress_file('MOTHER_INFANT_pBI143_POP_GEN.tar.gz',
                     keep_original=False)

# and finally untar it so we have a clean directory:
tar_extract_file('MOTHER_INFANT_pBI143_POP_GEN.tar',
                 output_file_path='.',
                 keep_original=False)
```

Running the lines above must have created a new directory called `MOTHER_INFANT_pBI143_POP_GEN` in our work directory. Run the `ls` command to confirm:


```python
ls
```

    MOTHER_INFANT_pBI143_POP_GEN/


Where the data pack directory should contain four folders as promised:


```python
ls MOTHER_INFANT_pBI143_POP_GEN
```

    README.txt  fin/  ita/	swe/  usa/


With anvi'o contigs-db and profile-db files in them:


```python
ls MOTHER_INFANT_pBI143_POP_GEN/*/
```

    MOTHER_INFANT_pBI143_POP_GEN/fin/:
    AUXILIARY-DATA.db  CONTIGS.db  PROFILE.db
    
    MOTHER_INFANT_pBI143_POP_GEN/ita/:
    AUXILIARY-DATA.db  CONTIGS.db  PROFILE.db
    
    MOTHER_INFANT_pBI143_POP_GEN/swe/:
    AUXILIARY-DATA.db  CONTIGS.db  PROFILE.db
    
    MOTHER_INFANT_pBI143_POP_GEN/usa/:
    AUXILIARY-DATA.db  CONTIGS.db  PROFILE.db


Good. Now we have the raw recruitment results described as anvi'o project files. Now we will change our work directory to the root of the data pack, and start playing with these data:


```python
import os
os.chdir('MOTHER_INFANT_pBI143_POP_GEN')
```

### Analysis

Unlike the vast majority of data analyses we do with anvi'o using [anvi'o programs](https://anvio.org/help), here we will perform this analysis by accessing anvi'o libraries directly from within Python scripts we will implement below. Many of these steps use anvi'o functionality that is accessible through two anvi'o programs:

* {% include PROGRAM name="anvi-gen-variability-profile" %} (which is a program that gives us access to nucleotide, codon, or amino acid variants in metagenomic read recruitment results),
* {% include PROGRAM name="anvi-gen-variability-network" %}  (which is a program that turns an anvi'o nucleotide variability report into a Gephi compatible XML network file).

But by accessing anvi'o libraries directly, we get to apply filtering rules that are dependent on sample names, such as removing infants from the dataset whose mother does not have a metagenome and vice versa. One could indeed implement those steps after getting the necessary output files from {% include PROGRAM name="anvi-gen-variability-profile" %} using R, EXCEL or by manually selecting samples to be considered for the analysis. But a Pythonic approach helps with reproducibility, and reduces human error.

The actual bottom line is the following: these analyses can be done with anvi'o without writing a single line of Python code, too. In addition to reproducibility, a side purpose of this workflow is to show those who might be interested in exploring anvi'o deeper what else can be done with it. So here we are.

Our analysis starts with importing some libraries that will be necessary later (of course, at the beginning of the analysis we didn't know which libraries were necessary, but we expanded this section as we made progress with the code):


```python
import argparse
import pandas as pd

import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

from anvio.variabilityops import NucleotidesEngine
from anvio.utils import store_dataframe_as_TAB_delimited_file

run = terminal.Run(width=27)
```

Here we define a Python function that takes a  set of anvi'o proifle-db and contigs-db files, and returns a comprehensive dictionary for the nucleotide variation per position:


```python
def get_snvs(contigs_db_path, profile_db_path):
    args = argparse.Namespace(contigs_db=contigs_db_path,
                              profile_db=profile_db_path,
                              gene_caller_ids='0',
                              compute_gene_coverage_stats=True)

    n = NucleotidesEngine(args, r=terminal.Run(verbose=False), p=terminal.Progress(verbose=False))
    n.process()

    return n.data
```

Please note that this is done by passing a set of arguments to the `NucleotidesEngine` class that we imported from `anvio.variabilityops`. We learned what to import from ths library and how to format the `args` from the source code of {% include PROGRAM name="anvi-gen-variability-profile" %} program. Please also note the `gene_caller_ids=0` directive among the arguments passed to the class. The gene caller id `0` corresponds to the mobA gene in pBI143. Since our intention was to focus on mobA, here we limit all data we will acquire from the `NucleotidesEngine` class to that gene.

The next function is a more complex one, but its sole purpose is to work with the sample names (which include information such as the unique family identifier, or whether a given sample belongs to the mother or an infant, etc) to summarize information about the mother-infant pairs remaining in our dataset, which is passed to the function as `udf`:


```python
def mother_infant_pair_summary(udf, drop_incomplete_families=False):
    # learn all the sample names that are present in the dataframe
    sample_names = list(udf.sample_id.unique())

    def get_countries(family_names):
        # learn which countries are present in the datsaet by
        # splitting that piece of information from sample names
        country_names = [s.split('_')[0] for s in family_names]
        return dict(zip(list(country_names), [country_names.count(i) for i in country_names]))

    # we will create this intermediate dictionary to resolve sample names into
    # some meaningful information about mother-infant pairs
    mother_infant_pairs = {}
    for sample_name in sample_names:
        family_name = '_'.join(sample_name.split('_')[0:2])[:-1]
        participant = sample_name.split('_')[1][-1]
        day = sample_name.split('_')[2]

        if family_name not in mother_infant_pairs:
            mother_infant_pairs[family_name] = {'M': [], 'C': []}

        mother_infant_pairs[family_name][participant].append(day)

    # now we can learn about family names that include at least one mother
    # and one infant:
    family_names_with_both_M_and_C = [p for p in mother_infant_pairs if mother_infant_pairs[p]['M'] and mother_infant_pairs[p]['C']]

    # now we subset the broader dictionary to have one that only includes
    # complete families (this is largely for reporting purposes)
    mother_infant_pairs_with_both_mother_and_infant = {}
    for family_name in family_names_with_both_M_and_C:
        mother_infant_pairs_with_both_mother_and_infant[family_name] = mother_infant_pairs[family_name]

    # report some reports .. these will be printed to the user's screen
    # when this function is called
    run.info('Num entires', f"{len(udf.index)}")
    run.info('Num samples', f"{udf.sample_id.nunique()}")
    run.info('Num families', f"{len(mother_infant_pairs)} / {get_countries(mother_infant_pairs.keys())}")
    run.info('   w/both members', f"{len(mother_infant_pairs_with_both_mother_and_infant)} / {get_countries(mother_infant_pairs_with_both_mother_and_infant.keys())}")

    if drop_incomplete_families:
        # reconstruct the original sample names for complete families:
        sample_names_for_families_with_both_M_and_C = []
        for family_name in mother_infant_pairs_with_both_mother_and_infant:
            for day in mother_infant_pairs_with_both_mother_and_infant[family_name]['M']:
                sample_names_for_families_with_both_M_and_C.append(f"{family_name}M_{day}")
            for day in mother_infant_pairs_with_both_mother_and_infant[family_name]['C']:
                sample_names_for_families_with_both_M_and_C.append(f"{family_name}C_{day}")

        # return a new dataframe after dropping all sample names from incomplete families
        # so we have a dataframe that is clean and reprsent only complete families (sad)
        return udf[udf['sample_id'].isin(sample_names_for_families_with_both_M_and_C)]
    else:
        pass
```

Now it is time to get single-nucleotide variant data for the mother-infant pairs from Finland, Sweden, USA, and Italy using the anvi'o project files we downloaded using the data pack before, using the fancy function `get_snvs` we defined above.

It is important to note that the resulting data frames will contain information only for samples that include at least one SNV in the mobA gene of the plasmid. This means, samples that have no SNVs will not be reported in the following data frames even though they were in the original dataset.


```python
df_fin = get_snvs('fin/CONTIGS.db', 'fin/PROFILE.db')

df_swe = get_snvs('swe/CONTIGS.db', 'swe/PROFILE.db')

df_usa = get_snvs('usa/CONTIGS.db', 'usa/PROFILE.db')

df_ita = get_snvs('ita/CONTIGS.db', 'ita/PROFILE.db')

# combine all data frames
df = pd.concat([df_fin, df_swe, df_usa, df_ita])

# reset the `entry_id` column to make sure each entry has a
# unique identifier:
df['entry_id'] = range(0, len(df))
```

Now we take a quick look at the resulting data frame that combines all data from all samples by simply sending this data frame to the function `mother_infant_pair_summary` implemented above, so you can appreciate its true utility:


```python
mother_infant_pair_summary(df)
```

    Num entires ................: 8183
    Num samples ................: 309
    Num families ...............: 102 / {'FIN': 15, 'SWE': 52, 'USA': 24, 'ITA': 11}
       w/both members ..........: 57 / {'FIN': 3, 'SWE': 36, 'USA': 16, 'ITA': 2}


The data frame `df` is quite a comprehensive one, as anvi'o generates an extremely rich output file to describe variants observed in metagenomes. We can take a very quick look, and go through the columns to have an idea (they scroll right):


```python
df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>entry_id</th>
      <th>sample_id</th>
      <th>split_name</th>
      <th>pos</th>
      <th>pos_in_contig</th>
      <th>corresponding_gene_call</th>
      <th>in_noncoding_gene_call</th>
      <th>in_coding_gene_call</th>
      <th>base_pos_in_codon</th>
      <th>codon_order_in_gene</th>
      <th>coverage</th>
      <th>cov_outlier_in_split</th>
      <th>cov_outlier_in_contig</th>
      <th>departure_from_reference</th>
      <th>competing_nts</th>
      <th>reference</th>
      <th>A</th>
      <th>C</th>
      <th>G</th>
      <th>T</th>
      <th>N</th>
      <th>codon_number</th>
      <th>unique_pos_identifier_str</th>
      <th>gene_length</th>
      <th>unique_pos_identifier</th>
      <th>contig_name</th>
      <th>consensus</th>
      <th>departure_from_consensus</th>
      <th>n2n1ratio</th>
      <th>entropy</th>
      <th>gene_coverage</th>
      <th>non_outlier_gene_coverage</th>
      <th>non_outlier_gene_coverage_std</th>
      <th>mean_normalized_coverage</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>FIN_0018M_000A</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001...</td>
      <td>569</td>
      <td>569</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>3</td>
      <td>89</td>
      <td>60</td>
      <td>0</td>
      <td>0</td>
      <td>1.000000</td>
      <td>TT</td>
      <td>C</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>60</td>
      <td>0</td>
      <td>90</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001...</td>
      <td>1104</td>
      <td>0</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001</td>
      <td>T</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>0.000000</td>
      <td>54.063406</td>
      <td>52.346307</td>
      <td>6.167071</td>
      <td>1.109808</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>FIN_0072M_000A</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001...</td>
      <td>344</td>
      <td>344</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>3</td>
      <td>14</td>
      <td>115</td>
      <td>0</td>
      <td>0</td>
      <td>0.947826</td>
      <td>CT</td>
      <td>T</td>
      <td>0</td>
      <td>109</td>
      <td>0</td>
      <td>6</td>
      <td>0</td>
      <td>15</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001...</td>
      <td>1104</td>
      <td>1</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001</td>
      <td>C</td>
      <td>0.052174</td>
      <td>0.055046</td>
      <td>0.204867</td>
      <td>131.563406</td>
      <td>127.259220</td>
      <td>16.885726</td>
      <td>0.874103</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2</td>
      <td>FIN_0072M_000A</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001...</td>
      <td>389</td>
      <td>389</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>3</td>
      <td>29</td>
      <td>135</td>
      <td>0</td>
      <td>0</td>
      <td>0.955556</td>
      <td>CG</td>
      <td>C</td>
      <td>0</td>
      <td>6</td>
      <td>129</td>
      <td>0</td>
      <td>0</td>
      <td>30</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001...</td>
      <td>1104</td>
      <td>2</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001</td>
      <td>G</td>
      <td>0.044444</td>
      <td>0.046512</td>
      <td>0.181820</td>
      <td>131.563406</td>
      <td>127.259220</td>
      <td>16.885726</td>
      <td>1.026121</td>
    </tr>
    <tr>
      <th>3</th>
      <td>3</td>
      <td>FIN_0072M_000A</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001...</td>
      <td>446</td>
      <td>446</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>3</td>
      <td>48</td>
      <td>137</td>
      <td>0</td>
      <td>0</td>
      <td>0.970803</td>
      <td>AG</td>
      <td>A</td>
      <td>4</td>
      <td>0</td>
      <td>133</td>
      <td>0</td>
      <td>0</td>
      <td>49</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001...</td>
      <td>1104</td>
      <td>3</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001</td>
      <td>G</td>
      <td>0.029197</td>
      <td>0.030075</td>
      <td>0.131940</td>
      <td>131.563406</td>
      <td>127.259220</td>
      <td>16.885726</td>
      <td>1.041323</td>
    </tr>
    <tr>
      <th>4</th>
      <td>4</td>
      <td>FIN_0072M_000A</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001...</td>
      <td>458</td>
      <td>458</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>3</td>
      <td>52</td>
      <td>114</td>
      <td>0</td>
      <td>0</td>
      <td>0.956140</td>
      <td>AG</td>
      <td>A</td>
      <td>5</td>
      <td>0</td>
      <td>109</td>
      <td>0</td>
      <td>0</td>
      <td>53</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001...</td>
      <td>1104</td>
      <td>4</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001</td>
      <td>G</td>
      <td>0.043860</td>
      <td>0.045872</td>
      <td>0.180022</td>
      <td>131.563406</td>
      <td>127.259220</td>
      <td>16.885726</td>
      <td>0.866502</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>257</th>
      <td>8178</td>
      <td>ITA_006C_001D</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001...</td>
      <td>1400</td>
      <td>1400</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>3</td>
      <td>366</td>
      <td>25</td>
      <td>1</td>
      <td>1</td>
      <td>0.960000</td>
      <td>AG</td>
      <td>G</td>
      <td>24</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>367</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001...</td>
      <td>1104</td>
      <td>89</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001</td>
      <td>A</td>
      <td>0.040000</td>
      <td>0.041667</td>
      <td>0.167944</td>
      <td>134.579710</td>
      <td>134.309496</td>
      <td>29.018131</td>
      <td>0.185764</td>
    </tr>
    <tr>
      <th>258</th>
      <td>8179</td>
      <td>ITA_024M_000D</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001...</td>
      <td>782</td>
      <td>782</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>3</td>
      <td>160</td>
      <td>21</td>
      <td>0</td>
      <td>0</td>
      <td>0.238095</td>
      <td>CT</td>
      <td>C</td>
      <td>0</td>
      <td>16</td>
      <td>0</td>
      <td>5</td>
      <td>0</td>
      <td>161</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001...</td>
      <td>1104</td>
      <td>29</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001</td>
      <td>C</td>
      <td>0.238095</td>
      <td>0.312500</td>
      <td>0.548874</td>
      <td>13.516304</td>
      <td>13.027230</td>
      <td>5.382654</td>
      <td>1.553679</td>
    </tr>
    <tr>
      <th>259</th>
      <td>8180</td>
      <td>ITA_024M_000D</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001...</td>
      <td>860</td>
      <td>860</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>3</td>
      <td>186</td>
      <td>27</td>
      <td>1</td>
      <td>1</td>
      <td>0.185185</td>
      <td>AC</td>
      <td>A</td>
      <td>22</td>
      <td>5</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>187</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001...</td>
      <td>1104</td>
      <td>34</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001</td>
      <td>A</td>
      <td>0.185185</td>
      <td>0.227273</td>
      <td>0.479166</td>
      <td>13.516304</td>
      <td>13.027230</td>
      <td>5.382654</td>
      <td>1.997587</td>
    </tr>
    <tr>
      <th>260</th>
      <td>8181</td>
      <td>ITA_024M_000D</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001...</td>
      <td>1050</td>
      <td>1050</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>250</td>
      <td>15</td>
      <td>0</td>
      <td>0</td>
      <td>0.266667</td>
      <td>AG</td>
      <td>G</td>
      <td>4</td>
      <td>0</td>
      <td>11</td>
      <td>0</td>
      <td>0</td>
      <td>251</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001...</td>
      <td>1104</td>
      <td>47</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001</td>
      <td>G</td>
      <td>0.266667</td>
      <td>0.363636</td>
      <td>0.579915</td>
      <td>13.516304</td>
      <td>13.027230</td>
      <td>5.382654</td>
      <td>1.109771</td>
    </tr>
    <tr>
      <th>261</th>
      <td>8182</td>
      <td>ITA_017C_003D</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001...</td>
      <td>654</td>
      <td>654</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>118</td>
      <td>67</td>
      <td>1</td>
      <td>1</td>
      <td>0.089552</td>
      <td>AG</td>
      <td>A</td>
      <td>61</td>
      <td>0</td>
      <td>6</td>
      <td>0</td>
      <td>0</td>
      <td>119</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001...</td>
      <td>1104</td>
      <td>90</td>
      <td>USA0006_updated_plasmid_May_08_18_000000000001</td>
      <td>A</td>
      <td>0.089552</td>
      <td>0.098361</td>
      <td>0.301501</td>
      <td>52.745471</td>
      <td>50.215139</td>
      <td>10.296381</td>
      <td>1.270251</td>
    </tr>
  </tbody>
</table>
<p>8183 rows × 34 columns</p>
</div>



While this is our raw data frame that contains all variants from all samples that had at least one variable nucleotide position in the mobA gene, it needs some cleaning.

For instance, in some samples the mobA gene will have a pretty low coverage to perform a robust analysis of SNV patterns, and such samples should be removed first. Here are some of the samples that have the lowest coverage values:


```python
sorted(df.gene_coverage.unique())[:20]
```

    [3.2454710144927534,
     4.181159420289855,
     4.6231884057971016,
     4.684782608695652,
     5.2853260869565215,
     6.155797101449275,
     6.183876811594203,
     6.631340579710145,
     6.665760869565218,
     6.7644927536231885,
     7.306159420289855,
     7.406702898550725,
     7.54981884057971,
     7.729166666666667,
     7.913949275362318,
     8.008152173913043,
     8.177536231884059,
     9.181159420289855,
     10.353260869565217,
     10.72463768115942]


For a very stringent analysis, here we can drop samples from our data frame where the mobA gene has less than 50X coverage:


```python
dfx = df.drop(df[df.gene_coverage < 50].index)
```

The new data frame `dfx` only contains samples that have enough reads to cover mobA at 50X or more. But as you can imagine, this removal step does not include any logic to 'maintain' families in the dataset. Following the removal of samples based on the coverage of mobA, now there will be some infants without mother samples, and some mother samples with no infants. Since the purpose of this analysis to investigate plasmid transfer through SNVs, we don't need those samples.

And `mother_infant_pair_summary` comes to our rescue once again, as we request this function to further drop samples from `dfx` that belong to 'incomplete' families:


```python
dfx = mother_infant_pair_summary(dfx, drop_incomplete_families=True)
```

    Num entires ................: 6076
    Num samples ................: 245
    Num families ...............: 87 / {'FIN': 12, 'SWE': 48, 'USA': 21, 'ITA': 6}
       w/both members ..........: 49 / {'FIN': 2, 'SWE': 33, 'USA': 13, 'ITA': 1}


In our final dataset, we have 49 families for which at least one mother metagenome and one infant metagenome was present.

At this point we can report this final data frame as a TAB-delimited file, an output that will be identical to an output {% include PROGRAM name="anvi-gen-variability-profile" %} would have generated:


```python
# this is the file name in which we will store all the
# final varaibility data:
variability_profile_path = "pBI143_SNVs.txt"

dfx.reset_index(drop=True, inplace=True)
dfx["entry_id"] = dfx.index

# order by [corresponding_gene_call, codon_order_in_gene]
dfx = dfx.sort_values(by = ["corresponding_gene_call", "codon_order_in_gene"])

# ask anvi'o to store it:
store_dataframe_as_TAB_delimited_file(dfx,
                                      variability_profile_path,
                                      columns=dfx.columns.tolist())
```




    'pBI143_SNVs.txt'



If you wish, you can take a look at it in your terminal.

The next step is to represent the information in this file as a network so we can visualize the relationships between all samples in the dataset with respect to their shared SNVs using Gephi. But first, we would like to generate a dictionary with sample information using sample names, so we can highlight samples that belong to the same family, etc.

For this, we will parse the sample names, and create a dictionary, `sample_information_dict`, to pass to the `VariabiltyNetwork` class below:


```python
# an empty dictionary
sample_information_dict = {}

# get the final sample names from `dfx`
sample_names = list(dfx.sample_id.unique())

# go through each sample name, split it into pieces,
# based on the `_` character, and fill in the
# dictionary
for sample_name in sample_names:
    family_name = '_'.join(sample_name.split('_')[0:2])[:-1]
    participant = sample_name.split('_')[1][-1]
    day = sample_name.split('_')[2]
    country = family_name.split('_')[0]

    sample_information_dict[sample_name] = {'family_name': family_name,
                                            'participant': participant,
                                            'country': country,
                                            'sample_name': sample_name,
                                            'coverage': dfx[dfx.sample_id == sample_name].gene_coverage.tolist()[0]}
```

The resulting is a very simple data structure that looks like this:


```python
sample_information_dict
```

```python
    {'SWE_0100C_360D': {'family_name': 'SWE_0100',
                        'participant': 'C',
                        'country': 'SWE',
                        'sample_name': 'SWE_0100C_360D',
                        'coverage': 209.56159420289856},
     'FIN_0226C_014D': {'family_name': 'FIN_0226',
                        'participant': 'C',
                        'country': 'FIN',
                        'sample_name': 'FIN_0226C_014D',
                        'coverage': 103.84510869565217},
     'FIN_0226C_030D': {'family_name': 'FIN_0226',
                        'participant': 'C',
                        'country': 'FIN',
                        'sample_name': 'FIN_0226C_030D',
                        'coverage': 101.2038043478261},
     'FIN_0226C_090D': {'family_name': 'FIN_0226',
                        'participant': 'C',
                        'country': 'FIN',
                        'sample_name': 'FIN_0226C_090D',
                        'coverage': 58.583333333333336},
     'FIN_0226M_000A': {'family_name': 'FIN_0226',
                        'participant': 'M',
                        'country': 'FIN',
                        'sample_name': 'FIN_0226M_000A',
                        'coverage': 250.69927536231884},
     'FIN_0226M_000B': {'family_name': 'FIN_0226',
                        'participant': 'M',
                        'country': 'FIN',
                        'sample_name': 'FIN_0226M_000B',
                        'coverage': 444.6847826086956},
     'FIN_0226M_090D': {'family_name': 'FIN_0226',
                        'participant': 'M',
                        'country': 'FIN',
                        'sample_name': 'FIN_0226M_090D',
                        'coverage': 551.0742753623189},
     'SWE_0005C_000D': {'family_name': 'SWE_0005',
                        'participant': 'C',
                        'country': 'SWE',
                        'sample_name': 'SWE_0005C_000D',
                        'coverage': 6812.400362318841},
     'SWE_0005C_120D': {'family_name': 'SWE_0005',
                        'participant': 'C',
                        'country': 'SWE',
                        'sample_name': 'SWE_0005C_120D',
                        'coverage': 2413.664855072464},
     'SWE_0005M_000D': {'family_name': 'SWE_0005',
                        'participant': 'M',
                        'country': 'SWE',
                        'sample_name': 'SWE_0005M_000D',
                        'coverage': 315.9284420289855},
     'SWE_0007C_360D': {'family_name': 'SWE_0007',
                        'participant': 'C',
                        'country': 'SWE',
                        'sample_name': 'SWE_0007C_360D',
                        'coverage': 2025.1358695652175},
     'SWE_0007M_000D': {'family_name': 'SWE_0007',
                        'participant': 'M',
                        'country': 'SWE',
                        'sample_name': 'SWE_0007M_000D',
                        'coverage': 552.3215579710145},
     'SWE_0008C_000D': {'family_name': 'SWE_0008',
                        'participant': 'C',
                        'country': 'SWE',
                        'sample_name': 'SWE_0008C_000D',
                        'coverage': 4070.4184782608695},
     'SWE_0008C_120D': {'family_name': 'SWE_0008',
                        'participant': 'C',
                        'country': 'SWE',
                        'sample_name': 'SWE_0008C_120D',
                        'coverage': 3511.5190217391305},
     'SWE_0008C_360D': {'family_name': 'SWE_0008',
                        'participant': 'C',
                        'country': 'SWE',
                        'sample_name': 'SWE_0008C_360D',
                        'coverage': 3447.538043478261},
     'SWE_0008M_000D': {'family_name': 'SWE_0008',
                        'participant': 'M',
                        'country': 'SWE',
                        'sample_name': 'SWE_0008M_000D',
                        'coverage': 256.16847826086956},
     'SWE_0011C_120D': {'family_name': 'SWE_0011',
                        'participant': 'C',
                        'country': 'SWE',
                        'sample_name': 'SWE_0011C_120D',
                        'coverage': 7212.425724637681},
      (...)
}
```

Now it is time to ask anvi'o to convert all this into an XML file:

```python
from anvio.variabilityops import VariabilityNetwork

variability_network_path = "pBI143_SNVs.gexf"

args = argparse.Namespace(input_file=variability_profile_path,
                          include_competing_NTs='noise-robust',
                          output_file=variability_network_path)

variability_network = VariabilityNetwork(args, r=terminal.Run(verbose=False), p=terminal.Progress(verbose=False))
variability_network.samples_information_dict = sample_information_dict
variability_network.generate()
```

You can now open the resulting XML file, [`pBI143_SNVs.gexf`](files/pBI143_SNVs.gexf), in [Gephi](https://gephi.org/) to visualize this network. And that's exactly what we did. We opened the file `pBI143_SNVs.gexf` in Gephi:

{% include IMAGE path="images/pBI143-SNVs-Gephi-01.png" width=80 caption="Upon opening the GEXF file of mother-infant variants in Gephi" %}

Then we run the ForceAtlas2 algorithm with the options to dissuade hubs and LinLog mode with non-overlapping nodes while partitioning all nodes into 'family' bins (a piece of information that was available to Gephi thanks to the `sample_information_dict` we passed to the class `VariabiltyNetwork`) so we can see which samples belong to the same family:

{% include IMAGE path="images/pBI143-SNVs-Gephi-02.png" width=80 caption="Following the convergence of ForceAtlas2 with some minor adjustments" %}

We then switched to the Preview panel to export an SVG,

{% include IMAGE path="images/pBI143-SNVs-Gephi-03.png" width=80 caption="Jumping to the Preview panel to get an SVG with desired visual qualities (such as keeping the sample names big on the figure so we can later identify which sample was the mother sample for a given family, etc)" %}

And we finalized this image for publication using Inkscape.

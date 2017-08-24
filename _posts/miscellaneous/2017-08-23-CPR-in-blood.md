---
layout: post
title: "Candidate Phyla Radiation in Human Blood?"
excerpt: "Astonishing if its true. Astonishing if its false. Astonishments all over."
modified: 2017-08-23
tags: [cpr]
categories: [miscellaneous]
comments: true
authors: [meren]
image:
  feature: http://merenlab.org/images/miscellaneous/2017-08-23-CPR-in-blood/Kowarsky_et_al.png
  display: false
---

{% include _toc.html %}

{% capture images %}{{site.url}}/images/miscellaneous/2017-08-23-CPR-in-blood{% endcapture %}

Twitter is bad. I mostly follow scientists, and end up running into interesting findings from other groups, and end up wanting to take a quick look at their data. Although most of these sessions of procrastination don't end up on the blog, sometimes they do: [1](http://merenlab.org/2015/12/09/musings-over-commamox/){:target="_blank"}, [2](http://merenlab.org/2017/01/03/loki-the-link-archaea-eukaryota/){:target="_blank"}, [3](http://merenlab.org/2017/02/05/zhou-salmonella/){:target="_blank"}). Well, today I run into another interesting paper:

[![BLOOD]({{images}}/Kowarsky_et_al.png)](http://www.pnas.org/content/early/2017/08/21/1707009114.full){:target="_blank"}{:.center-img .width-70}

A very interesting finding that has the potential to start some debates about contamination.

Simply Kowarsky et al. takes hundres of blood samples collected from tens of patients, and use shotgun seqeuncing assembly strategies to recover contigs for cell-free DNA. They remove sequences that match to human genome, and they look into what is there. The authors also validate some of their findings by performing independent bench experiments.

While people were asking each other on Twitter whether or not the findings were driven by contamination, I was looking at the authors' use of the term "*dark matter*" in their manucript and was asking myself [when should I move to the countryside and raise chickens](http://merenlab.org/2017/06/22/microbial-dark-matter/){:target="_blank"}. But that's me.

Since the findings are so interesting, I wanted to take a quick look at their contigs. I was very glad to learn that the authors did make available not only the short reads for metagenomes, but also their assembly results. I thank them very much for their observance of open science practices!

## Downloading and characterizing Kowarsky et al.  contigs

I first downloaded the assembled contigs and created an anvi'o contigs database. Everything took about 2-3 minutes:

``` bash
# download the assembled contigs
wget http://www.pnas.org/content/suppl/2017/08/21/1707009114.DCSupplemental/pnas.1707009114.sd05.txt -O Kowarsky_et_al.fa

# generate an anvi'o contigs database
anvi-gen-contigs-database -f Kowarsky_et_al.fa -o Kowarsky_et_al.db

# run anvi'o's default HMM profiles on this contigs database
# using 10 cores
anvi-run-hmms -c Kowarsky_et_al.db -T 10
```

As the authors mentioned in their study, the output from the `anvi-run-hmms` showed 0 rRNA gene hits. But the same output indicated that there were some bacterial single-copy core genes in the assembly.

Whenever there is a new contigs database, I take a quick look at the occurrence of the number bacterial and archaeal single-copy core genes, since this is a great way to roughly estimate [the number of bacterial or archaeal genomes one should expect to find in an assembly]({% post_url anvio/2015-12-07-predicting-number-of-genomes %}).

Running these two commands,

``` bash
anvi-script-gen_stats_for_single_copy_genes.py Kowarsky_et_al.db
anvi-script-gen_stats_for_single_copy_genes.R Kowarsky_et_al.db.hits Kowarsky_et_al.db.genes
```

tells me that the assembly does not explain any highly complete bacterial genome:

[![BLOOD]({{images}}/scgs.png)]({{images}}/scgs.png){:.center-img .width-70}

Fine.

But a closer inspection via good'ol [squinting the eyes technique](https://www.google.com/search?tbm=isch&sa=1&q=trump+looking+at+the+eclipse+without+protection), reveals the occurrence of a good number of Ribosomal proteins occurring in the assembly (which is also pointed out by the authors in their paper).

## Sequences of ribosomal proteins 

Next, I ask anvi'o to give me amino acid sequences of all ribosomal proteins:

``` bash
# first list all the available gene names in Campbell et al.
# single-copy core gene collection (the output is not shown):
anvi-get-sequences-for-hmm-hits -c Kowarsky_et_al.db \
                                --hmm-source Campbell_et_al \
                                --list-available-gene-names

# and from there I select all the ribosomal genes, and ask
# for their amino acid seqeunces to be stored in 
anvi-get-sequences-for-hmm-hits -c Kowarsky_et_al.db \
                                --hmm-source Campbell_et_al \
                                -o Kowarsky_et_al_SCG_hits.fa \
                                --get-aa-sequences \
                                --gene-names "Ribosom_S12_S23, Ribosomal_L1, Ribosomal_L10, Ribosomal_L11, Ribosomal_L11_N, \
                                              Ribosomal_L12, Ribosomal_L13, Ribosomal_L14, Ribosomal_L16, Ribosomal_L17, \
                                              Ribosomal_L18e, Ribosomal_L18p, Ribosomal_L19, Ribosomal_L2, Ribosomal_L20, \
                                              Ribosomal_L21p, Ribosomal_L22, Ribosomal_L23, Ribosomal_L27, Ribosomal_L28, \
                                              Ribosomal_L29, Ribosomal_L2_C, Ribosomal_L3, Ribosomal_L32p, Ribosomal_L35p, \
                                              Ribosomal_L4, Ribosomal_L5, Ribosomal_L5_C, Ribosomal_L6, Ribosomal_L9_C, \
                                              Ribosomal_L9_N, Ribosomal_S10, Ribosomal_S11, Ribosomal_S13, Ribosomal_S15, \
                                              Ribosomal_S16, Ribosomal_S17, Ribosomal_S18, Ribosomal_S19, Ribosomal_S2, \
                                              Ribosomal_S20p, Ribosomal_S3_C, Ribosomal_S4, Ribosomal_S5, Ribosomal_S5_C, \
                                              Ribosomal_S6, Ribosomal_S7"
```

The resulting FASTA file is [here]({{ site.url }}/files/Kowarsky_et_al_SCG_hits.fa). But I summarized the output for you a bit:

|Found in the assembly|*e*-value|Contig|Hit start|Hit stop|
|:--|:--:|:--:|:--:|:--:|
|Ribosomal_L20|5.5e-36|PR_node_48|6269|6617|
|Ribosomal_L32p|1.1e-09|PR_node_60|1388|1538|
|Ribosomal_L9_N|4.7e-18|PR_node_197|1959|2403|
|Ribosomal_L27|7.3e-34|PR_node_197|2470|2749|
|Ribosomal_S16|8.7e-15|PR_node_199|1751|2105|
|Ribosomal_L13|1.4e-33|PR_node_260|1670|2027|
|Ribosomal_L17|1.9e-36|PR_node_260|2032|2383|
|Ribosomal_L13|7.3e-32|PR_node_287|1680|2406|
|Ribosomal_S2|3.8e-62|PR_node_296|1014|1740|
|Ribosomal_L21p|4.2e-28|PR_node_340|1693|2065|
|Ribosomal_L2|1.1e-28|PR_node_444|1|280|
|Ribosomal_L23|1.3e-21|PR_node_444|325|619|
|Ribosomal_L4|1.1e-60|PR_node_444|628|1270|
|Ribosomal_S7|8.3e-56|PR_node_1196|0|477|
|Ribosomal_S20p|7.7e-15|PR_node_1285|146|410|
|Ribosomal_L10|8.6e-23|PR_node_1290|148|637|
|Ribosomal_L14|5.1e-30|HT_node_2492|11|407|
|Ribosomal_S13|2.9e-25|LT_node_1007|1|250|

Authors mention that they found ribosomal proteins in 14 contigs, which perfectly matches to the table above: there are 14 distinct contigs that carry ribosomal proteins in this analysis as well (it is always a great feeling to find perfect agreement with published stuff).

As you can see from this output, multiple of these hits are occurring in the same contig (i.e., PR_node_260 contains two ribosomal proteins one after another, etc), however, we don't know whether any of these distinct contigs are coming from the same population.

The [supplementary methods](http://www.pnas.org/content/suppl/2017/08/21/1707009114.DCSupplemental/pnas.1707009114.sapp.pdf) suggest that these contigs are assembled independently from each dataset. As a final note on this, in Kowarsky et al.'s notation in the study PR, HT, and LT stand for samples originate from 'preagnent', 'heart transplant', and 'lung transplant' individuals.

So one can say most of the contigs with ribosomal proteins are coming from blood samples of preganent individuals.

Fine.

Since there isn't much information on the taxonomic origins of these bacterial hits in Kowarsky et al.'s study, I decided to focus on these ribosomal proteins a bit more. 

## Searching ribosomal proteins in NCBI's nr database

I simply took the FASTA file of amino acid sequences, went to NCBI and used `blastp` to search these amino acid seqeunces in the NCBI's database of non-redundant protein sequences.

Once the serach was done, I clicked the `Download` menu, and selected `Multiple-file JSON` option.

I had never tried to make sense of JSON files form NCBI BLAST output, so I quickly (and dirtily) wrote the following Python program to get the best hit for each query sequence:


``` python
#!/usr/bin/env python

import sys
import json
import glob

# poor man's whatever:
QUERY    = lambda: hits['BlastOutput2']['report']['results']['search']['query_title'].split('___')[0]
QLEN     = lambda: hits['BlastOutput2']['report']['results']['search']['query_len']
HIT      = lambda: hits['BlastOutput2']['report']['results']['search']['hits'][index]
DESC     = lambda: HIT()['description'][index]
TITLE    = lambda: DESC()['title']
SCINAME  = lambda: DESC()['sciname']
ACC      = lambda: '[%(desc)s](https://www.ncbi.nlm.nih.gov/protein/%(desc)s)' % {'desc': DESC()['accession']}
HSPS     = lambda: HIT()['hsps'][0]
PCTALIGN = lambda: str(HSPS()['align_len']* 100 / QLEN()) + '%'
PCTID    = lambda: str(HSPS()['identity'] * 100 / HSPS()['align_len']) + '%'


print('|'.join(['', 'Found in the assembly', 'Best hit on NCBI', 'Percent alignment', 'Percent identity', 'Accession', '']))
print('|'.join(['', ':--', ':--', ':--:', ':--:', ':--:', '']))

# go through every json file in the directory:
for j in glob.glob('*.json'):
    hits = json.load(open(j))

    # skip the poop file
    if 'BlastOutput2' not in hits:
        continue

    # report the best hit:
    index = 0

    # unless the best hit resolves to a multispecies .. if it does, increment
    # index
    while 1:
        if TITLE().find('MULTISPECIES') == -1:
            break

        index += 1

    print('|'.join(['', QUERY(), SCINAME(), PCTALIGN(), PCTID(), ACC(), '']))
```

When it is run in a directory of JSON files downloaded form the NCBI, this program produces a markdown-formatted output table, so lazy people like myself can copy-paste it on their blogs:

|Found in the assembly|Best hit on NCBI|Percent alignment|Percent identity|Accession|
|:--|:--|:--:|:--:|:--:|
|Ribosomal_S2|Candidatus Taylorbacteria bacterium RIFCSPHIGHO2_01_FULL_46_22b|96%|58%|[OHA18677](https://www.ncbi.nlm.nih.gov/protein/OHA18677){:target="_blank"}|
|Ribosomal_L21p|Parcubacteria (Nomurabacteria) bacterium GW2011_GWD2_39_12|83%|58%|[KKR01970](https://www.ncbi.nlm.nih.gov/protein/KKR01970){:target="_blank"}|
|Ribosomal_L2|Parcubacteria (Nomurabacteria) bacterium GW2011_GWE2_36_115|100%|79%|[KKP91491](https://www.ncbi.nlm.nih.gov/protein/KKP91491){:target="_blank"}|
|Ribosomal_L23|Candidatus Zambryskibacteria bacterium RIFCSPHIGHO2_02_FULL_38_22|97%|52%|[OHA94501](https://www.ncbi.nlm.nih.gov/protein/OHA94501){:target="_blank"}|
|Ribosomal_L4|Candidatus Lloydbacteria bacterium RIFOXYC12_FULL_46_25|95%|63%|[OGZ17743](https://www.ncbi.nlm.nih.gov/protein/OGZ17743){:target="_blank"}|
|Ribosomal_S7|Candidatus Zambryskibacteria bacterium RIFCSPHIGHO2_01_FULL_39_63|98%|72%|[OHA86974](https://www.ncbi.nlm.nih.gov/protein/OHA86974){:target="_blank"}|
|Ribosomal_S20p|Parcubacteria bacterium GW2011_GWA2_47_12|97%|61%|[KKU59139](https://www.ncbi.nlm.nih.gov/protein/KKU59139){:target="_blank"}|
|Ribosomal_L10|Candidatus Zambryskibacteria bacterium RIFCSPHIGHO2_02_38_10.5|99%|60%|[OHA92782](https://www.ncbi.nlm.nih.gov/protein/OHA92782){:target="_blank"}|
|Ribosomal_L20|Candidatus Zambryskibacteria bacterium RIFCSPHIGHO2_01_FULL_43_27|96%|64%|[OHA89065](https://www.ncbi.nlm.nih.gov/protein/OHA89065){:target="_blank"}|
|Ribosomal_L32p|Candidatus Lloydbacteria bacterium RIFCSPHIGHO2_01_FULL_41_20|95%|44%|[OGZ03679](https://www.ncbi.nlm.nih.gov/protein/OGZ03679){:target="_blank"}|
|Ribosomal_L9_N|Parcubacteria bacterium GW2011_GWA2_49_9|99%|55%|[KKW11623](https://www.ncbi.nlm.nih.gov/protein/KKW11623){:target="_blank"}|
|Ribosomal_L27|Parcubacteria (Nomurabacteria) bacterium GW2011_GWF1_34_20|96%|71%|[KKP60887](https://www.ncbi.nlm.nih.gov/protein/KKP60887){:target="_blank"}|
|Ribosomal_S16|Parcubacteria (Nomurabacteria) bacterium GW2011_GWF2_36_126|91%|65%|[KKP93401](https://www.ncbi.nlm.nih.gov/protein/KKP93401){:target="_blank"}|
|Ribosomal_L13|Candidatus Zambryskibacteria bacterium RIFCSPLOWO2_12_FULL_39_23|97%|66%|[OHB12744](https://www.ncbi.nlm.nih.gov/protein/OHB12744){:target="_blank"}|
|Ribosomal_L17|Candidatus Lloydbacteria bacterium RIFCSPHIGHO2_01_FULL_41_20|100%|64%|[OGZ03809](https://www.ncbi.nlm.nih.gov/protein/OGZ03809){:target="_blank"}|
|Ribosomal_L13|Polaromonas naphthalenivorans CJ2|38%|96%|[A1VJJ3](https://www.ncbi.nlm.nih.gov/protein/A1VJJ3){:target="_blank"}|
|Ribosomal_L14|Rhodospirillales bacterium 69-11|81%|36%|[OJW27693](https://www.ncbi.nlm.nih.gov/protein/OJW27693){:target="_blank"}|
|Ribosomal_S13|Enterococcus faecalis TX2141|100%|100%|[EFT89376](https://www.ncbi.nlm.nih.gov/protein/EFT89376){:target="_blank"}|

This is a bit crazy.

Amino acid sequences of these ribosomal proteins have reasaonable hits in the NCBI, and most of these hits are most closely related to Candidate Phyla Radiation genomes released in [Christopher T. Brown et al.'s 2015 paper](https://www.nature.com/nature/journal/v523/n7559/full/nature14486.html){:target="_blank"} and [Karthik Anantharaman](https://twitter.com/KarthikGeomicro) et al.'s [recent publication](https://www.nature.com/articles/ncomms13219){:target="_blank"}. Another interesting things is that except the last two lines in this table, all these genes are coming from contigs found in assemblies of pregnant samples.

---

Before making any suggestions, here are some points to take into consideration:

* All these contigs may be coming from a **single microbial population**. In fact, in my opinion, it is very likely since there are only a small number of ribosomal proteins that occur only once.

* If this is truly a single population, why do taxonomic annotations spread all over the place and not a single hit? The likely explanation could be that it is so distant from any member of the CPR described so far by the Banfield group and others, the most close hits vary a lot for no good reason (I'm sure there is a good name for this phenomenon).

* Just like the fact that there may be a single population, this population may be in a single individual. Becasue CPR genomes are often quite lonley in the sequence space, in our experience, they get assembled very nicely even in very complex samples, or even when they are not quite abundant (i.e., [our recovery of CPRs form TARA Oceans metagenomes](http://merenlab.org/data/2017_Delmont_et_al_HBDs/){:target="_blank"}). The next very important step is to map all metagenomes to these contigs to see whether they recruit reads from multiple individuals. If you would like to do that, you can get all metagenomic short reads from 'pregnant' data set, and map those to the contigs PR_node_60, PR_node_48, PR_node_444, PR_node_340, PR_node_296, PR_node_287, PR_node_260, PR_node_199, PR_node_197, PR_node_1290, PR_node_1285, and PR_node_1196 you will find in [this file](http://www.pnas.org/content/suppl/2017/08/21/1707009114.DCSupplemental/pnas.1707009114.sd05.txt){:target="_blank"}. If there is a single individual, it may be possible for someone to go back to that sample, and sequence it deeper to recover a good CPR genome bin.


But regardless, considering that nothing went wrong with the experimental procedures, sampling, library preparation, sequencing, and assembly steps, and even if we are in fact talking about one microbial population in a single pregnant woman, it is probably safe to say for the first time "**CPR in human blood!**".

There you have it.

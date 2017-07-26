---
layout: post
title: "Recovering ribosomal RNA gene sequences with anvi'o"
excerpt: "Use anvi'o to get rRNA genes out of single-cells, clultivars, metagenome-assembled genomes, or even from entire metagenomic assemblies because why not."
modified: 2017-07-09
tags: [16S]
categories: [anvio]
comments: true
authors: [meren]
image:
    feature: http://merenlab.org/images/anvio/2017-07-09-recovering-rRNAs/16S.png
    display: false
---

{% capture images %}{{site.url}}/images/anvio/2017-07-09-recovering-rRNAs{% endcapture %}

{% include _project-anvio-version.html %}

{: .notice}
This feature is available in anvi'o <b>`v2.3.3` or later</b>. You can learn which version you have on your computer by typing `anvi-profile --version` in your terminal.

{% include _toc.html %}

In most cases, getting ribosomal RNA genes out of isolate or metagneome-assembled genomes is not as straightforward as one would like. Gene callers usually don't perform well when it comes to identifying 16S or 23S rRNA genes, and using primer sequences for these regions is not exactly 2017. We added a new feature in anvi'o that reduces down all these to two command lines:

``` bash
anvi-script-FASTA-to-contigs-db GENOME.fa
anvi-get-sequences-for-hmm-hits -c GENOME.db --hmm-source Ribosomal_RNAs
```

So you can spend more time on the next thing on your bottomless TODO list.


## Preface

We always wanted to be able to show rRNA genes in our anvi'o displays. 

Although genes with many conserved regions such as the 16S rRNA gene often fail to remain in our assemblies, some of our bins did contain these genes occasionally. In those cases, while you work with your bins, of course it would have been lovely to see contigs that happen to contain rRNA genes, and maybe get those sequences and quickly BLAST them, etc.

Thanks to [Torsten Seemann](https://twitter.com/torstenseemann){:target="_blank"}'s previous efforts for [Barrnap](https://github.com/tseemann/barrnap/tree/master/db){:target="_blank"}, we now have an operational HMM collection for rRNA genes.

You will find in this post how to work with the HMM profile Ribosomal RNAs in anvi'o to get rRNA genes from a single FASTA file, or from anvi'o projects. There is also some technical blurb at the end of the post if you have too much time in your hands.

## Getting rRNA gene hits from a FASTA file

Let's assume you have a FASTA file. I will download a *B. fragilis* genome from European Nucleotide Archive just to give a reproducible example, and generate a contigs database for it:

``` bash
# download the genome
wget "http://www.ebi.ac.uk/ena/data/view/CR626927.1,CR626928.1&download=fasta&display=fasta" -O B_fragilis_ATCC_25285.fa

# clean up the crappy defline
anvi-script-reformat-fasta B_fragilis_ATCC_25285.fa \
                           --simplify-names \
                           -o B_fragilis.fa

# create an anvi'o contigs database
anvi-gen-contigs-database -f B_fragilis.fa \
                          -o B_fragilis.db
```

Fine. Now we have a contigs database. The next step is to run HMM profiles, which, if you are using anvi'o `v2.3.3` or later, will include the profile for ribosomal RNAs:

``` bash
$ anvi-run-hmms -c B_fragilis.db --num-threads 10
(...)
HMM Profiling for Ribosomal_RNAs
===============================================
Reference ....................................: Seeman T, https://github.com/tseemann/barrnap
Kind .........................................: Ribosomal_RNAs
Alphabet .....................................: RNA
Context ......................................: CONTIG
Domain .......................................: N\A
HMM model path ...............................: /Users/meren/github/anvio/anvio/data/hmm/Ribosomal_RNAs/genes.hmm.gz
Number of genes ..............................: 12
Number of CPUs will be used for search .......: 10
Temporary work dir ...........................: /var/folders/x5/gt4031w53fs63csv1fp0r_3w0000gn/T/tmpbf6ur4j2
HMM scan output ..............................: /var/folders/x5/gt4031w53fs63csv1fp0r_3w0000gn/T/tmpbf6ur4j2/hmm.output
HMM scan hits ................................: /var/folders/x5/gt4031w53fs63csv1fp0r_3w0000gn/T/tmpbf6ur4j2/hmm.hits
Log file .....................................: /var/folders/x5/gt4031w53fs63csv1fp0r_3w0000gn/T/tmpbf6ur4j2/00_log.txt
Number of raw hits ...........................: 18

Psst. Your fancy HMM profile 'Ribosomal_RNAs' speaking
===============================================
Alright! You just called an HMM profile that runs on contigs. Becuse it is not
working with anvio gene calls directly, the resulting hits will need tobe added
as 'new gene calls' into the contigs database. This is a new feature, and if it
starts screwing things up for you please let us know. Other than that you are
pretty much golden. Carry on.

Contigs with at least one gene call ..........: 2 of 2 (100.0%)
Gene calls added to db .......................: 18 (from source "Ribosomal_RNAs")

(...)
```

Done! Now you can do all sorts of fancy stuff. Let's get the 16S rRNA gene sequences:

``` bash
# listing all HMMs in the contigs database
$ anvi-get-sequences-for-hmm-hits -c B_fragilis.db \
                                  --list-hmm-sources
                                  
* Campbell_et_al [type: singlecopy] [num genes: 139]
* Ribosomal_RNAs [type: Ribosomal_RNAs] [num genes: 12]
* Rinke_et_al [type: singlecopy] [num genes: 162]

# taking a look at the available gene names in Ribosomal_RNAs
$ anvi-get-sequences-for-hmm-hits -c B_fragilis.db \
                                  --hmm-source Ribosomal_RNAs \
                                  --list-available-gene-names
                                  
* Ribosomal_RNAs [type: Ribosomal_RNAs]: Archaeal_23S_rRNA, Archaeal_5S_rRNA,
Archaeal_5_8S_rRNA, Bacterial_16S_rRNA, Bacterial_23S_rRNA, Bacterial_5S_rRNA,
Eukaryotic_28S_rRNA, Eukaryotic_5S_rRNA, Eukaryotic_5_8S_rRNA,
Mitochondrial_12S_rRNA, Mitochondrial_16S_rRNA, Archaeal_16S_rRNA

# getting back sequences for 16S rRNA genes
$ anvi-get-sequences-for-hmm-hits -c B_fragilis.db \
                                  --hmm-source Ribosomal_RNAs \
                                  --gene Bacterial_16S_rRNA \
                                  -o sequences.fa
Auxiliary Data ...............................: Found: B_fragilis.h5 (v. 1)
Contigs DB ...................................: Initialized: B_fragilis.db (v. 8)
Hits .........................................: 18 hits for 1 source(s)
Filtered hits ................................: 6 hits remain after filtering for 1 gene(s)
Mode .........................................: DNA seqeunces
Genes are concatenated .......................: False
Output .......................................: sequences.fa
```
 
So far so good. One of the sequences in the file is this:

```
>Bacterial_16S_rRNA___Ribosomal_RNAs___823984|bin_id:B_fragilis|source:Ribosomal_RNAs|e_value:0|contig:c_000000000001|gene_callers_id:4360|start:3205534|stop:3207058|length:1524
AGGAGGTGTTCCAGCCGCACCTTCCGGTACGGCTACCTTGTTACGACTTAGCCCCAGTCACCAGTTTTACCCTAGGACGATCCTTGCGGTTACGTACTTCAGGTACCCCCGGCTCCCATG
GCTTGACGGGCGGTGTGTACAAGGCCCGGGAACGTATTCACCGCGCCGTGGCTGATGCGCGATTACTAGCGAATCCAGCTTCACGAAGTCGGGTTGCAGACTTCGATCCGAACTGAGAGA
GGATTTTGGGATTAGCATACGGTCACCCGCTAGCTGCCTTCTGTACCCCCCATTGTAACACGTGTGTAGCCCCGGACGTAAGGGCCGTGCTGATTTGACGTCATCCCCACCTTCCTCACA
TCTTACGACGGCAGTCTCTCTAGAGTCCTCAGCATAACCTGTTAGTAACTAAAGATAAGGGTTGCGCTCGTTATGGCACTTAAGCCGACACCTCACGGCACGAGCTGACGACAACCATGC
AGCACCTTCACAGCGGTGATTGCTCACTGACATGTTTCCACATCATTCCACTGCAATTTAAGCCCGGGTAAGGTTCCTCGCGTATCATCGAATTAAACCACATGTTCCTCCGCTTGTGCG
GGCCCCCGTCAATTCCTTTGAGTTTCACCGTTGCCGGCGTACTCCCCAGGTGGAATACTTAATGCTTTCGCTTGGCCGCTTACTGTATATCGCAAACAGCGAGTATTCATCGTTTACTGT
GTGGACTACCAGGGTATCTAATCCTGTTTGATACCCACACTTTCGAGCATCAGTGTCAGTTGCAGTCCAGTGAGCTGCCTTCGCAATCGGAGTTCTTCGTGATATCTAAGCATTTCACCG
CTACACCACGAATTCCGCCCACCTCTACTGTACTCAAGACTGACAGTATCAACTGCAATTTTACGGTTGAGCCGCAAACTTTCACAACTGACTTACCAGTCCACCTACGCTCCCTTTAAA
CCCAATAAATCCGGATAACGCTCGGATCCTCCGTATTACCGCGGCTGCTGGCACGGAGTTAGCCGATCCTTATTCATATAATACATACAAAACAGTATACATACTGCACTTTATTCTTAT
ATAAAAGAAGTTTACGACCCATAGAGCCTTCATCCTTCACGCTACTTGGCTGGTTCAGGCTAGCGCCCATTGACCAATATTCCTCACTGCTGCCTCCCGTAGGAGTTTGGACCGTGTCTC
AGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGAGCCGTTACCTCACCAACAACCTAATGGAACGCATCCCCATCCTTTACCGGAATCCTTTAATAAT
GAAACCATGCGGAATCATTATGCTATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTAAAGGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCAGCAAAGAAAG
CAAGCTTTCTTCCTGATGCCCCTCGACTTGCATGTGTTAAGCCTGTAGCTAGCGTTCATCCTGAGCCAGGATCAAACTCTTCAT
```

And a BLAST search indicates that it matches 100% over its entire lenght to the 'Bacteroides fragilis strain NCTC 9343 16S ribosomal RNA gene, complete sequence'.

As you can see, things may occasionally work even when you are dealing with bioinformatics.


## Working with anvi'o standard files

The beauty is that you don't really need to do anything if you already have been using anvi'o to study your metagenome. All you need to do is to re-run the program `anvi-run-hmms` on your contigs database, and you are golden. When you visualize your metagenomes with the program `anvi-intearctive`, you will automatically see where all the rRNA genes are among your contigs. Or you can access to them via the program `anvi-get-sequences-for-hmm-hits` as demonstrated above.

Here I will download the Infant Gut data pack from [this tutorial]({{site.url}}/tutorials/infant-gut){:target="_blank"} to demonstrate how things run for an already existing anvi'o project. Let's go.

``` bash
# download the data, unpack, and get into the directory:
wget https://ndownloader.figshare.com/files/8252861 -O INFANT-GUT-TUTORIAL.tar.gz
tar -zxvf INFANT-GUT-TUTORIAL.tar.gz
cd INFANT-GUT-TUTORIAL

# run HMMs for the new profile
anvi-run-hmms -c CONTIGS.db --installed-hmm-profile Ribosomal_RNAs
```

Alright! 

Now we run the interactive interface,

``` bash
anvi-interactive -p PROFILE.db -c CONTIGS.db
```

You will see the addition of a new layer on the display:

[![rRNA]({{images}}/infant-gut.png)]({{images}}/infant-gut.png){:.center-img .width-70}

The marks in most outer layer identify contigs with rRNA genes. If you zoom in to one of those, and click 'Inspect' from your right-click menu, you can see them clearly marked in red:

[![16S]({{images}}/16S.png)]({{images}}/16S.png){:.center-img .width-70}

You can get the sequence for it, or go to NCBI with it directly from the interface!

---

We hope this is useful to you, and please feel free to send any questions to our way as usual.


## Under the hood

Here are some technical stuff for archival purposes. You clearly don't need to read this if you are not interested.

### Why not before?

We always wanted to report rRNA genes in contigs, but our major problem was not being able to identify them as proper gene calls (even when they were assembled properly, which is quite problematic within itself). For instance, if you have a 16S rRNA gene somewhere, Prodigal often reports multiple short and inaccurate gene fragments for it. Which makes it impossible to annotate those fragments accurately, or use them for anything.

In his [Barrnap repository](https://github.com/tseemann/barrnap/tree/master/db){:target="_blank"}, [Torsten Seemann](https://twitter.com/torstenseemann){:target="_blank"} has a nice collection of HMMs for rRNAs for bacteria, archaea, and eukarya (thank you very much for putting these together, Torsten!). So we could essentially use his collections to extend anvi'o's HMM profiles, which is quite flexible. We did it, and [added the new profile into our repository](https://github.com/merenlab/anvio/tree/master/anvio/data/hmm/Ribosomal_RNAs){:target="_blank"}, but there were multiple other problems along the way.

### Design issues

We had designed anvi'o's HMM framework to work with gene calls. For instance, when a new contigs database is generated, first Prodigal identifies all the gene calls (unless the user has used `--external-gene-calls` directive to provide their own gene calls), and then AA sequences for gene calls are handed over to another anvi'o module `anvi-run-hmms` uses to search for HMM profiles, say, to identify single-copy core genes. When we wanted to work with rRNA HMMs, we run into three major problems;

1. Torsten's HMM profiles for rRNAs were obviously in RNA alphabet, so we needed to hand over sequences in RNA alphabet to the anvi'o module that was used to working with AA sequences,

2. Because Prodigal couldn't identify those genes properly, we no longer could rely on gene calls to run these HMMs, so instead of genes, we now needed to send contigs to this module,

3. Since we no longer are working with gene calls, we all of a sudden had a problem of 'follow-up': when you hand over a bunch of gene calls with unique IDs to another module, it is much more easier to make sense of the results when they are back, and link them to the information that is already in the database. But working with contigs was going to take that away.

We solved these things in multiple ways without braking other things in the repository.

First, we added a new property to the standard anvi'o HMM profile files: `target.txt`. For [Campbell et al.](https://github.com/merenlab/anvio/tree/master/anvio/data/hmm/Campbell_et_al){:target="_blank"}, which is the good old bacterial single-copy core gene HMMs, the content of `target.txt` is `AA:GENE`, telling anvi'o that it wants *amino acid* sequences for *genes*. For the [new rRNA profile](https://github.com/merenlab/anvio/tree/master/anvio/data/hmm/Ribosomal_RNAs){:target="_blank"}, the same file says `RNA:CONTIG`. As you can gess, this one tells anvi'o that it is supposed ot be working with *contig sequences* in *RNA* alphabet.

We also extended our design to do create new gene calls in the contigs db based on results that are returned from an HMM profile that uses contig sequences. These additions are fully extensible, so if you come up with your own collection of, say, antibiotic-resistance genes, they should work seamlessly. All you need to do is to mimic one of the already-working HMM profiles in the database.

Then we decided to display rRNA gene calls in 'red', instad of 'green' or 'gray', in the inspection pages of the interactive interface. These take some serious coding, but anvi'o is very flexible, and adding features more fun than a nightmare.

### Torsten HMMs to anvi'o HMMs

Just for archival purposes, these are the steps I used to reformat Torsten's files to convert them into a proper anvi'o HMM profile on my Mac computer:


``` bash
# checkout the repo, and make sure you are on the right time
# and right place for everything to work:
cd /tmp
git clone https://github.com/tseemann/barrnap.git
git checkout ebfdc202842b4ec16ac6b3a380b17f2e00ab6b68


# Go into the databases directory. You can double check
# every change that will follow with `git diff`:
cd barrnap/db

# Remove all DESC lines. HMMER converts those spaces
# into TABs, extending number of TAB characters arbitrarily:
sed -i '' '/^DESC / d' *.hmm

# Make those names look better and distinguishable:
sed -i '' 's/NAME  /NAME  Bacterial_/g' bac.hmm
sed -i '' 's/NAME  /NAME  Archaeal_/g' arc.hmm
sed -i '' 's/NAME  /NAME  Mitochondrial_/g' mito.hmm
sed -i '' 's/NAME  /NAME  Eukaryotic_/g' euk.hmm

# Add noise cutoffs. totally arbitrary :/ maybe at some 
# point we will fix this behavior and the HMM framework
# in anvi'o will be able to use HMM profiles without any
# model noise cutoffs. but for now, this is it:
sed -i '' '/CKSUM /a \
GA    750 750;\
TC    750 750;\
NC    750 750;\
' *.hmm

# create a new directory (which will become all the anvi'o data):
mkdir Ribosomal_RNAs

# concatenate all into one file:
cat *.hmm > Ribosomal_RNAs/genes.hmm

# Create a `genes.txt` file:
echo "gene accession hmmsource" > Ribosomal_RNAs/genes.txt

for i in `grep NAME Ribosomal_RNAs/genes.hmm | awk '{print $2}'`
do
    echo "$i None barrnap"
done >> Ribosomal_RNAs/genes.txt

perl -p -i -e 's/ /\t/g' Ribosomal_RNAs/genes.txt

gzip Ribosomal_RNAs/genes.hmm

# create other necessary files:
echo "Ribosomal_RNAs" > Ribosomal_RNAs/kind.txt
echo "Seeman T, https://github.com/tseemann/barrnap" > Ribosomal_RNAs/reference.txt
echo "RNA:CONTIG" > Ribosomal_RNAs/target.txt

# done!
```

We hope these are useful.
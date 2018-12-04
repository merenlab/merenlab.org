---
layout: post
title: "Importing GhostKOALA/KEGG annotations into anvi'o"
excerpt: "KEGG modules, meet anvi'o. Anvi'o, meet KEGG modules."
modified: 2018-01-17
tags: [tutorial, GhostKoala]
categories: [anvio]
comments: true
authors: [elaina]
---

{% capture images %}{{site.url}}/images/anvio/2018-01-17-working-with-GhostKoala{% endcapture %}

{% include _toc.html %}

{:.notice}
**A note from the Meren Lab**: We are very thankful to Elaina for sharing her expertise on behalf of all anvi’o users who wished to import KEGG annotations into their workflows. Elaina is currently a third-year PhD student at the University of Southern California, where she uses her background in microbiology and molecular biology, and her learnings from genome-resolved metagenomics, to generate testable hypotheses regarding the diversity and functioning of microbial populations. She is the developer of [BinSanity](https://peerj.com/articles/3035/), which is an automated binning algorithm that was recently used in another study she co-authored to recover [more than two thousand metagenome-assembled genomes](https://www.nature.com/articles/sdata2017203) from the TARA Oceans Project.

Genome annotation is a key step in analyzing bioinformatic data, but with a variety of available databases it can be difficult to decide where to start. One useful database is the [Kyoto Encyclopedia of Genes and Genomes (KEGG)](https://doi.org/10.1093/nar/gkv1070). KEGG integrates functional information, biological pathways, and sequence similarity. GhostKOALA is an automatic annotation tool to assign KEGG Identifiers to metagenomes. GhostKOALA is a web server though, so is it really worth it to annotate your contigs using this service over a program you can run locally, such as interproscan or the NCBI COGs database? Well one reason to use KEGG is that the database contains a variety of well described and environmentally important metabolic pathways with associated pathway maps available. Another is that GhostKOALA can handle metagenomes containing eukaryotes. KEGG is also widely used as a reference so KEGG identifiers can be easily cross linked in published literature. In my own experience using this workflow I found that making initial observations about the functional roles of the MAGs I generated was expedited because of the ability to quickly take KEGG Identifiers and import them into the BRITE mapping function of KEGG or using something like the [KEGG Decoder](https://github.com/bjtully/BioData/tree/master/KEGGDecoder) to get a visualization of functional potential. 

This tutorial walks you through annotating your contigs with [GhostKOALA](http://www.kegg.jp/ghostkoala/), and compiling the results into a format easily importable to your anvi'o contigs database using `KEGG-to-anvio`. You can download all of the associated scripts for this workflow [here](https://github.com/edgraham/GhostKoalaParser).

All you'll need to run this tutorial is an installation of anvi'o, python, and internet access. In python you'll need the pandas and BioPython modules. I was running [Pandas](https://pandas.pydata.org/) `v0.22.0`, and [BioPython](http://biopython.org/) `v1.70`, both of which are dependencies of anvi'o (hence you should have them if you have anvi'o).

So first lets download all of the scripts and files we'll need for this tutorial. They are all on github so we can clone the repository as shown below:


``` bash
 $ git clone https://github.com/edgraham/GhostKoalaParser.git
```

This will produce a directory containing all of the files and scripts we'll use throughout this tutorial, which assumes that you already have an anvi'o contigs database (as explained in the [metagenomic workflow tutorial](http://merenlab.org/2016/06/22/anvio-tutorial-v2/)).

Now that you have your contigs database, you may want to annotate your gene calls. There are so many ways to do this, but I am going to show you how to pull KEGG annotations into yours using [GhostKOALA](https://doi.org/10.1016/j.jmb.2015.11.006), which is KEGGs metagenome annotating service.

GhostKOALA is one of the few ways I have found to access KEGG gene calls for large datasets without our lab subscribing to use the downloadable KEGG database. But bfore I start I just wanted to recommend that everyone gives [the paper](https://doi.org/10.1016/j.jmb.2015.11.006) a quick glance over as one thing you'll notice is that GhostKOALA will not output bit scores or e-values, and it has a built-in way of determining whether a particular annotation is accurate or not.


### Export Anvio Gene Calls

The first step to this is exporting amino acid sequences from your anvi'o contigs database.

``` bash
 $ anvi-get-sequences-for-gene-calls -c CONTIGS.db \
                                     --get-aa-sequences \
                                     -o protein-sequences.fa
```

This will make sure anvi'o will be able associate annotations back to the genes described in the contigs database.

### Run GhostKOALA

{:.notice}
GhostKOALA does have an upload limit of 300MB. If your anvi'o database is rather large and you find that your `protein-sequences.fa` file is larger than that, you can split the file into multiple FASTA files, and concatenate the GhostKOALA outputs you will be generating during this step.

Before we run GhostKOALA we need to make a few modifications to our 'protein-sequences.fa' file. The parser that GhostKOALA uses for fasta sequences appears to have an unfriendly relationship with sequence id's that begin with a number and will send you angry error messages if you try and submit the file directly exported from anvi'o which looks like the one below:

```
>0
MAEYQNIFTQVQVQGPAEMGVDPAGTLSRERTNGTSFSKLAGLFGNAQLGPIYLGTFGLI
SLVTGFAWFFMVGLSFWDQVDYSPALFLRELFWLALEPPAEEYGLSIPPMAEGGYFLLAS
FFLLISVICWWVRTYLRAEELGMGKHVAWAFASAIWLFLVLGLFRPILMGSWSEMVPYGI
FPHLDWTNLFSLTYGNLFYNPFHALSIVFLYGSALLFAMHGATILAVSRYGGEREIEQIV
DRGTASERAALFWRWTM
>1
TYLRAEELGMGKHVAWAFASAIWLFLVLGLFRPILMGSWSEMVPYGI
FPHLDWTNLFSLTYGNLFYNPFHALSIVFLYGSALLFAMHGATILAVSRYGGEREIEQIV
DRGTASERAALFWRWTMGFNATMEGIHRWAWWFAVLTTLTGGIGILLTGTVVDNWFIWAQ
DHGYAPLN
(...)
```

To remedy this, run this command on your FASTA file to add `genecall_` prefix to every defline in your FASTA file:

```bash
 $ sed -i '' '%s/\>/\>genecall_/g' test.fa
```

Your FASTA file is now ready to be submitted to GhostKOALA for annotation.

To do this go to the [GhostKOALA webserver](http://www.kegg.jp/ghostkoala/), and click the "Choose File" button underneath the section that says **Upload query amino acid sequences in FASTA format**. From the menu you will upload your `protein-sequences.fa`, and will be asked to provide an email address.

{:.notice}
You can only run one instance of GhostKOALA at a time per email address.

Before the run starts you'll get an email that will look like this:

```
Your GhostKOALA job request
Query dataset: 1 entries
KEGG database to be searched: c_family_euk+genus_prok

Please click on the link below to either submit or cancel your job.

   (you will see your links here)

If no action is taken within 24 hours, your request will be deleted.
```

Click the link that says 'Submit', and your run will begin processing.

{:.notice}
The KEGG parser script also has the option to combine your GhostKOALA results with interproscan. If you want to incorporate both annotations from interproscan and KEGG in the same table, follow the tutorial [here](http://merenlab.org/2016/06/18/importing-functions/) to run interproscan, **but run interproscan with the flags** `-f tsv`, `--goterms`, `--iprlookup`, and `--pathways`. 

### Generate the KEGG orthology table

{:.notice}
This step is **optional** and is here for anyone who is curious about how I generated the KEGG Orthology file :)

In the repository you cloned earlier there is a file called `KO_Orthology_ko00001.txt`.

If you take a peak at this file it looks like this:

|Metabolism|Overview|01200 Carbon metabolism [PATH:ko01200]|K00844  HK; hexokinase [EC:2.7.1.1]|
|Metabolism|Overview|01200 Carbon metabolism [PATH:ko01200]|K12407  GCK; glucokinase [EC:2.7.1.2]|
|Metabolism|Overview|01200 Carbon metabolism [PATH:ko01200]|K00845  glk; glucokinase [EC:2.7.1.2]|
|Metabolism|Overview|01200 Carbon metabolism [PATH:ko01200]|K00886  ppgK; polyphosphate glucokinase [EC:2.7.1.63]|
|Metabolism|Overview|01200 Carbon metabolism [PATH:ko01200]|K08074  ADPGK; ADP-dependent glucokinase [EC:2.7.1.147]|
|Metabolism|Overview|01200 Carbon metabolism [PATH:ko01200]|K00918  pfkC; ADP-dependent phosphofructokinase/glucokinase [EC:2.7.1.146 2.7.1.147]|
|Metabolism|Overview|01200 Carbon metabolism [PATH:ko01200]|K01810  GPI, pgi; glucose-6-phosphate isomerase [EC:5.3.1.9]|
|Metabolism|Overview|01200 Carbon metabolism [PATH:ko01200]|K06859  pgi1; glucose-6-phosphate isomerase, archaeal [EC:5.3.1.9]|
|(...)|(...)|(...)|

I will now show you how to generate this file, which contains the necessary information to convert the KEGG Orthology assignments to functions.

Because the KEGG database is currently working under a subscription model, I had to find a workaround to access the information to match the orthologies with function. To do this you can run the following command (or click the Download htext link in [this page](http://www.genome.jp/kegg-bin/get_htext?ko00001.keg)) to download the htext file:

```
wget 'https://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=htext&filedir=' -O ko00001.keg
```

Once it is finished, you can use this (not-so-beautiful) code snippet to parse that file into the table above. It should work if you simply copy-paste it into your terminal and have the file `ko00001.keg` in your working directory:

``` bash
kegfile="ko00001.keg"

while read -r prefix content
do
    case "$prefix" in A) col1="$content";; \
                      B) col2="$content" ;; \
                      C) col3="$content";; \
                      D) echo -e "$col1\t$col2\t$col3\t$content";;
    esac 
done < <(sed '/^[#!+]/d;s/<[^>]*>//g;s/^./& /' < "$kegfile") > KO_Orthology_ko00001.txt
```

This simply goes through the hierarchical ‘.keg’ file you downloaded, and extracts different layers. The resulting output file should be a tab delimited file where the first column corresponds to the broadest classification, the fifth corresponds to the gene itself.

You may notice that some of the keg identifiers appear multiple times in this parsed folder. This is because many of the metabolism genes are constituents of multiple pathways (a piece of this puzzle that I have yet to determine an efficient way to incorporate into the anvi'o annotations).

Now that you know how I generated the `KO_Orthology_ko00001.txt` lets get back to the fun part and import our functions!


### Parsing the results from GhostKOALA

Once GhostKOALA has finished running it will send you an email with a link to the results.

Download the annotation file, which will be called `user_ko.txt`, and make sure the results look like this:

``` bash
 $ head user_ko.txt
 genecall_0       K01923
 genecall_1
 genecall_2       K03611
 genecall_3
 genecall_4       K01952
 genecall_5       K01693
 genecall_6       K00817
 genecall_7       K00013
 genecall_8       K00765
(...)
```

Now run this command on your terminal to add the necessary header line to this file:

```
echo -e "contig\taccession_id" > .temp && cat user_ko.txt >> .temp && mv .temp user_ko.txt
```

Now your file is ready to be imported!

For this, we will use the program `KEGG-to-anvio` (which is in the directory you clone from GitHub early on).

The help menu is shown below:

```bash
 $ python KEGG-to-anvio -h
 
 usage: KEGG-to-anvio [-h] [--KeggDB KEGGDB] [-i I]
                        [--interproscan INTERPROSCAN] [-o O]

Combines annotation Data for input to anvio

optional arguments:
  -h, --help            show this help message and exit
  --KeggDB KEGGDB       identify the Kegg Orthology file (modified from htext
                        using given bash script)
  -i I                  specify the file containing GhostKoala Results
  --interproscan INTERPROSCAN
                        interproscan results
  -o O                  Specify an output file
```

Now the next step it pretty simple. If you don't have interproscan results, you can run it as:

```bash
$ python KEGG-to-anvio --KeggDB KO_Orthology_ko00001.txt \
                          -i user_ko.txt \
                          -o KeggAnnotations-AnviImportable.txt
```

if you have interproscan results, lets say in a output file called `interproscan-results.txt` (assuming interproscan was run with flags `-f tsv`, `--goterms`,  `--iprlookup`, and `--pathways`), then you can run this command as:

```bash
$ python KEGG-to-anvio --KeggDB KO_Orthology_ko00001.txt \
                          -i user_ko.txt \
                          -o KeggAnnotations-AnviImportable.txt \
                          --interproscan interproscan-results.txt
```

Now you should have a file called `KeggAnnotations-AnviImportable.txt` in your work directory, which can be imported into anvi'o using the program `anvi-import-functions`!! 

Here is how you can do it:

```bash
$ anvi-import-functions -c CONTIGS.db \
                        -i KeggAnnotations-AnviImportable.txt
```

<div class="extra-info" markdown="1">

<span class="extra-info-header">Two birds one stone</span>

If you decide you want to knock two birds out with one stone, you can also take the taxonomy data produced in GhostKOALA and convert that into an anvio importable format using the script `GhostKOALA-taxonomy-to-anvio` found in the GitHub repository you cloned.

The taxonomy file will download as a file called `user.out.top`. You can then run the parser like this:

``` bash
python GhostKOALA-taxonomy-to-anvio user.out.top KeggTaxonomy.txt
```

Which then can be imported into your anvi'o contigs database this way:

``` bash
anvi-import-taxonomy -c CONTIGS.db \
                     -i KeggTaxonomy.txt
```
</div>


You now now have KEGG functions in your anvi'o contigs database!

Hopefully (in the near future) we can find a better way to do this without us having to do all the convuluted steps with GhostKOALA.

{% include _join-anvio-slack.html %}

<div style="margin:50px">&nbsp;</div>

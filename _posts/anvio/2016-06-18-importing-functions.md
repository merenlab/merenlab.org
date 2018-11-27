---
layout: post
title: "Importing functions into contigs database"
excerpt: "Making those functions in the summary output bloom with stuff!"
modified: 2016-06-18
categories: [anvio]
comments: true
---

This article describes how to import functional annotations of genes into anvi'o.

Anvi'o also provides a simpler way to annotate your genes using NCBI COGs. If you would like to explore that option, please read [this article]({% post_url anvio/2016-10-25-cog-annotation %}).

{:.notice}
Have you been struggling with importing your functions into anvi'o? You want to use XXX to annotate your funcitons because everyone says it is awesome, but the output from XXX is not compatible with anvi'o? Let us know, and we will see what we can do!

Anvi'o accepts functional annotation at the gene-level. Gene IDs in your input files should correspond to gene IDs in the contigs database, hence, functional annotation should be done *after* exporting genes from the contigs database. The basic workflow goes like this: (1) generate your contigs database, (2) export your gene sequences, (3) assign functions to them, and (4) import results back into your contigs database.

{:.notice}
**Important note**: There are many ways to have your genes annotated with functions, but, there is only one way to make sure the gene IDs in the files you generate using external tools correspond to the gene IDs in the database: export your DNA or AA sequences from the anvi'o contigs database you wish to annotate using anvi'o program `anvi-get-sequences-for-gene-calls`.

{% include _toc.html %}

## Introduction

You will be using `anvi-import-functions` to do things described here, but before we start with the input file formats here is a little information about the program. You can import functions into anvi'o *incrementally*. For isntance, let's assume you have annotated your genes with PFAMs, imported those results into anvi'o. Then you decided to do it with TIGRFAMs, too. In fact you *can* import TIGRFAMs as an additional data into the database without erasing PFAMs. Functions can be imported multiple times from multiple sources without overwriting previous imports. This is the default behavior of the program, but alternatively, you can use the `--drop-previous-annotations` flag to make a fresh start, in which case all previous imports would be erased from the db, and only the final input will be stored in it.

When you import your gene calls you get more comprehensive outputs during summary. Also when you inspect your contigs, you get to click on gene calls and see how they are annotated:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/misc/gene_call.png"><img src="{{ site.url }}/images/anvio/misc/gene_call.png" width="50%" /></a>
</div>


## Simple matrix

This is the simplest way to get the functional annotation of genes into anvi'o. The TAB-delimited input matrix should follow this format:

|gene_callers_id|source|accession|function|e_value|
|:--|:--:|:--:|:--|:--:|
|1|Pfam|PF01132|Elongation factor P (EF-P) OB domain|4e-23|
|1|Pfam|PF08207|Elongation factor P (EF-P) KOW-like domain|3e-25|
|1|TIGRFAM|TIGR00038|efp: translation elongation factor P|1.5e-75|
|2|Pfam|PF01029|NusB family|2.5e-30|
|2|TIGRFAM|TIGR01951|nusB: transcription antitermination factor NusB|1.5e-36|
|3|Pfam|PF00117|Glutamine amidotransferase class-I|2e-36|
|3|Pfam|PF00988|Carbamoyl-phosphate synthase small chain, CPSase domain|1.2e-48|
|3|TIGRFAM|TIGR01368|CPSaseIIsmall: carbamoyl-phosphate synthase, small subunit|1.5e-132|
|4|Pfam|PF02787|Carbamoyl-phosphate synthetase large chain, oligomerisation domain|1.4e-31|
|4|TIGRFAM|TIGR01369|CPSaseII_lrg: carbamoyl-phosphate synthase, large subunit|0|
|5|TIGRFAM|TIGR02127|pyrF_sub2: orotidine 5'-phosphate decarboxylase|1.9e-59|
|6|Pfam|PF00625|Guanylate kinase|5.7e-39|
|6|TIGRFAM|TIGR03263|guanyl_kin: guanylate kinase|3.5e-62|
|8|Pfam|PF01192|RNA polymerase Rpb6|4.9e-13|
|8|TIGRFAM|TIGR00690|rpoZ: DNA-directed RNA polymerase, omega subunit|1.7e-20|
|9|TIGRFAM|TIGR01034|metK: methionine adenosyltransferase|2.5e-169|
|11|Pfam|PF13419|Haloacid dehalogenase-like hydrolase|2.8e-27|
|11|TIGRFAM|TIGR01509|HAD-SF-IA-v3: HAD hydrolase, family IA, variant 3|1.2e-11|
|12|Pfam|PF00551|Formyl transferase|1.4e-34|
|12|TIGRFAM|TIGR00460|fmt: methionyl-tRNA formyltransferase|2.9e-70|
|13|Pfam|PF12710|haloacid dehalogenase-like hydrolase|2.3e-14|
|13|TIGRFAM|TIGR00338|serB: phosphoserine phosphatase SerB|4.9e-76|
|13|TIGRFAM|TIGR01488|HAD-SF-IB: HAD phosphoserine phosphatase-like hydrolase, family IB|6e-29|
|14|Pfam|PF00004|ATPase family associated with various cellular activities (AAA)|7.7e-45|
|14|Pfam|PF16450|Proteasomal ATPase OB/ID domain|1.8e-34|
|14|TIGRFAM|TIGR03689|pup_AAA: proteasome ATPase|1e-206|
|(...)|(...)|(...)|(...)|(...)|

Once you have your matrix ready, this is the command line to import it:

``` bash
$ anvi-import-functions -c contigs.db -i input_matrix.txt
```

As you can see,

* Not every gene call has to be present in the matrix,
* It is OK if there are multiple annotations from the same source for a given gene call,
* It is OK if a give gene is annotated only by a single source.

If the **accession** information is not available to you, it is OK to leave it blank. If you have no e-values associated with your annotations, it is OK to put `0` for every entry. If there are multiple annotations from a single source for a single gene call, anvi'o uses e-values to use only the most significant one to show in interfaces.

## EggNOG database + emapper

Anvi'o has a parser for emapper output for [EggNOG database](http://eggnogdb.embl.de). After exporting your amino acid seqeunces from your contigs database,

``` bash
$ anvi-get-sequences-for-gene-calls -c CONTIGS.db \
                                    --get-aa-sequences \
                                    -o amino-acid-sequences.fa
```

You can either do it online, through the [interface on EggNOG](http://eggnogdb.embl.de/#/app/emapper), or locally by following the relevant [tutorial on emapper](https://github.com/jhcepas/eggnog-mapper).

Once you have your output file, say `amino-acid-sequences.fa.emapper.annotations`, you can import it back into the contigs database from which you exported the amino acid sequences:

``` bash
anvi-script-run-eggnog-mapper -c CONTIGS.db \
                              --annotation amino-acid-sequences.fa.emapper.annotations \
                              --use-version 1.0.3
```

{:.notice}
The script is called 'run' eggnog mapper, but what it does is to import its results. Well, historically, this script could run `emapper` on a given contigs database automatically and import results back into it seamlessly without leaving any room for user errors. However, due to the `SIGINT` and `SIGTERM` calls in the [emapper server code](https://github.com/jhcepas/eggnog-mapper/blob/master/eggnogmapper/server.py), the software is not quite software friendly and can't be called from within other software in a straightforward manner. Hence the `--annotation` flag as a workaround, so you can import your own results, and the weird name for the anvi'o script. Becasue bioinformatics.

The version in the command line above is the version number of `emapper` you used for annotation (because the number of fields change from one version to the other, and anvi'o needs to know exactly which version you used).


## InterProScan

Anvi'o has a parser for [InterProScan](http://www.ebi.ac.uk/interpro/download.html). To use InterProScan you should first export AA sequences for all your gene calls:

``` bash
$ anvi-get-sequences-for-gene-calls -c CONTIGS.db \
                                    --get-aa-sequences \
                                    -o amino-acid-sequences.fa
```

After running InterProScan on this file like this (assuming you have it [downloaded](http://www.ebi.ac.uk/interpro/download.html)):

``` bash
$ ./interproscan.sh -i amino-acid-sequences.fa \
                    -f tsv \
                    -o interpro-output.tsv
```

You can import results into the database:

``` bash
$ anvi-import-functions -c contigs.db \
                        -i interpro-output.tsv \
                        -p interproscan
```

That's it!


<div class="extra-info" markdown="1">

<span class="extra-info-header">Xabier Vázquez-Campos to the rescue</span>


### It doesn't work?

If the InterProScan parser within `anvi-import-functions` is not working for you, chances are that you used the `--iprlookup`, `--goterms` and/or the `--pathways` options when you ran InterProScan. This happens because the parser expects a table with 11 columns, and you have 15. Columns 12 and 13 are the InterPro integrated annotations, while 14 and 15 correspond to the GOterms and pathway annotations respectively.

Fortunately, there are ways around it. The easiest way it is to remove these extra columns and keep going:

```
cut -f 1-11 interpro-output.tsv > interpro-output_fix.tsv
anvi-import-functions -c contigs.db \
                        -i interpro-output_fix.tsv \
                        -p interproscan
```

Although this doesn't handle your _precious_ GOterm annotations or the InterPro cross-references (IPR numbers).

Alternatively, you can use the _standard import_. Although, you'll likely have redundant annotations (some protein domains can appear more than once in a single protein).

```
echo -e "gene_callers_id\tsource\taccession\tfunction\te_value" > table_to_import.tsv
cut -f 1,4-6,9 interpro-output.tsv >> table_to_import.tsv
anvi-import-functions -c contigs.db \
                        -i table_to_import.tsv
```

And, unfortunately, this may not work in some cases. E.g. if you got TMHMM or Phobius annotations, the e-value column will be always a dash `-`, not a number and the standard `anvi-import-functions` doesn't like it.

#### The final solution

To avoid any and all of these problems I created [`iprs2anvio.sh`](https://github.com/xvazquezc/stuff/blob/master/iprs2anvio.sh), a little BASH script that handles:

* InterPro cross-referenced annotations (the 12th and 13th columns) (`-r` or `--ipr` option).
* Dereplicates multiple matches of the same motif profile to a single protein, keeping the lowest e-value.
  + Also handles score-based annotations: ProSiteProfiles and HAMAP, and keeps the annotation with the higher score.
* Handles the dashes `-` in the e-value column.
* Allows to specify the annotation databases of interest and even split the output, i.e. a single file per database (`-s` or `--split`).
* Can extract optional annotations from columns 14 and 15:
  + GOterms (`-g` or `--go_terms`).
  + Pathways (`-p` or `--pathways`).

```
iprs2anvio.sh -i interpro-output.tsv -o interpro-output -g -p -r
anvi-import-functions -c contigs.db \
                      -i interpro-output_iprs2anvio.tsv
```

It is not the most efficient script in the world, but does the job.

{:.notice}
Note that `iprs2anvio.sh` has been developed on Linux and it hasn't been tested in Mac yet. It requires GNU [`getopt`](https://www.gnu.org/software/libc/manual/html_node/Getopt.html) to handle options, likely in your system already if you use Linux. For Mac, you can install it easily using [MacPorts](https://www.macports.org/install.php) with `sudo port install getopt`


[Xabier Vázquez-Campos](https://twitter.com/XabiVC).
</div>


## Prokka

You can use annotations produced by [Prokka](https://github.com/tseemann/prokka) with anvi'o, however, in order to do that you will need to create a new contigs database, since Prokka has to run on your contigs, and does not have a workaround to work only with gene calls. 

A recipe to work with Prokka by Antti Karkman is [available here]({% post_url anvio/2017-05-18-working-with-prokka %}).

## KEGG Modules

A tutorial on how to import KEGG modules into an anvi'o contigs database by Elaina Graham is [available here]({% post_url anvio/2018-01-17-importing-ghostkoala-annotations %}).

## A better way?

You have better / faster / more accurate ways to do it? Let us know!

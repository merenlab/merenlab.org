---
layout: page
title: anvi-db-info [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Access self tables, display values, or set new ones totally on your own risk.

See **[program help menu](../../../../vignette#anvi-db-info)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[pan-db](../../artifacts/pan-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[genomes-storage-db](../../artifacts/genomes-storage-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[structure-db](../../artifacts/structure-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[genes-db](../../artifacts/genes-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>

## Usage


Displays information about an anvi'o database, and allows users to modify that information when absolutely necessary.

This program is particularly useful for debugging, but also handy in a pinch if you want to check some facts about your database - to answer questions like "did I run HMMs on this <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> yet?" or "is this a merged <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span>?" This program can also be very dangerous when used to inappropriately modify database information, so if you want to change something, please proceed with caution.

### What information will I see?

All anvi'o databases contain a table of self-describing information known as the "self" table. It helps anvi'o keep track of critical facts such as the type of the database, its version number, and the date it was created. It also saves information about how the database was generated, what sorts of data it contains, what programs have been run on it, and so on. In general, this table exists so that anvi'o can make sure you are doing the right things with your data and that nothing will blow up. `anvi-db-info` will show you the contents of the self table when you run this program on an anvi'o database.

The information in the self table will be different depending on the kind of database you are looking at. For example, a <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> self table will indicate the number of contigs (and splits) in the database, whether or not gene calling was done (and with what gene callers), and which functional annotation sources have been used to annotate the genes. A <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span> self table will list which samples it contains mapping information for, how many reads where mapped from each sample, and whether or not SNVs have been profiled. A <span class="artifact-n">[modules-db](/software/anvio/help/7/artifacts/modules-db)</span> (see also <span class="artifact-n">[kegg-data](/software/anvio/help/7/artifacts/kegg-data)</span>) self table will tell you how many KEGG modules are saved in the database and what is the hash value of the database contents. We could go on, but you probably get the picture.

### View information about a database

This is the only way that most people will use this program, and it is very simple. Just provide the path to any anvi'o database to this program, and check the output on your terminal screen:

<div class="codeblock" markdown="1">
anvi&#45;db&#45;info path&#45;to&#45;DB.db
</div>

Let's be even more specific and say you have a <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> called `CONTIGS.db`. To look at its self table, you would run the following:
<div class="codeblock" markdown="1">
anvi&#45;db&#45;info CONTIGS.db
</div>

That's it! Easy-peasy lemon-squeezy.

### Example output

Here is an example of what you might see for a <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span>.

```
DB Info (no touch)
===============================================
Database Path ................................: CONTIGS.db
Description ..................................: [Not found, but it's OK]
Type .........................................: contigs
Version ......................................: 19


DB Info (no touch also)
===============================================
project_name .................................: TEST
contigs_db_hash ..............................: hashe451b636
split_length .................................: 20000
kmer_size ....................................: 4
num_contigs ..................................: 1
total_length .................................: 2089645
num_splits ...................................: 104
gene_level_taxonomy_source ...................: None
genes_are_called .............................: 1
splits_consider_gene_calls ...................: 1
creation_date ................................: 1584050598.57349
scg_taxonomy_was_run .........................: 1
external_gene_calls ..........................: 0
external_gene_amino_acid_seqs ................: 0
skip_predict_frame ...........................: 0
scg_taxonomy_database_version ................: v89
trna_taxonomy_was_run ........................: 0
trna_taxonomy_database_version ...............: None
modules_db_hash ..............................: 72700e4db2bc
gene_function_sources ........................: KEGG_Module,KOfam,Pfam,Transfer_RNAs,KEGG_Class

* Please remember that it is never a good idea to change these values. But in some
cases it may be absolutely necessary to update something here, and a programmer
may ask you to run this program and do it. But even then, you should be
extremely careful.

AVAILABLE GENE CALLERS
===============================================
* 'prodigal' (1,678 gene calls)
* 'Ribosomal_RNAs' (10 gene calls)


AVAILABLE HMM SOURCES
===============================================
* 'Bacteria_71' (type: singlecopy; num genes: 71)
* 'Archaea_76' (type: singlecopy; num genes: 76)
* 'Protista_83' (type: singlecopy; num genes: 83)
* 'Ribosomal_RNAs' (type: Ribosomal_RNAs; num genes: 12)

```

Most of this output is self-explanatory. But one thing that may not be quite obvious to some is that in many cases we use `0` to indicate 'False' and `1` to indicate 'True'. So for this example, you will see that SCG taxonomy was run on this database, but tRNA taxonomy was not.

### Modifying database information
We just need to start by saying - you probably shouldn't do this. Manually changing the values in the self table has the potential to break things downstream because it lets you avoid some of anvi'o's internal sanity checks which prevent you from doing things you shouldn't. If you change things and start running into ugly errors, do not be surprised.

That being said, sometimes you just need to live on the edge and do some hacking, and `anvi-db-info` will let you do that. If a programmer sent you here to update a value in the self table or if you are just foraging ahead on your own, this is how you would do it. Let's change the `project_name` value as an example because it is mostly descriptive and seems fairly safe:

<div class="codeblock" markdown="1">
anvi&#45;db&#45;info &#45;&#45;self&#45;key project_name &#45;&#45;self&#45;value "test" CONTIGS.db
</div>

If you run this, you will see a warning telling you what the current value of `project_name` is and what it will be changed to, but the value will not actually be changed just yet. If you are sure you want to do this, you then need to run:

<div class="codeblock" markdown="1">
anvi&#45;db&#45;info &#45;&#45;self&#45;key project_name &#45;&#45;self&#45;value "test" CONTIGS.db  &#45;&#45;just&#45;do&#45;it
</div>

Then go on your merry adventuring way.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-db-info.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-db-info) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

---
layout: post
title: "Importing taxonomy into contigs database"
excerpt: "Various ways to add the taxonomic annotations into anvi'o"
modified: 2016-06-18
categories: [anvio]
comments: true
---

This article describes various ways to import *gene-level* taxonomy into anvi'o (i.e., you have taxonomical annotations for each gene, and you want them to be in your database). Which can be very useful during the manual curation of metagenome-assembled genomes.

This post describes ways to import gene-level taxonomy data. We have some parsers to make sense of the output of some programs, but you can also import raw taxonomy for your genes you curated. If what you find here does not solve your problem, please feel free to suggest new ways to deal with gene-level taxonomy, [by entering an issue](http://github.com/meren/anvio).

The basic workflow goes like this: (1) generate your contigs database from your FASTA file, (2) export your gene sequences, (3) annotate them with taxonomy, and (4) import results back into your contigs database using `anvi-import-taxonomy` program.

{:.notice}
**Important note**: There are many ways to have your genes annotated with taxonomy. But, there is **only one way** to make sure the gene IDs in your taxonomy files correspond to the gene caller IDs in the database: export your DNA or AA sequences from the anvi'o contigs database you wish to annotate using the anvi'o program `anvi-get-sequences-for-gene-calls`. 

{% include _toc.html %}


## Kaiju

{:.notice}
[Jarrod J. Scott](https://scholar.google.com/citations?user=QaUurTgAAAAJ&hl=en),  kindly helped us to develop a parser for Kaiju and add it in anvi'o, and kindly provided the following tutorial to demonstrate how to import Kaiju output into anvi'o. We thank him very much for his patience and help throughout this entire process.

Anvi'o has a parser for the `Kaiju` classifier. The approach is described in in Menzel, P. et al. (2016) Fast and sensitive taxonomic classification for metagenomics with Kaiju. Nat. Commun. 7:11257 (open access). (http://www.nature.com/ncomms/2016/160413/ncomms11257/full/ncomms11257.html)

See the `Kaiju` GitHub page for detailed instructions on installing `Kaiju` as well as obtaining and formatting databases. (https://github.com/bioinformatics-centre/kaiju)

Assuming you have an `anvi'o` `contigs` database with gene calls, getting the file(s) for `Kaiju` is easy. Simply run:

`anvi-get-sequences-for-gene-calls -c CONTIGS.db -o gene_calls.fa`

Sweet, you are now ready to run `Kaiju` on your `fasta` file!

The `Kaiju` package comes with many handy scripts for analyzing the taxonomic content of you samples. For the purpose of the `anvi’o` parser and getting the `kaiju` output into you `contigs` database you will need to run three scripts:

First you need to run `makeDB.sh` to create a reference database.
Next, run the `kaiju` command to classify you genes.
Finally run `addTaxonNames` to amend the `kaiju` output with taxon names.  

First, download and format a reference database. You only need to run this once unless  you want multiple databases  or you want to update an existing database. `Kaiju` offers **four** main dbs—`NCBI RefSeq`, `proGenomes`, `NCBI nr`, and `Marine Metagenomic Portal`. See the `Kaiju` documentation or simply type `makeDB.sh -h` for a list of available dbs.

`makeDB.sh -e  -t 20`

This will download and format the NCBI BLAST non-redundant protein database “nr” with the addition of fungi and microbial eukaryotes using 20 parallel threads. This is the largest db offered by Kaiju and once  downloaded and formatted, is `~65GB`. By contrast the proGenome db is `13GB` and the Marine Portal is `~8GB`. See the `Kaiju` github page for detailed computational requirements for each database. 

Chose a suitable directory to house your database. Once database downloading and formatting is complete, `Kaiju` only needs the files `kaiju_db_nr_euk.fmi`, `nodes.dmp`, and `names.dmp` files. You can erase everything else including the `genome` directory, the `*.bwt` file, and the `*.sa` file.  

Again. you only need to run `makeDB.sh` once. 

Now its time to do some classifying. `Kaiju` offers many options to tweak your classification parameters. You can explore these further by typing `kaiju -h` or checking out the documentation. For the purposes of the `anvi’o` parser however you must, must, must enable the verbose output with the `-v` flag. The parser  expects a defined number of columns which are obtained with the `-v` flag   

To  run the classifier, `kaiju` needs the location of the `nodes.dmp` file, the database file (`.fmi`), your `fasta` file, and an output filename (PLUS the `-v` option), like so :

`kaiju -t /path/to/nodes.dmp -f /path/to/kaiju_db.fmi -i gene_calls.fa -o gene_calls_nr.out -z 16 -v`

The `-z` option indicates the number of parallel threads. For reference,~250,000 genes took about  3 minutes with 16 threads and 5GB of RAM (80GB total) plus another ~15 minutes to read and load the db. 

The output is a 7 column tab-delimited file  that contain **1)** whether the read is classified (C) or unclassified (U), **2**) the read ID, **3**) NCBI taxon ID of closest hit, **4**) best match score, **5**) comma-separated list of taxon IDs with best match, **6**) accession numbers of all database sequences with the best match, and **7**) matching amino acid fragment sequence. 

Cool? Now its time to add taxon names so that the `anvi’o` parser can add the taxonomy to the `contigs` database. The `addTaxonNames` script needs location of the `nodes.dmp` file, location of the `names.dmp` file, the output from the `kaiju` classification step, and an output name. Most importantly you **must** specify the taxonomic ranks to be included in the output as follows using the `-r` flag. **Comma-separated, NO SPACES**. 

`-r superkingdom,phylum,order,class,family,genus,species`

So the full command will look something like this:

`addTaxonNames -t path/to/nodes.dmp -n path/to/names.dmp -i  gene_calls_nr.out -o gene_calls_nr.names -r superkingdom,phylum,order,class,family,genus,species`

Whew. Now you are ready to run the `anvi’o` parser for `kaiju`.

At this point its not a bad idea to make a copy of your `contigs.db`---just in case. In order to get the `kaiju` taxonomic profile into your `contigs.db` you will run the following with the `-p kaiju` flag (indicating the name of the parser) and the input file, which is the output from the `addTaxonNames` step.

`anvi-import-taxonomy-for-genes -p kaiju -i gene_calls_nr.names`

`anvi’o` will immediately throw the following error--not because `anvi’o` is mad at you-- but because `anvi’o` is concerned about your well-being and doesn’t want you to ruin a perfectly good `contigs.db`. `Anvi'o` wants to make certain you are giving it a properly formatted input file :

`Config Error: Anvi'o assumes you used this exact parameter during your kaiju run: '-r        
              superkingdom,phylum,order,class,family,genus,species'. If you haven't, you will
              run into trouble later. If you are positive that you did include that parameter
              to your run, re-run this program with `--just-do-it` flag`.

So lets rerun with the `--just-do-it`` flag                     

This time anvi’o 
anvi-import-taxonomy-for-genes -p kaiju -i gene_calls_nr.names --just-do-it`

This time `anvi'o` will present a message indicating that it scanned your input file and gives you a second chance to back out of the arrangement. If the `phylum` level names `anvi'o` presents don't look correct you can `CNTRL C` out of the whole thing--otherwise, the taxonomy will be added to your `contigs.db`

Thats it! you can make sure everything looks ok by running:

 `anvi-export-table CONTIGS.db --table genes_taxonomy -l`
and hopefully you get something like this:

`Database .....................................: "CONTIGS.db" has been initiated with its 21 tables.`
`Table ........................................: "genes_taxonomy" has been read with 214509 entries and 2 columns.`
`Table columns ................................: "gene_callers_id, taxon_id"`


## Centrifuge

Anvi'o also has a parser for [centrifuge](https://github.com/infphilo/centrifuge). Please see [this recipe]({% post_url anvio/2016-06-18-installing-third-party-software %}#centrifuge) to download and set up centrifuge. 

Assuming you generated an anvi'o contigs database. To import taxonomy into this contigs database, first you will export all gene calls: 

{% highlight bash %}
$ anvi-get-sequences-for-gene-calls -c CONTIGS.db -o gene-calls.fa
{% endhighlight %}

Then you will run the following command:

{% highlight bash %}
$ centrifuge -f -x $CENTRIFUGE_BASE/p+h+v/p+h+v gene-calls.fa -S centrifuge_hits.tsv
{% endhighlight %}

{:.notice}
If the environment variable `$CENTRIFUGE_BASE` is not properly set, you will get an error. See the export instructions [here]({% post_url anvio/2016-06-18-installing-third-party-software %}#centrifuge) to try again.

This step takes about one minute on my laptop for 40,000 genes.

When centrifuge is done running, you should find two files in your work directory, which you will import into anvi'o. These files are `centrifuge_report.tsv`, and `centrifuge_hits.tsv`. Just to make sure that they are not empty, feel free to run this command:

{% highlight bash %}
$ wc -l centrifuge_report.tsv centrifuge_hits.tsv
    4666 centrifuge_report.tsv
  215625 centrifuge_hits.tsv
  220291 total
{% endhighlight %}

Fine. It is time to import these results! To do this, you will use the program `anvi-import-taxonomy` with the parser for `centrifuge`:


{% highlight bash %}
$ anvi-import-taxonomy -c CONTIGS.db -i centrifuge_report.tsv centrifuge_hits.tsv -p centrifuge
{% endhighlight %}

{:.notice}
You can use any file name you like, however, the order of input files is important: following the parameter `-i`, you should first declare the `report` file, and then the `hits` file.

This is it. If everything went alright, the interactive interface and anvi'o summary results should contain taxonomy information.


## Simple matrix

If you have taxonomy information for your genes, but none of the parsers listed here helps you to import them, the following is the simplest way to get the taxonomical annotation of genes into a contigs database.

Basically first you can create The TAB-delimited input matrix that follows this format:

|gene_callers_id|t_domain|t_phylum|t_class|t_order|t_family|t_genus|t_species|
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
|1|||||||Bacteroides fragilis|
|2|||||||Bacteroides fragilis|
|3|||||||Bifidobacterium longum|
|5|||||||Bifidobacterium longum|
|7|||||||Bifidobacterium longum|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|(...)|


And then you can use the following command to import it using the parser `default_matrix`:

{% highlight bash %}
$ anvi-import-taxonomy -c CONTIGS.db -i input_matrix.txt -p default_matrix
{% endhighlight %}

{:notice}
Not every gene call has to be in the matrix, and not every level of taxonomy has to be present, anvi'o will find a way to deal with that, but the more the merrier.

---


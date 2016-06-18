---
layout: post
title: "Importing taxonomy into contigs database"
excerpt: "Various ways to add the taxonomic annotations into anvi'o"
modified: 2016-06-18
categories: [anvio]
comments: true
---

This article describes various ways to import taxonomy into anvi'o.

Please feel free to suggest new ones, [by entering an issue](http://github.com/meren/anvio), and we will be happy to try to do our best to write a parser for it.

Anvi'o accepts taxonomic annotation at the gene-level. Annotations in your files should correspond to open reading frames in your contigs, hence, taxonomical annotation should be done *after* exporting DNA or AA sequences from the contigs database to maintain the link. The basic workflow goes like this: (1) generate your contigs database, (2) export your gene sequences and annotate them with taxonomy, and (3) import results back into the same contigs database. This article describes the third step.

{:.notice}
**Important note**: There are many ways to have your genes annotated with taxonomy. But, there is only one way to make sure the gene IDs in your taxonomy files correspond to the gene caller IDs in the database: export your DNA or AA sequences from the anvi'o contigs database you wish to annotate using anvi'o programs `anvi-get-dna-sequences-for-gene-calls` or `anvi-get-aa-sequences-for-gene-calls`. 

{% include _toc.html %}

## Simple matrix

This is the simplest way to get the taxonomical annotation of genes into the contigs database. The TAB-delimited input matrix should follow this format:

|gene_callers_id|t_phylum|t_class|t_order|t_family|t_genus|t_species|
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
|1||||||Bacteroides fragilis|
|2||||||Bacteroides fragilis|
|3||||||Bifidobacterium longum|
|5||||||Bifidobacterium longum|
|7||||||Bifidobacterium longum|
|(...)|(...)|(...)|(...)|(...)|(...)|(...)|

Not every gene call has to be in the matrix, and not every level of taxonomy has to be present.

Once you have your matrix ready, this is the command line to import it using the parser `default_matrix`:

{% highlight bash %}
$ anvi-import-taxonomy -c CONTIGS.db -i input_matrix.txt -p default_matrix
{% endhighlight %}

## Centrifuge Output

Anvi'o has a parser for [centrifuge](https://github.com/infphilo/centrifuge). Please see [this recipe]({% post_url anvio/2016-06-18-installing-third-party-software %}#centrifuge) to download and set up centrifuge. 

Assuming you generated an anvi'o contigs database. To import taxonomy into this contigs database, first you will export all gene calls: 

{% highlight bash %}
$ anvi-get-dna-sequences-for-gene-calls -c CONTIGS.db -o gene-calls.fa
{% endhighlight %}

Then you will run the following command:

{% highlight bash %}
$ centrifuge -f -x $CENTRIFUGE_BASE/b+h+v/b+h+v gene-calls.fa -S centrifuge_hits.tsv
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

This is it. If everything went alright, the interactive interface and anvi'o summary results should contain taxonomy information.

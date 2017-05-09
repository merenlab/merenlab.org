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
**Important note**: There are many ways to have your genes annotated with functions, but, there is only one way to make sure the gene IDs in the files you generate using external tools correspond to the gene IDs in the database: export your DNA or AA sequences from the anvi'o contigs database you wish to annotate using anvi'o programs `anvi-get-dna-sequences-for-gene-calls` or `anvi-get-aa-sequences-for-gene-calls`.

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

{% highlight bash %}
$ anvi-import-functions -c contigs.db -i input_matrix.txt
{% endhighlight %}

As you can see,

* Not every gene call has to be present in the matrix,
* It is OK if there are multiple annotations from the same source for a given gene call,
* It is OK if a give gene is annotated only by a single source.

If the **accession** information is not available to you, it is OK to leave it blank. If you have no e-values associated with your annotations, it is OK to put `0` for every entry. If there are multiple annotations from a single source for a single gene call, anvi'o uses e-values to use only the most significant one to show in interfaces.

## InterProScan

Anvi'o has a parser for [InterProScan](http://www.ebi.ac.uk/interpro/download.html). To use InterProScan you should first export AA sequences for all your gene calls:

{% highlight bash %}
$ anvi-get-aa-sequences-for-gene-calls -c CONTIGS.db -o protein-sequences.fa
{% endhighlight %}

After running InterProScan on this file like this (assuming you have it [downloaded](http://www.ebi.ac.uk/interpro/download.html)):

{% highlight bash %}
$ ./interproscan.sh -i protein-sequences.fa -f tsv -o interpro-output.tsv
{% endhighlight %}

You can import results into the database:

{% highlight bash %}
$ anvi-import-functions -c contigs.db -i interpro-output.tsv -p interproscan
{% endhighlight %}

That's it!

You have better / faster / more accurate ways to do it? Let us know!

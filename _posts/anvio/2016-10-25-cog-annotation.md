---
layout: post
title: "Annotating an anvi'o contigs database with COGs"
excerpt: "Yes. Good ol' COGs. Into your contigs db. Just like that. 60% of the time, every time."
modified: 2016-10-25
categories: [anvio]
comments: true
authors: [meren]
---

{: .notice}
The COG workflow is for `2.1.0` and later versions of anvi'o. You can learn which version you have on your computer by typing `anvi-profile --version` in your terminal.

This article describes how to setup COGs on your system, and how to annotate your gene calls in an anvi'o contigs database with COGs.

This is different than [importing functions into anvi'o]({{ post_url anvio/2016-06-18-importing-functions }}), because this workflow is run by anvi'o programs that use blastp or DIAMOND to search NCBI's now-quite-outdated-and-not-maintained-but-still-awesome COG database, and provides a one-step-solution for the functional annotation problem.

**Citation** information! If you are using this workflow, along with the search algorithm you will elect to use, you should also cite the following work: [Expanded microbial genome coverage and improved protein family annotation in the COG database](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4383993/).

{% include _toc.html %}

{:.notice}
{% include _fixthispage.html source="anvio/2016-11-25-cog-annotation.md" %}

## Annotating an anvi'o contigs database with COGs

It is quite simple. All you need is the program `anvi-run-ncbi-cogs`. Here is an example:

{% highlight bash %}
$ anvi-run-ncbi-cogs -c CONTIGS.db --num-threads 20
Directory to store temporary files ...........: /tmp/tmpVtIKcc
Directory will be removed after the run ......: True
Sequences ....................................: 51 sequences reported.
FASTA ........................................: /tmp/tmpVtIKcc/aa_sequences.fa
BLASTP results ...............................: /tmp/tmpVtIKcc/blast-search-results.txt

WARNING
===============================================
It seems there have already been functions in the db from 11 sources. Fine.
Anvio will populate this database INCREMENTALLY with the additional 94 entries
originating from 2 sources in the incoming data: COG_CATEGORY, COG_FUNCTION

Gene functions ...............................: 94 function calls from 2 sources for 47 unique gene calls has been added to the contigs database.
{% endhighlight %}

{:.notice}
Just to give a ballpark idea: annotations for about 60,000 amino acid sequences were added into the contigs database in about 5 hours using 50 cores.

Well, if this is your first time, then probably it will not go as smooth as it here as you will need to let anvi'o set up the COG distribution on your system.


## Setting up the COG distribution

This is something you will do only once (unless you have to do it again later for various reasons). During this step anvi'o will download necessary files from [ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/](ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/), and *reformat* them so they can be used later by `anvi-run-ncbi-cogs` program and other tools. The formatting step includes changing reorganizing information in raw files, serializing a very large text file into binary Python object for fast access while converting protein IDs to COGs, and finally generating BLAST and DIAMOND search databases.

Luckily, all you need to do is to run `anvi-setup-ncbi-cogs`. Here is an example:

{:.notice}
Depending on your installation, you may have to run the following command with superuser privileges! If you want to avoid that, you can use the `--cogs-data-dir` parameter to set up a different path.

{% highlight bash %}
$ anvi-setup-ncbi-cogs --num-threads 20
COG data dir .................................: /groups/merenlab/anvio/anvio/data/misc/COG

WARNING
===============================================
This program will first check whether you have all the raw files, and then will
attempt to regenerate everything that is necessary from them.

Downloaded succesfully .......................: /groups/merenlab/anvio/anvio/data/misc/COG/RAW_DATA_FROM_NCBI/cog2003-2014.csv
Downloaded succesfully .......................: /groups/merenlab/anvio/anvio/data/misc/COG/RAW_DATA_FROM_NCBI/prot2003-2014.fa.gz
Downloaded succesfully .......................: /groups/merenlab/anvio/anvio/data/misc/COG/RAW_DATA_FROM_NCBI/cognames2003-2014.tab
Downloaded succesfully .......................: /groups/merenlab/anvio/anvio/data/misc/COG/RAW_DATA_FROM_NCBI/fun2003-2014.tab
Diamond log ..................................: /groups/merenlab/anvio/anvio/data/misc/COG/DB_DIAMOND/log.txt
Diamond search db ............................: /groups/merenlab/anvio/anvio/data/misc/COG/DB_DIAMOND/COG.dmnd
BLAST log ....................................: /groups/merenlab/anvio/anvio/data/misc/COG/DB_BLAST/log.txt
BLAST search db ..............................: /groups/merenlab/anvio/anvio/data/misc/COG/DB_BLAST/COG
{% endhighlight %}

{:.notice}
If you are on a server system and you don't want the interface to ask you any questions, add `--just-do-it` flag to your command line.

{:.notice}
If something is wrong with your setup and you can't figure out what is the problem, try running `anvi-setup-ncbi-cogs` with the flag `--reset`.

If you are still reading this page it probably is there is a problem for which you couldn't find an answer here. Sigh. OK. How about sending an e-mail? :(

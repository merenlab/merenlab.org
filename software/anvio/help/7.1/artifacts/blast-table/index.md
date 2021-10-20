---
layout: page
title: blast-table [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/blast-table
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact is typically provided **by the user** for anvi'o to import into its databases, process, and/or use.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


There are no anvi'o tools that generate this artifact, which means it is most likely provided to the anvi'o ecosystem by the user.


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-script-filter-fasta-by-blast](../../programs/anvi-script-filter-fasta-by-blast)</span></p>


## Description

This describes the BLAST table that is outputted when you run [Protein BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins) from the terminal. 

When given to <span class="artifact-n">[anvi-script-filter-fasta-by-blast](/software/anvio/help/7.1/programs/anvi-script-filter-fasta-by-blast)</span>, which is currently the only program that uses this artifact, it expects output form 6. By default, this incldues the following data columns: 

    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen
    
However, you'll have to provide the columns in your file and their order to the program wirth the flag `--outfmt`. For the program to work properly, your table must at least include the columns `qseqid`, `bitscore`, `length`, `qlen`, and `pident`.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/blast-table.md) to update this information.


---
layout: page
title: anvi-export-locus [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

This program helps you cut a &#x27;locus&#x27; from a larger genetic context (e.g., contigs, genomes). By default, anvi&#x27;o will locate a user-defined anchor gene, extend its selection upstream and downstream based on the --num-genes argument, then extract the locus to create a new contigs database. The anchor gene must be provided as --search-term, --gene-caller-ids, or --hmm-sources. If --flank-mode is designated, you MUST provide TWO flanking genes that define the locus region (Please see --flank-mode help for more information). If everything goes as plan, anvi&#x27;o will give you individual locus contigs databases for every matching anchor gene found in the original contigs database provided. Enjoy your mini contigs databases!.

See **[program help menu](../../../../vignette#anvi-export-locus)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[locus-fasta](../../artifacts/locus-fasta) <img src="../../images/icons/FASTA.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>

## Usage


This program lets you export selections of your <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> around all occurances of a user-defined anchor gene. 

The output of this is a folder that contains a separate <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> for the region around each hit of the anchor gene. (In fact, you'll get a FASTA file, <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span>, <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span>, and a copy of the runlog).

For example, you could specify the recognition site for a specific enzyme and use this program to pull out all potential sites where that enzyme could bind. 

### Required Parameters

You'll need to provide a <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> (of course), as well as the name of the output directory and a prefix to use when naming all of the output databases. 

You can define the region of interest either by defining the two flanking genes or by searching for an anchor gene and defining a number of genes around this gene that you want to look at. For example, if you set `num-genes` as 1, then each locus will contain the gene of interest, a gene upstream of it, and a gene downstream of it, for a total of three genes. 

### Defining the region of interest

There are four ways to indicate the desired anchor gene:

1. Provide a search term in the functional annotations of all of your genes. (If you're trying to find a gene with a vague function, you might want to use <span class="artifact-n">[anvi-search-functions](/software/anvio/help/7/programs/anvi-search-functions)</span> to find out which genes will show up first. Alternatively, you can you <span class="artifact-n">[anvi-export-functions](/software/anvio/help/7/programs/anvi-export-functions)</span> to look at a full list of the functional annotaitons in this database). 

    <div class="codeblock" markdown="1">
    anvi&#45;export&#45;locus &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> \
                      &#45;&#45;num&#45;genes 2 \
                      &#45;o GLYCO_DIRECTORY \
                      &#45;O Glyco \
                      &#45;&#45;search&#45;term "Glycosyltransferase involved in cell wall bisynthesis" \ 
    </div>
    
    You also have the option to specify an annotation source with the flag `--annotation source`

2.  Provide a specific gene caller ID. 

    <div class="codeblock" markdown="1">
    anvi&#45;export&#45;locus &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> \
                      &#45;&#45;num&#45;genes 2 \
                      &#45;o output_directory \
                      &#45;O GENE_1 \
                      &#45;&#45;gene&#45;caller&#45;ids 1
    </div>

3. Provide a search term for the HMM source annotations. To do this, you must also specify an hmm-source. (You can use the flag `--list-hmm-sources` to list the available sources). 

    <div class="codeblock" markdown="1">
    anvi&#45;export&#45;locus &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> \
                      &#45;&#45;num&#45;genes 2 \
                      &#45;o Ribosomal_S20p \
                      &#45;O Ribosomal_S20p \
                      &#45;&#45;use&#45;hmm \
                      &#45;&#45;hmm&#45;source Bacteria_71 \
                      &#45;&#45;search&#45;term Ribosomal_S20p
    </div>
    
    4. Run in `flank-mode` and provide two flanking genes that define the locus region.
    
    <div class="codeblock" markdown="1">
    anvi&#45;export&#45;locus &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> \
                      &#45;&#45;flank&#45;mode \
                      &#45;o locus_output \
                      &#45;O gyclo_to_acyl \
                      &#45;&#45;search&#45;term "Glycosyltransferase involved in cell wall bisynthesis","Acyl carrier protein" \ 
    </div>

### Additional Options 

You can also remove partial hits, ignore reverse complement hits, or overwrite all files in a pre-existing output. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-export-locus.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-export-locus) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

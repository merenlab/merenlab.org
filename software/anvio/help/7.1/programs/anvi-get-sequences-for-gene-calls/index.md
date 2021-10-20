---
layout: page
title: anvi-get-sequences-for-gene-calls [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-get-sequences-for-gene-calls
image:
  featurerelative: ../../../images/header.png
  display: true
---

A script to get back sequences for gene calls.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[genomes-storage-db](../../artifacts/genomes-storage-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[genes-fasta](../../artifacts/genes-fasta) <img src="../../images/icons/FASTA.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[external-gene-calls](../../artifacts/external-gene-calls) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Usage


This program allows you to **export the sequences of your gene calls** from a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> or <span class="artifact-n">[genomes-storage-db](/software/anvio/help/7.1/artifacts/genomes-storage-db)</span> in the form of a <span class="artifact-n">[genes-fasta](/software/anvio/help/7.1/artifacts/genes-fasta)</span>. 

If you want other information about your gene calls from a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>, you can run <span class="artifact-n">[anvi-export-gene-calls](/software/anvio/help/7.1/programs/anvi-export-gene-calls)</span> (which outputs a <span class="artifact-n">[gene-calls-txt](/software/anvio/help/7.1/artifacts/gene-calls-txt)</span>) or get the coverage and detection information with <span class="artifact-n">[anvi-export-gene-coverage-and-detection](/software/anvio/help/7.1/programs/anvi-export-gene-coverage-and-detection)</span>.

### Running on a contigs database

You can run this program on a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> like so:

<div class="codeblock" markdown="1">
anvi&#45;get&#45;sequences&#45;for&#45;gene&#45;calls &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                                  &#45;o path/to/output
</div>

This is create a <span class="artifact-n">[genes-fasta](/software/anvio/help/7.1/artifacts/genes-fasta)</span> that contains every gene in your contigs database. If you only want a specific subset of genes, you can run the following: 

<div class="codeblock" markdown="1">
anvi&#45;get&#45;sequences&#45;for&#45;gene&#45;calls &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                                  &#45;o path/to/output \
                                  &#45;&#45;gene&#45;caller&#45;ids 897,898,1312 \
                                  &#45;&#45;delimiter ,
</div>

Now the resulting <span class="artifact-n">[genes-fasta](/software/anvio/help/7.1/artifacts/genes-fasta)</span> will contain only those three genes. 

You also have the option to report the output in [gff3 format](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md), report extended deflines for each gene, or report amino acid sequences instead of nucleotide sequences.

### Running on a genomes storage database

You can also get the sequences from gene calls in a <span class="artifact-n">[genomes-storage-db](/software/anvio/help/7.1/artifacts/genomes-storage-db)</span>, like so:

<div class="codeblock" markdown="1">
anvi&#45;get&#45;sequences&#45;for&#45;gene&#45;calls &#45;g <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/7.1/artifacts/genomes&#45;storage&#45;db)</span> \
                                  &#45;o path/to/output
</div>

This will create a <span class="artifact-n">[genes-fasta](/software/anvio/help/7.1/artifacts/genes-fasta)</span> that contains every gene in your genomes storage database. To focus on only a subset of the genomes contained in your database, use the flag `--genome-names`. You can provide a comma-delimited list of genome names or a flat text file that contains one genome per line. Alternatively, you could provide a list of gene-caller-ids as specified above. 

You also have the option to report the output in [gff3 format](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md), report extended deflines for each gene, or report amino acid sequences instead of nucleotide sequences.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-get-sequences-for-gene-calls.md) to update this information.


## Additional Resources


* [A tutorial on getting gene-level taxonomy for a contigs-db](http://merenlab.org/2016/06/18/importing-taxonomy/)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-get-sequences-for-gene-calls) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

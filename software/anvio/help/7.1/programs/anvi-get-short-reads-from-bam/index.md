---
layout: page
title: anvi-get-short-reads-from-bam [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-get-short-reads-from-bam
image:
  featurerelative: ../../../images/header.png
  display: true
---

Get short reads back from a BAM file with options for compression, splitting of forward and reverse reads, etc.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/ekiefl.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Evan Kiefl</span><div class="page-author-social-box"><a href="mailto:kiefl.evan@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/evankiefl" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/ekiefl" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[bin](../../artifacts/bin) <img src="../../images/icons/BIN.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[bam-file](../../artifacts/bam-file) <img src="../../images/icons/BAM.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[short-reads-fasta](../../artifacts/short-reads-fasta) <img src="../../images/icons/FASTA.png" class="artifact-icon-mini" /></span></p>


## Usage


This script get the short reads (in the form of a <span class="artifact-n">[short-reads-fasta](/software/anvio/help/7.1/artifacts/short-reads-fasta)</span>) out of a <span class="artifact-n">[bam-file](/software/anvio/help/7.1/artifacts/bam-file)</span>.  

A basic run of this program is as follows: 

<div class="codeblock" markdown="1">
anvi&#45;get&#45;short&#45;reads&#45;from&#45;bam &#45;o path/to/output \ 
                              BAM_FILE_1.bam BAM_FILE_2.bam
</div>

This will get all of the short reads out of the provided bam files (`BAM_FILE_1.bam` and `BAM_FILE_2.bam`) and put them into a single file. 

### Narrowing the input 

You can choose to only return the short reads that are contained within a <span class="artifact-n">[collection](/software/anvio/help/7.1/artifacts/collection)</span> or <span class="artifact-n">[bin](/software/anvio/help/7.1/artifacts/bin)</span>, as so:

<div class="codeblock" markdown="1">
anvi&#45;get&#45;short&#45;reads&#45;from&#45;bam &#45;o path/to/output \ 
                              &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                              &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7.1/artifacts/profile&#45;db)</span> \
                              &#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/7.1/artifacts/collection)</span> \
                              BAM_FILE_1.bam BAM_FILE_2.bam
</div>

### Changing the output format

You can split the output based on the directionality of paired-end reads. Adding the tag `--split-R1-and-R2` causes the program to create three separate output files: one for R1 (sequences in the forward direction), one for R2 (sequences in the reverse direction; i.e. reverse complement of R1 sequences), and one for unparied reads. When doing this, you can name these three files with a prefix by using the flag `-O`.  

<div class="codeblock" markdown="1">
anvi&#45;get&#45;short&#45;reads&#45;from&#45;bam &#45;o path/to/output \ 
                              &#45;&#45;split&#45;R1&#45;and&#45;R2 \ 
                              &#45;O BAM_1_and_BAM_2 \
                              BAM_FILE_1.bam BAM_FILE_2.bam
</div>

You can also compress the output by adding the flag `--gzip-output`


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-get-short-reads-from-bam.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-get-short-reads-from-bam) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

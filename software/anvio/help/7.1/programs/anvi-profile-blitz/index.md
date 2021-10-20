---
layout: page
title: anvi-profile-blitz [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-profile-blitz
image:
  featurerelative: ../../../images/header.png
  display: true
---

FAST profiling of BAM files to get contig- or gene-level coverage and detection stats. Unlike `anvi-profile`, which is another anvi&#x27;o program that can profile BAM files, this program is designed to be very quick and only report long-format files for various read recruitment statistics per item. Plase also see the program `anvi-script-get-coverage-from-bam` for recovery of data from BAM files without an anvi&#x27;o contigs database.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[bam-file](../../artifacts/bam-file) <img src="../../images/icons/BAM.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[bam-stats-txt](../../artifacts/bam-stats-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Usage


This program **produces a <span class="artifact-n">[bam-stats-txt](/software/anvio/help/7.1/artifacts/bam-stats-txt)</span> from one or more <span class="artifact-n">[bam-file](/software/anvio/help/7.1/artifacts/bam-file)</span> given a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>**. It is designed to serve people who only need to process read recruitment data stored in a <span class="artifact-n">[bam-file](/software/anvio/help/7.1/artifacts/bam-file)</span> to recover coverage and detection statistics (along with others) for their genes and/or contigs, and will report what's going on nicely with memory usage information and estimated time of completion:

[![anvi-profile-blitz](../../images/anvi-profile-blitz.png){:.center-img}](../../images/anvi-profile-blitz.png)

There are other programs in anvi'o software ecosystem that are similar to this one:

* <span class="artifact-n">[anvi-profile](/software/anvio/help/7.1/programs/anvi-profile)</span> also takes a <span class="artifact-n">[bam-file](/software/anvio/help/7.1/artifacts/bam-file)</span> and profiles it. **They both require a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>**. But while <span class="artifact-n">[anvi-profile](/software/anvio/help/7.1/programs/anvi-profile)</span> produces a <span class="artifact-n">[single-profile-db](/software/anvio/help/7.1/artifacts/single-profile-db)</span> for downstream analyses in anvi'o, <span class="artifact-n">[anvi-profile-blitz](/software/anvio/help/7.1/programs/anvi-profile-blitz)</span> produces text files for downstream analyses by the user (via R, Python, or other solutions). In contrast to <span class="artifact-n">[anvi-profile](/software/anvio/help/7.1/programs/anvi-profile)</span>, <span class="artifact-n">[anvi-profile-blitz](/software/anvio/help/7.1/programs/anvi-profile-blitz)</span> is orders of magnitude faster with similar memory usage.

* <span class="artifact-n">[anvi-script-get-coverage-from-bam](/software/anvio/help/7.1/programs/anvi-script-get-coverage-from-bam)</span> also takes a <span class="artifact-n">[bam-file](/software/anvio/help/7.1/artifacts/bam-file)</span> and profiles it. **They both produce text output files.** But while <span class="artifact-n">[anvi-script-get-coverage-from-bam](/software/anvio/help/7.1/programs/anvi-script-get-coverage-from-bam)</span> does not require a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>, <span class="artifact-n">[anvi-profile-blitz](/software/anvio/help/7.1/programs/anvi-profile-blitz)</span> requires one to work. They will both run very rapidly, <span class="artifact-n">[anvi-script-get-coverage-from-bam](/software/anvio/help/7.1/programs/anvi-script-get-coverage-from-bam)</span> will work with much smaller amount of memory.

## Output files

For output file formats, please see <span class="artifact-n">[bam-stats-txt](/software/anvio/help/7.1/artifacts/bam-stats-txt)</span>.

## Running

You can use this program with one or more BAM files to recover minimal or extended statistics for contigs or genes in a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>.

{:.warning}
Since the program will not be able to ensure the <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> was generated from the same <span class="artifact-n">[contigs-fasta](/software/anvio/help/7.1/artifacts/contigs-fasta)</span> that was used for read recruitment that resulted in <span class="artifact-n">[bam-file](/software/anvio/help/7.1/artifacts/bam-file)</span>s for analysis, you can make serious mistakes unless you mix up your workflow and start profiling BAM files that have nothing to do with a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>. If you make a mistake like that, in the best case scenario you will get an empty output file because the program will skip all contigs with non-matching name. In the worst case scenario you will get a file if some names in <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> incorrectly matches to some names in the <span class="artifact-n">[bam-file](/software/anvio/help/7.1/artifacts/bam-file)</span>. While this warning may be confusing, you can avoid all these if you use the SAME FASTA FILE both as reference for read recruitment and as input for <span class="artifact-n">[anvi-gen-contigs-database](/software/anvio/help/7.1/programs/anvi-gen-contigs-database)</span>.

### Contigs mode, default output

Profile contigs, produce a default output:

<div class="codeblock" markdown="1">
anvi&#45;profile&#45;blitz <span class="artifact&#45;n">[bam&#45;file](/software/anvio/help/7.1/artifacts/bam&#45;file)</span> \
                   &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                   &#45;o OUTPUT.txt
</div>

This example is with a single BAM file, but you can also have multiple BAM files as a parameter by using wildcards,

<div class="codeblock" markdown="1">
anvi&#45;profile&#45;blitz &#42;.bam \
                   &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                   &#45;o OUTPUT.txt
</div>

or by providing multiple paths:

<div class="codeblock" markdown="1">
anvi&#45;profile&#45;blitz /path/to/SAMPLE&#45;01.bam \
                   /path/to/SAMPLE&#45;02.bam \
                   /another/path/to/SAMPLE&#45;03.bam
                   &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                   &#45;o OUTPUT.txt
</div>

### Contigs mode, minimal output

Profile contigs, produce a minimal output. This is the fastest option:

<div class="codeblock" markdown="1">
anvi&#45;profile&#45;blitz <span class="artifact&#45;n">[bam&#45;file](/software/anvio/help/7.1/artifacts/bam&#45;file)</span> \
                   &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                   &#45;&#45;report&#45;minimal \
                   &#45;o OUTPUT.txt
</div>

### Genes mode, default output

Profile genes, produce a default output:

<div class="codeblock" markdown="1">
anvi&#45;profile&#45;blitz <span class="artifact&#45;n">[bam&#45;file](/software/anvio/help/7.1/artifacts/bam&#45;file)</span> \
                   &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                   &#45;&#45;gene&#45;mode \
                   &#45;o OUTPUT.txt
</div>

### Genes mode, minimal output

Profile genes, produce a default output:

<div class="codeblock" markdown="1">
anvi&#45;profile&#45;blitz <span class="artifact&#45;n">[bam&#45;file](/software/anvio/help/7.1/artifacts/bam&#45;file)</span> \
                   &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                   &#45;&#45;gene&#45;mode \
                   &#45;&#45;report&#45;minimal \
                   &#45;o OUTPUT.txt
</div>


## Performance

The memory use will be correlated linaerly with the size of the <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>, but once everything is loaded, the memory usage will not increase substantially over time.

With the flag `--report-minimal`, <span class="artifact-n">[anvi-profile-blitz](/software/anvio/help/7.1/programs/anvi-profile-blitz)</span> profiled on a laptop computer 100,000 contigs that contained 1 billion nts in 6 minutes and used  ~300 Mb memory. This contigs database had 1.5 million genes, and memory usage increased to 1.7 Gb when <span class="artifact-n">[anvi-profile-blitz](/software/anvio/help/7.1/programs/anvi-profile-blitz)</span> run in `--gene-mode`. The flag `--gene-mode` does not change time complexity dramatically.

Anvi'o has this program because [Emile Faure](https://twitter.com/faureemile) presented us with a [challenge](https://anvio.slack.com/archives/C8SFMGYF3/p1631723790065300): Emile had a ~140 Gb anvi'o <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> that contained nearly 70 million contig sequences from over 200 single-assembled metagenomes, and wanted to learn the coverages of each gene in the contigs database in 200 metagenomes individually. Yet the combination of <span class="artifact-n">[anvi-profile](/software/anvio/help/7.1/programs/anvi-profile)</span> and <span class="artifact-n">[anvi-summarize](/software/anvio/help/7.1/programs/anvi-summarize)</span> jobs would take **more than 40 days** to complete. Since all Emile needed was to learn the coverages from BAM files, we implemented <span class="artifact-n">[anvi-profile-blitz](/software/anvio/help/7.1/programs/anvi-profile-blitz)</span> to skip the profiling step. The run took **8 hours to compute and report coverage values for 175 million genes in 70 million contigs**, and the memory use remained below 200 Gb.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-profile-blitz.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-profile-blitz) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

---
layout: page
title: anvi-script-get-primer-matches [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-script-get-primer-matches
image:
  featurerelative: ../../../images/header.png
  display: true
---

You provide this program with FASTQ files for one or more samples AND one or more short sequences, and it collects reads from FASTQ files that matches to your sequences. This tool can be most powerful if you want to collect all short reads from one or more metagenomes that are downstream to a known sequence. Using the comprehensive output files you can analyze the diversity of seuqences visually, manually, or using established strategies such as oligotyping..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[samples-txt](../../artifacts/samples-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[short-reads-fasta](../../artifacts/short-reads-fasta) <img src="../../images/icons/FASTA.png" class="artifact-icon-mini" /></span></p>


## Usage


This program finds all reads in a given set of FASTQ files based on user-provided primer sequences.

The primary utility of this program is to get back short reads that may be extending into hypervariable regions of genomes that often suffer from significant drops in coverage in conventional read-recruitment analyses, thus preventing any meaningful insights into coverage or variability patterns.

In these situations, one can identify downstream conserved sequences (typically 15 to 25 nucleotides long) using the anvi'o interactive interface or through other means, and then provide those sequences to this program so it can find all matching sequences in a set of FASTQ files without any mapping.

{:.notice}
To instead get short reads mapping to a gene, use <span class="artifact-n">[anvi-get-short-reads-mapping-to-a-gene](/software/anvio/help/7.1/programs/anvi-get-short-reads-mapping-to-a-gene)</span>.

Here is a typical command line to run it:

<div class="codeblock" markdown="1">
anvi&#45;script&#45;get&#45;primer&#45;matches &#45;&#45;<span class="artifact&#45;n">[samples&#45;txt](/software/anvio/help/7.1/artifacts/samples&#45;txt)</span> samples.txt  \
                               &#45;&#45;primer&#45;sequences sequences.txt \
                               &#45;&#45;output&#45;dir OUTPUT
</div>

The <span class="artifact-n">[samples-txt](/software/anvio/help/7.1/artifacts/samples-txt)</span> file is to list all the samples one is interested in, and the primer sequences file lists each primer sequence of interest. Each of these files can contain a single entry, or multiple ones.

This will output all of the matching sequences into three <span class="artifact-n">[fasta](/software/anvio/help/7.1/artifacts/fasta)</span> files in the directory `OUTPUT`. These <span class="artifact-n">[fasta](/software/anvio/help/7.1/artifacts/fasta)</span> files differ in their format and will include those that describe,

* Raw sequences: sequences from the FASTQ files that matched to a primer where each sequence reported as is with no processing.
* Trimmed sequences: Raw sequences where the upstream of the primer sequence trimmed, as a result all matching sequences will start at the same position, and
* Gapped sequences: Trimmed sequences padded with gap characters to eliminate length variation artificially.

The last two formats provide downstream possibilities to generate <span class="artifact-n">[oligotypes](/software/anvio/help/7.1/artifacts/oligotypes)</span> and cluster short reads from an hypervariable region to estimate their diversity and oligotype proportion. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-script-get-primer-matches.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-get-primer-matches) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

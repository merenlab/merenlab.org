---
layout: page
title: anvi-script-get-coverage-from-bam [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Get nucleotide-level, contig-level, or bin-level coverage values from a BAM file.

See **[program help menu](../../../../vignette#anvi-script-get-coverage-from-bam)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[coverages-txt](../../artifacts/coverages-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[bam-file](../../artifacts/bam-file) <img src="../../images/icons/BAM.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[collection-txt](../../artifacts/collection-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>

## Usage


This program gets the coverage values from a <span class="artifact-n">[bam-file](/software/anvio/help/7/artifacts/bam-file)</span>, and puts them into a <span class="artifact-n">[coverages-txt](/software/anvio/help/7/artifacts/coverages-txt)</span>. 

You must provide a BAM file, but there are three ways you can choose contigs to analyze within that file: 
1. Give a contig name. Here, you can only report coverage per nucleotide position (In this example, the user is specifically asking for this anyway with the `-m` flag)

    <div class="codeblock" markdown="1">
    anvi&#45;script&#45;get&#45;coverage&#45;from&#45;bam &#45;b <span class="artifact&#45;n">[bam&#45;file](/software/anvio/help/7/artifacts/bam&#45;file)</span> \ 
                                     &#45;c NAME_OF_CONTIG \ 
                                     &#45;m pos
    </div>

2. Give a file that contains a list of contigs (one per line; same format as the `--contigs-of-interest` tag for <span class="artifact-n">[anvi-profile](/software/anvio/help/7/programs/anvi-profile)</span>). Here, you can ask for the contig averages (as in this example) or nucleotide position coverage. 

    <div class="codeblock" markdown="1">
    anvi&#45;script&#45;get&#45;coverage&#45;from&#45;bam &#45;b <span class="artifact&#45;n">[bam&#45;file](/software/anvio/help/7/artifacts/bam&#45;file)</span> \ 
                                     &#45;l NAME_OF_FILE \
                                     &#45;m contig
    </div>

3. Give a <span class="artifact-n">[collection-txt](/software/anvio/help/7/artifacts/collection-txt)</span> file for the program to determine the coverage for all contigs in those bins. Here, you can ask for the contig averages, nucleotide position coverage or coverage per bin (as in this example). 

    <div class="codeblock" markdown="1">
    anvi&#45;script&#45;get&#45;coverage&#45;from&#45;bam &#45;b <span class="artifact&#45;n">[bam&#45;file](/software/anvio/help/7/artifacts/bam&#45;file)</span> \ 
                                     &#45;C <span class="artifact&#45;n">[collection&#45;txt](/software/anvio/help/7/artifacts/collection&#45;txt)</span> \
                                     &#45;m bin
    </div>



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-script-get-coverage-from-bam.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-get-coverage-from-bam) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

---
layout: page
title: anvi-script-get-coverage-from-bam [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-script-get-coverage-from-bam
image:
  featurerelative: ../../../images/header.png
  display: true
---

Get nucleotide-level, contig-level, or bin-level coverage values from a BAM file very rapidly. For other anvi&#x27;o programs that are designed to profile BAM files, see `anvi-profile` and `anvi-profile-blitz`.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/ekiefl.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Evan Kiefl</span><div class="page-author-social-box"><a href="mailto:kiefl.evan@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/evankiefl" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/ekiefl" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[bam-file](../../artifacts/bam-file) <img src="../../images/icons/BAM.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[collection-txt](../../artifacts/collection-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[coverages-txt](../../artifacts/coverages-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Usage


This program gets the coverage values from a <span class="artifact-n">[bam-file](/software/anvio/help/7.1/artifacts/bam-file)</span>, and puts them into a <span class="artifact-n">[coverages-txt](/software/anvio/help/7.1/artifacts/coverages-txt)</span>. 

You must provide a BAM file, but there are three ways you can choose contigs to analyze within that file: 
1. Give a contig name. Here, you can only report coverage per nucleotide position (In this example, the user is specifically asking for this anyway with the `-m` flag)

    <div class="codeblock" markdown="1">
    anvi&#45;script&#45;get&#45;coverage&#45;from&#45;bam &#45;b <span class="artifact&#45;n">[bam&#45;file](/software/anvio/help/7.1/artifacts/bam&#45;file)</span> \ 
                                     &#45;c NAME_OF_CONTIG \ 
                                     &#45;m pos
    </div>

2. Give a file that contains a list of contigs (one per line; same format as the `--contigs-of-interest` tag for <span class="artifact-n">[anvi-profile](/software/anvio/help/7.1/programs/anvi-profile)</span>). Here, you can ask for the contig averages (as in this example) or nucleotide position coverage. 

    <div class="codeblock" markdown="1">
    anvi&#45;script&#45;get&#45;coverage&#45;from&#45;bam &#45;b <span class="artifact&#45;n">[bam&#45;file](/software/anvio/help/7.1/artifacts/bam&#45;file)</span> \ 
                                     &#45;l NAME_OF_FILE \
                                     &#45;m contig
    </div>

3. Give a <span class="artifact-n">[collection-txt](/software/anvio/help/7.1/artifacts/collection-txt)</span> file for the program to determine the coverage for all contigs in those bins. Here, you can ask for the contig averages, nucleotide position coverage or coverage per bin (as in this example). 

    <div class="codeblock" markdown="1">
    anvi&#45;script&#45;get&#45;coverage&#45;from&#45;bam &#45;b <span class="artifact&#45;n">[bam&#45;file](/software/anvio/help/7.1/artifacts/bam&#45;file)</span> \ 
                                     &#45;C <span class="artifact&#45;n">[collection&#45;txt](/software/anvio/help/7.1/artifacts/collection&#45;txt)</span> \
                                     &#45;m bin
    </div>



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-script-get-coverage-from-bam.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-get-coverage-from-bam) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

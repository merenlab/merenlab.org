---
layout: page
title: anvi-script-reformat-fasta [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-script-reformat-fasta
image:
  featurerelative: ../../../images/header.png
  display: true
---

Reformat FASTA file (remove contigs based on length, or based on a given list of deflines, and/or generate an output with simpler names).

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[fasta](../../artifacts/fasta) <img src="../../images/icons/FASTA.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[contigs-fasta](../../artifacts/contigs-fasta) <img src="../../images/icons/FASTA.png" class="artifact-icon-mini" /></span></p>


## Usage


This program **converts a <span class="artifact-n">[fasta](/software/anvio/help/7.1/artifacts/fasta)</span> file to a <span class="artifact-n">[contigs-fasta](/software/anvio/help/7.1/artifacts/contigs-fasta)</span>.** In other words, it reformats your FASTA formatted file to meet the conditions required of a <span class="artifact-n">[contigs-fasta](/software/anvio/help/7.1/artifacts/contigs-fasta)</span>, which is able to be used by other anvi'o programs.

<div class="codeblock" markdown="1">
anvi&#45;script&#45;reformat&#45;fasta <span class="artifact&#45;n">[fasta](/software/anvio/help/7.1/artifacts/fasta)</span> \
                           &#45;o <span class="artifact&#45;n">[contigs&#45;fasta](/software/anvio/help/7.1/artifacts/contigs&#45;fasta)</span> \
                           &#45;&#45;simplify&#45;names
</div>

{:.notice}
If you use the flag *--report-file*, it will also create a TAB-delimited file for you to keep track of which defline in the new file corresponds to which defline in the original file.

### Removing the short reads

Removing short contigs from a FASTA file will improve the performance of the <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> later. The example below runs the same command while also removing sequences that are shorter than 1,000 nts:

<div class="codeblock" markdown="1">
anvi&#45;script&#45;reformat&#45;fasta <span class="artifact&#45;n">[fasta](/software/anvio/help/7.1/artifacts/fasta)</span> \
                           &#45;o <span class="artifact&#45;n">[contigs&#45;fasta](/software/anvio/help/7.1/artifacts/contigs&#45;fasta)</span> \
                           &#45;l 1000 \
                           &#45;&#45;simplify&#45;names
</div>



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-script-reformat-fasta.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-reformat-fasta) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

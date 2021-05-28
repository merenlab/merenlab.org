---
layout: page
title: anvi-export-contigs [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Export contigs (or splits) from an anvi&#x27;o contigs database.

See **[program help menu](../../../../vignette#anvi-export-contigs)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[contigs-fasta](../../artifacts/contigs-fasta) <img src="../../images/icons/FASTA.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>

## Usage


This program **exports the contig sequences from a <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span>**, outputting them as a <span class="artifact-n">[contigs-fasta](/software/anvio/help/7/artifacts/contigs-fasta)</span>. It also has the ability to output the sequences of your splits instead. 

You can run this program as follows: 

<div class="codeblock" markdown="1">
anvi&#45;export&#45;contigs &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> \
                    &#45;o path/to/<span class="artifact&#45;n">[contigs&#45;fasta](/software/anvio/help/7/artifacts/contigs&#45;fasta)</span>
</div>

To run it on only a named subset of your contigs, you can provide a list of contigs as a separate file (in the same format as a <span class="artifact-n">[splits-txt](/software/anvio/help/7/artifacts/splits-txt)</span>). For example: 

<div class="codeblock" markdown="1">
anvi&#45;export&#45;contigs &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> \
                    &#45;o path/to/<span class="artifact&#45;n">[contigs&#45;fasta](/software/anvio/help/7/artifacts/contigs&#45;fasta)</span> \
                    &#45;&#45;contigs&#45;of&#45;interest my_favorite_contigs.txt 
</div>

where `my_favorite_contigs.txt` looks like this:

    contig_0001
    contig_0005
    contig_0035
    
### Splits mode

Want to look at your splits instead of your contigs? Just run with the flag `splits-mode` attached. 

<div class="codeblock" markdown="1">
anvi&#45;export&#45;contigs &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> \
                    &#45;o path/to/<span class="artifact&#45;n">[contigs&#45;fasta](/software/anvio/help/7/artifacts/contigs&#45;fasta)</span> \
                    &#45;&#45;splits&#45;mode
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-export-contigs.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-export-contigs) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

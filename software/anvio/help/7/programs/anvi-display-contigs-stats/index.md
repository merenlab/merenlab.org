---
layout: page
title: anvi-display-contigs-stats [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Start the anvi&#x27;o interactive interactive for viewing or comparing contigs statistics.

See **[program help menu](../../../../vignette#anvi-display-contigs-stats)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[contigs-stats](../../artifacts/contigs-stats) <img src="../../images/icons/STATS.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[interactive](../../artifacts/interactive) <img src="../../images/icons/DISPLAY.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[svg](../../artifacts/svg) <img src="../../images/icons/SVG.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>

## Usage


This program **lets you look at various stats in your <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span>(s)**. 

When you run with all default parameters, as so,

<div class="codeblock" markdown="1">
anvi&#45;display&#45;contigs&#45;stats <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span>
</div>

it will open the interactive interface. From there, you'll be able to see a breakdown of the genes in your contigs database. This will look something like this: 

![An example of the anvi'o interface for contigs stats](../../images/contigs-stats-interface-example.png)

Let's walk through the stats displayed here: 

At the top of the page are two graphs: 
* The bars in the top graph represent every integer N and L statistic from 1 to 100. The y-axis is the respective N length and the x-axis is the percentage of the total dataset looked at (the exact L and N values can be seen by hovering over each bar). In other words, if you had sorted your contigs by length (from longest to shortest), and walked through each one, every time you had seen another 1 percent of your total dataset, you would add a bar to the graph showing the number of contigs that you had seen (the L statistic) and the length of the one you were looking at at the moment (the N statistic). 
* The lower part of the graph tells you about which HMM hits your contigs database has. Each column is a gene in a specific <span class="artifact-n">[hmm-source](/software/anvio/help/7/artifacts/hmm-source)</span>, and the graph tells you how many hits each gene has in your data. (Hover your mouse over the graph to see the specifics of each gene.) The sidebar shows you how many of the genes in this graph were seen exactly that many times. For example, in the graph above, for the Bacteria_71 <span class="artifact-n">[hmm-source](/software/anvio/help/7/artifacts/hmm-source)</span>, a lot of genes were detected 9-11 times, so those bars are longer. This helps you estimate about how many of these genomes there are in your contigs database (so here, there is likely around 9-11 bacteria genomes in this contigs database). 

Below this are the **contigs stats** which are displayed in the following order:
- The total length of your contigs in nucleotides
- The number of contigs in your database
- The number of contigs that are of varying lengths. (for example "Num Contigs > 2.5 kb" gives you the number of contigs that are longer than 2500 base pairs)
- The length of the longest and shortest contig in your database in nucleotides
- The number of genes in your contigs (as predicted by [Prodigal](https://github.com/hyattpd/Prodigal))
- L50, L75, L90: If you ordered the contigs in your database from longest to shortest, these stats describe the *number of contigs* you would need to go through before you had looked at a certain percent of a genome. For example, L50 describes the number of contigs you would have to go through before you reached 50 percent of the entire dataset. 
- N50, N75, N90:  If you ordered the contigs in your database from longest to shortest, these stats describe the *length of the contig* you would be looking when you had looked at a certain percent of a genome. For example, N50 describes the length of contig you would be on when you reached 50 percent of the entire genome length. 
- The number of HMM hits in your contigs. This goes through every <span class="artifact-n">[hmm-source](/software/anvio/help/7/artifacts/hmm-source)</span> and gives the number of hits its genes had in all of your contigs. Basically, this is the number of hits that is given in the lower graph at the top of the page. 
- The number of genomes that anvi'o predicts are in your sample, based on how many hits the single copy core genes got from the various <span class="artifact-n">[hmm-source](/software/anvio/help/7/artifacts/hmm-source)</span>s. See the description of the lower graph above, or [this blog post](http://merenlab.org/2015/12/07/predicting-number-of-genomes/) for more information. 

This interface is espeically useful if you want to compare multiple databases, since you can view all of their stats stimultaneously.

You can also change various server configuration settings when you run this command. 

### The interactive interface is great, but I need text

Then you're still in the right place. Just add the tag `--report-as-text` and you'll get a lovely tab-deliminated output in the path provided:

<div class="codeblock" markdown="1">
anvi&#45;display&#45;contigs&#45;stats <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> \ 
                    &#45;&#45;report&#45;as&#45;text \
                    &#45;o path/to/my_cool_text_file.txt
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-display-contigs-stats.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-display-contigs-stats) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

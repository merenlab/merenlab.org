---
layout: page
title: anvi-gen-phylogenomic-tree [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Generate phylogenomic tree from aligment file.

See **[program help menu](../../../../vignette#anvi-gen-phylogenomic-tree)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[phylogeny](../../artifacts/phylogeny) <img src="../../images/icons/NEWICK.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[concatenated-gene-alignment-fasta](../../artifacts/concatenated-gene-alignment-fasta) <img src="../../images/icons/FASTA.png" class="artifact-icon-mini" /></span></p>

## Usage


This program generates a NEWICK-formatted phylogenomic tree (see <span class="artifact-n">[phylogeny](/software/anvio/help/7/artifacts/phylogeny)</span>) based on a given <span class="artifact-n">[concatenated-gene-alignment-fasta](/software/anvio/help/7/artifacts/concatenated-gene-alignment-fasta)</span>. 

As mentioned in the [phylogenetics tutorial](http://merenlab.org/2017/06/07/phylogenomics/), it currently only has the option to use [FastTree](http://microbesonline.org/fasttree/) to do so, but be aware that there are many other programs that you can do this with. Some of the options we are familiar with (and are not yet represented in `anvi-gen-phylogenomic-tree`) include [MrBayes](http://mrbayes.sourceforge.net/), [MEGA](http://www.megasoftware.net/), and PHYLIP, [among many others](http://evolution.genetics.washington.edu/phylip/software.html#methods), most of which will happily take a <span class="artifact-n">[concatenated-gene-alignment-fasta](/software/anvio/help/7/artifacts/concatenated-gene-alignment-fasta)</span>. 

Anyway, running this program is simple. Just provide the <span class="artifact-n">[concatenated-gene-alignment-fasta](/software/anvio/help/7/artifacts/concatenated-gene-alignment-fasta)</span> with all of the genes that you want to use and the output file path for your <span class="artifact-n">[phylogeny](/software/anvio/help/7/artifacts/phylogeny)</span>:

<div class="codeblock" markdown="1">
anvi&#45;gen&#45;phylogenomic&#45;tree &#45;f <span class="artifact&#45;n">[concatenated&#45;gene&#45;alignment&#45;fasta](/software/anvio/help/7/artifacts/concatenated&#45;gene&#45;alignment&#45;fasta)</span> \
                           &#45;o PATH/TO/<span class="artifact&#45;n">[phylogeny](/software/anvio/help/7/artifacts/phylogeny)</span>
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-gen-phylogenomic-tree.md) to update this information.


## Additional Resources


* [View this program in action in the anvi&#x27;o phylogenetics workflow](http://merenlab.org/2017/06/07/phylogenomics/)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-gen-phylogenomic-tree) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

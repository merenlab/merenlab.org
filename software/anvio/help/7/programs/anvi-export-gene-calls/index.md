---
layout: page
title: anvi-export-gene-calls [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Export gene calls from an anvi&#x27;o contigs database.

See **[program help menu](../../../../vignette#anvi-export-gene-calls)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[gene-calls-txt](../../artifacts/gene-calls-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>

## Usage


This program exports your gene calls given a <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> and a gene caller. The output of this is a <span class="artifact-n">[gene-calls-txt](/software/anvio/help/7/artifacts/gene-calls-txt)</span>. 

To see the gene callers available in your contigs database, run 

<div class="codeblock" markdown="1">
anvi&#45;export&#45;gene&#45;calls &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> \
                       &#45;&#45;list&#45;gene&#45;callers
</div>

By default, this list will probably include [Prodigal](https://github.com/hyattpd/Prodigal), which identifies genes when creating a <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span>. For in this example, we'll use export Prodigal-identified genes. Note that you can also get genes from more than one source by providing several gene-callers in a comma-delimited list.  

Then, you can export all of your gene callers in an orderly fashion by running 

<div class="codeblock" markdown="1">
anvi&#45;export&#45;gene&#45;calls &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> \
                       &#45;&#45;gene&#45;caller Prodigal \
                       &#45;o path/to/output
</div>

This will put a <span class="artifact-n">[gene-calls-txt](/software/anvio/help/7/artifacts/gene-calls-txt)</span> in the path you specified containing all of your Prodigal genes. 

If you don't want to display the amino acid sequences of each gene (they can crowd the file very quickly if you don't want to see them), just add the flag `--skip-sequence-reporting`


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-export-gene-calls.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-export-gene-calls) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

---
layout: page
title: anvi-matrix-to-newick [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Takes a distance matrix, returns a newick tree.

See **[program help menu](../../../vignette#anvi-matrix-to-newick)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[dendrogram](../../artifacts/dendrogram)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[view-data](../../artifacts/view-data)</span></p>

## Usage


This program converts a distance matrix (computed from a <span class="artifact-n">[view-data](/software/anvio/help/artifacts/view-data)</span> artifact) into a <span class="artifact-n">[dendrogram](/software/anvio/help/artifacts/dendrogram)</span>. 




{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-matrix-to-newick.md) to update this information.


## Additional Resources


* [See this program in action in the pangenomics tutorial](http://merenlab.org/2016/11/08/pangenomics-v2/#creating-a-quick-pangenome-with-functions)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-matrix-to-newick) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

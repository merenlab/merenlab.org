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

See **[program help menu](../../../../vignette#anvi-matrix-to-newick)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[dendrogram](../../artifacts/dendrogram) <img src="../../images/icons/NEWICK.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[view-data](../../artifacts/view-data) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>

## Usage


This program converts a distance matrix (computed from a <span class="artifact-n">[view-data](/software/anvio/help/7/artifacts/view-data)</span> artifact) into a <span class="artifact-n">[dendrogram](/software/anvio/help/7/artifacts/dendrogram)</span>. 

It uses the numerical data in a <span class="artifact-n">[view-data](/software/anvio/help/7/artifacts/view-data)</span> to compute a distance matrix behind the scenes, and then runs some hierarchical clustering to create a <span class="artifact-n">[dendrogram](/software/anvio/help/7/artifacts/dendrogram)</span> for all of your items. 

With all default parameters, a run would look like this:

<div class="codeblock" markdown="1">
anvi&#45;matrix&#45;to&#45;newick &#45;o path/for/<span class="artifact&#45;n">[dendrogram](/software/anvio/help/7/artifacts/dendrogram)</span> \ 
                      <span class="artifact&#45;n">[view&#45;data](/software/anvio/help/7/artifacts/view&#45;data)</span> 
</div>

If your input file has your samples as rows instead of columns, just add the flag `--transpose`. 

You can also ask for an additional output file: the order of the items in the resulting dendrogram as a <span class="artifact-n">[misc-data-items-order](/software/anvio/help/7/artifacts/misc-data-items-order)</span> in LIST format. To get this, simply provide a path to its desired location  with `--items-order-file`. 

Additionally, for hierarchical clustering, you can change the distance metric (a full list of the available metrics can be found [here](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html)) or the linkage method (though this is not recommended, the list of options can be found [here](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html)).


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-matrix-to-newick.md) to update this information.


## Additional Resources


* [See this program in action in the pangenomics tutorial](http://merenlab.org/2016/11/08/pangenomics-v2/#creating-a-quick-pangenome-with-functions)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-matrix-to-newick) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

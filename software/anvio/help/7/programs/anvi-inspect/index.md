---
layout: page
title: anvi-inspect [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Start an anvi&#x27;o inspect interactive interface.

See **[program help menu](../../../../vignette#anvi-inspect)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[interactive](../../artifacts/interactive) <img src="../../images/icons/DISPLAY.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[bin](../../artifacts/bin) <img src="../../images/icons/BIN.png" class="artifact-icon-mini" /></span></p>

## Usage


This lets you inspect a single split across your samples. This interface can also be opened from the <span class="artifact-n">[anvi-interactive](/software/anvio/help/7/programs/anvi-interactive)</span> interface by asking for details about a specific split.

From this view, you can clearly see the coverage and detection across your split, all SNVs, and the genes identified within your split and their functional annotations. You can also  easily compare all of this data across all of the samples that this split is present in.  

To run this program, just provide a <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span> and <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> pair and a single split name to inspect. 

<div class="codeblock" markdown="1">
anvi&#45;inspect &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7/artifacts/profile&#45;db)</span> \
             &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> \ 
             &#45;&#45;split&#45;name Day17a_QCcontig9_split_00003
</div>

You can also choose to hide SNVs marked as outliers or configure the server in various ways. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-inspect.md) to update this information.


## Additional Resources


* [Visualizing contig coverages](https://merenlab.org/2019/11/25/visualizing-coverages/)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-inspect) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

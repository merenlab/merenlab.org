---
layout: page
title: anvi-import-taxonomy-for-layers [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Import layers-level taxonomy into an anvi&#x27;o additional layer data table in an anvi&#x27;o single-profile database.

See **[program help menu](../../../../vignette#anvi-import-taxonomy-for-layers)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[layer-taxonomy](../../artifacts/layer-taxonomy) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[single-profile-db](../../artifacts/single-profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[layer-taxonomy-txt](../../artifacts/layer-taxonomy-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>

## Usage


This program lets you associate your layers with taxonomic information through a <span class="artifact-n">[single-profile-db](/software/anvio/help/7/artifacts/single-profile-db)</span>. 

This information is displayed in the interactive interface at the same place as <span class="artifact-n">[misc-data-layers](/software/anvio/help/7/artifacts/misc-data-layers)</span>, which is point (4) on [this page](http://merenlab.org/2017/12/11/additional-data-tables/#views-items-layers-orders-some-anvio-terminology). 

If instead you want the layers to *represent* taxonomic ranks, then you'll want to take a look at [this tutorial on phylogenomics](http://merenlab.org/2017/06/07/phylogenomics/).

Usually, the layers describe separate samples. However, when working with only one sample, you may break up different aspects of that sample to be represented in each layer, hence why you might want to associate them with taxonomic information. 

To run this program, simply provide a <span class="artifact-n">[layer-taxonomy-txt](/software/anvio/help/7/artifacts/layer-taxonomy-txt)</span>

<div class="codeblock" markdown="1">
anvi&#45;import&#45;taxonomy&#45;for&#45;layers &#45;p <span class="artifact&#45;n">[single&#45;profile&#45;db](/software/anvio/help/7/artifacts/single&#45;profile&#45;db)</span> \
                                &#45;i <span class="artifact&#45;n">[layer&#45;taxonomy&#45;txt](/software/anvio/help/7/artifacts/layer&#45;taxonomy&#45;txt)</span> 
</div>

You also have the option to change the minimum abundance cut off using `--min-abundance`. The default value is 0.1 percent. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-import-taxonomy-for-layers.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-import-taxonomy-for-layers) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

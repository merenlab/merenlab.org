---
layout: page
title: anvi-gen-variability-network [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

A program to generate a network description from an anvi&#x27;o variability profile (potentially outdated program).

See **[program help menu](../../../../vignette#anvi-gen-variability-network)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[variability-profile](../../artifacts/variability-profile) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>

## Usage


{:.notice}
**No one has described the usage of this program** :/ If you would like to contribute, please see previous examples [here](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs), and feel free to add a Markdown formatted file in that directory named "anvi-gen-variability-network.md". For a template, you can use the markdown file for `anvi-gen-contigs-database`. THANK YOU!


## Additional Resources


* [This program in the Infant Gut Tutorial](https://merenlab.org/tutorials/infant-gut/#visualizing-snv-profiles-as-a-network)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-gen-variability-network) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

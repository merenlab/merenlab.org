---
layout: page
title: anvi-export-gene-coverage-and-detection [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Export gene coverage and detection data for all genes associated with contigs described in a profile database.

See **[program help menu](../../../../vignette#anvi-export-gene-coverage-and-detection)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[coverages-txt](../../artifacts/coverages-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[detection-txt](../../artifacts/detection-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>

## Usage


This program gives you the **coverage and detection data** for all of the genes found in your <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span>, using the short reads data in your <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span>. 

<div class="codeblock" markdown="1">
anvi&#45;export&#45;gene&#45;coverage&#45;and&#45;detection &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> \
                                        &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7/artifacts/profile&#45;db)</span> \
                                        &#45;O MY_DATA
</div>

This will give you a <span class="artifact-n">[coverages-txt](/software/anvio/help/7/artifacts/coverages-txt)</span> and a <span class="artifact-n">[detection-txt](/software/anvio/help/7/artifacts/detection-txt)</span> whose file names will begin with `MY_DATA`


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-export-gene-coverage-and-detection.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-export-gene-coverage-and-detection) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

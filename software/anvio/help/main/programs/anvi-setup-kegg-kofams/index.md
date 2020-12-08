---
layout: page
title: anvi-setup-kegg-kofams [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Download and setup KEGG KOfam HMM profiles.

See **[program help menu](../../../../vignette#anvi-setup-kegg-kofams)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[kegg-data](../../artifacts/kegg-data) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[modules-db](../../artifacts/modules-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"></p>

## Usage


<span class="artifact-n">[anvi-setup-kegg-kofams](/software/anvio/help/main/programs/anvi-setup-kegg-kofams)</span> downloads and organizes data from KEGG for use by other programs, namely <span class="artifact-n">[anvi-run-kegg-kofams](/software/anvio/help/main/programs/anvi-run-kegg-kofams)</span> and <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/main/programs/anvi-estimate-metabolism)</span>. It downloads HMM profiles from the KOfams database as well as metabolism information such as that stored in the [KEGG MODULES resource](https://www.genome.jp/kegg/module.html). The KOfam profiles are prepared for later use by the HMMER software, and the metabolism information is put into a <span class="artifact-n">[modules-db](/software/anvio/help/main/artifacts/modules-db)</span>. The program generates a directory with these files (<span class="artifact-n">[kegg-data](/software/anvio/help/main/artifacts/kegg-data)</span>), which by default is located at `anvio/anvio/data/misc/KEGG/`.

### Set up KEGG data

<div class="codeblock" markdown="1">
anvi&#45;setup&#45;kegg&#45;kofams
</div>

### Set up KEGG data in non-default location

<div class="codeblock" markdown="1">
anvi&#45;setup&#45;kegg&#45;kofams &#45;&#45;kegg&#45;data&#45;dir /path/to/directory/KEGG
</div>

An important thing to note about this program is that it has rigid expectations for the format of the KEGG data that it works with. Future updates to KEGG may break things such that the data can no longer be directly obtained from KEGG or properly processed. In the event that this happens, this program still has you covered. You can provide an archived KEGG data directory to the script, and it will unpack that archive and make sure things are all in order.

### Set up from archived KEGG data

<div class="codeblock" markdown="1">
anvi&#45;setup&#45;kegg&#45;kofams &#45;&#45;kegg&#45;archive KEGG_archive.tar.gz
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-setup-kegg-kofams.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-setup-kegg-kofams) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

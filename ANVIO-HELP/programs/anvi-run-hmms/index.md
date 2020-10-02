---
layout: page
title: anvi-run-hmms [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

This program deals with populating tables that store HMM hits in an anvi&#x27;o contigs database..

See **[program help menu](../../../vignette#anvi-run-hmms)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[hmm-hits](../../artifacts/hmm-hits)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span> <span class="artifact-r">[hmm-source](../../artifacts/hmm-source)</span></p>

## Usage


Stores <span class="artifact-n">[hmm-hits](/software/anvio/help/artifacts/hmm-hits)</span> for a given <span class="artifact-n">[hmm-source](/software/anvio/help/artifacts/hmm-source)</span> in a <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span>.

One of the programs that users commonly run on newly generated <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span>, along with <span class="artifact-n">[anvi-scan-trnas](/software/anvio/help/programs/anvi-scan-trnas)</span>, <span class="artifact-n">[anvi-run-ncbi-cogs](/software/anvio/help/programs/anvi-run-ncbi-cogs)</span>, <span class="artifact-n">[anvi-run-scg-taxonomy](/software/anvio/help/programs/anvi-run-scg-taxonomy)</span>, and so on.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-run-hmms.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-run-hmms) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

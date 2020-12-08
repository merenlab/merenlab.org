---
layout: page
title: anvi-import-taxonomy-for-genes [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Import gene-level taxonomy into an anvi&#39;o contigs database.

See **[program help menu](../../../vignette#anvi-import-taxonomy-for-genes)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[gene-taxonomy](../../artifacts/gene-taxonomy)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span> <span class="artifact-r">[gene-taxonomy-txt](../../artifacts/gene-taxonomy-txt)</span></p>

## Usage


This program uses a <span class="artifact-n">[gene-taxonomy-txt](/software/anvio/help/main/artifacts/gene-taxonomy-txt)</span> to populate the taxonomic information for the genes in a <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span>. 

Once finished, your gene taxonomy will appear as an additional layer if you open the <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span> and an associated <span class="artifact-n">[profile-db](/software/anvio/help/main/artifacts/profile-db)</span> in <span class="artifact-n">[anvi-interactive](/software/anvio/help/main/programs/anvi-interactive)</span>. 

There is an entire blogpost about different ways to do this [here](http://merenlab.org/2016/06/18/importing-taxonomy/). It outlines how to get your sequences using <span class="artifact-n">[anvi-get-sequences-for-gene-calls](/software/anvio/help/main/programs/anvi-get-sequences-for-gene-calls)</span> than use either [Kaiju](https://github.com/bioinformatics-centre/kaiju) or [Centrifuge](https://github.com/infphilo/centrifuge) to get the taxonomy information for your genes. Finally, you bring that information back into anvi'o using this program.  


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-import-taxonomy-for-genes.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-import-taxonomy-for-genes) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

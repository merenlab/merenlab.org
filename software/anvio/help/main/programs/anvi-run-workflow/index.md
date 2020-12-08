---
layout: page
title: anvi-run-workflow [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Execute, manage, parallelize, and troubleshoot entire &#39;omics workflows and chain together anvi&#39;o and third party programs.

See **[program help menu](../../../vignette#anvi-run-workflow)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[contigs-workflow](../../artifacts/contigs-workflow)</span> <span class="artifact-p">[metagenomics-workflow](../../artifacts/metagenomics-workflow)</span> <span class="artifact-p">[pangenomics-workflow](../../artifacts/pangenomics-workflow)</span> <span class="artifact-p">[phylogenomics-workflow](../../artifacts/phylogenomics-workflow)</span> <span class="artifact-p">[trnaseq-workflow](../../artifacts/trnaseq-workflow)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[samples-txt](../../artifacts/samples-txt)</span> <span class="artifact-r">[fasta-txt](../../artifacts/fasta-txt)</span> <span class="artifact-r">[workflow-config](../../artifacts/workflow-config)</span></p>

## Usage


{:.notice}
**No one has described the usage of this program** :/ If you would like to contribute, please see previous examples [here](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs), and feel free to add a Markdown formatted file in that directory named "anvi-run-workflow.md". For a template, you can use the markdown file for `anvi-gen-contigs-database`. THANK YOU!


## Additional Resources


* [Tutorial](http://merenlab.org/2018/07/09/anvio-snakemake-workflows/)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-run-workflow) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

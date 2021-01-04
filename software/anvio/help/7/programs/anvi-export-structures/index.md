---
layout: page
title: anvi-export-structures [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Export .pdb structure files from a structure database.

See **[program help menu](../../../../vignette#anvi-export-structures)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[protein-structure-txt](../../artifacts/protein-structure-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[structure-db](../../artifacts/structure-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>

## Usage



This program exports the structures from a <span class="artifact-n">[structure-db](/software/anvio/help/7/artifacts/structure-db)</span> into the globally understood pdb format (<span class="artifact-n">[protein-structure-txt](/software/anvio/help/7/artifacts/protein-structure-txt)</span>), so they may be used for any follow-up analyses taking place outside of anvi'o.


To run, just provide a <span class="artifact-n">[structure-db](/software/anvio/help/7/artifacts/structure-db)</span> and an output path: 

<div class="codeblock" markdown="1">
anvi&#45;export&#45;structures &#45;s <span class="artifact&#45;n">[structure&#45;db](/software/anvio/help/7/artifacts/structure&#45;db)</span> \
                       &#45;o path/to/output
</div>

You can also provide a list of gene caller IDs, either directly through the parameter `--gene-caller-ids` or through a file with one gene caller ID per line.




{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-export-structures.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-export-structures) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

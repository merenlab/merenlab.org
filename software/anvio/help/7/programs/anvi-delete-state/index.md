---
layout: page
title: anvi-delete-state [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Delete an anvi&#x27;o state from a pan or profile database.

See **[program help menu](../../../../vignette#anvi-delete-state)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[pan-db](../../artifacts/pan-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[state](../../artifacts/state) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>

## Usage


This program, as implied by the name, is used to delete a <span class="artifact-n">[state](/software/anvio/help/7/artifacts/state)</span> from a <span class="artifact-n">[pan-db](/software/anvio/help/7/artifacts/pan-db)</span> or <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span>. This way, you can remove states that are clogging up the state list in the interface. 

It is generally a good idea to export your state before deleting it, just in case ((anvi-export-state)s).

To list available <span class="artifact-n">[state](/software/anvio/help/7/artifacts/state)</span>s in a database, call 

<div class="codeblock" markdown="1">
anvi&#45;delete&#45;state &#45;p <span class="artifact&#45;n">[pan&#45;db](/software/anvio/help/7/artifacts/pan&#45;db)</span> \
                 &#45;&#45;list&#45;states
</div>

Then, you can easily delete a <span class="artifact-n">[state](/software/anvio/help/7/artifacts/state)</span> with the command

<div class="codeblock" markdown="1">
anvi&#45;delete&#45;hmms &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7/artifacts/profile&#45;db)</span> \
                 &#45;s <span class="artifact&#45;n">[state](/software/anvio/help/7/artifacts/state)</span> 
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-delete-state.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-delete-state) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

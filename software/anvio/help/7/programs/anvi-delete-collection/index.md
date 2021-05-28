---
layout: page
title: anvi-delete-collection [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Remove a collection from a given profile database.

See **[program help menu](../../../../vignette#anvi-delete-collection)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[collection](../../artifacts/collection) <img src="../../images/icons/COLLECTION.png" class="artifact-icon-mini" /></span></p>

## Usage


This program, as implied by the name, is used to delete a <span class="artifact-n">[collection](/software/anvio/help/7/artifacts/collection)</span> from a <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span>. 

When you do this, you'll lose the collection forever, as well as the <span class="artifact-n">[bin](/software/anvio/help/7/artifacts/bin)</span>s within it. It is generally a good idea to export your binning effort into a <span class="artifact-n">[collection-txt](/software/anvio/help/7/artifacts/collection-txt)</span> using <span class="artifact-n">[anvi-export-collection](/software/anvio/help/7/programs/anvi-export-collection)</span> before deleting it, just to be safe. 

To list available collections in a database, call 

<div class="codeblock" markdown="1">
anvi&#45;delete&#45;collection &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7/artifacts/profile&#45;db)</span> \
                       &#45;&#45;list&#45;collections
</div>

Then, you can easily delete a collection with the command

<div class="codeblock" markdown="1">
anvi&#45;delete&#45;collection &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7/artifacts/profile&#45;db)</span> \
                       &#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/7/artifacts/collection)</span>
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-delete-collection.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-delete-collection) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

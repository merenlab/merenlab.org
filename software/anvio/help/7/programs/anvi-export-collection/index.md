---
layout: page
title: anvi-export-collection [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Export a collection from an anvi&#x27;o database.

See **[program help menu](../../../../vignette#anvi-export-collection)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[collection-txt](../../artifacts/collection-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[collection](../../artifacts/collection) <img src="../../images/icons/COLLECTION.png" class="artifact-icon-mini" /></span></p>

## Usage


This program, as one might think, allows you to export a <span class="artifact-n">[collection](/software/anvio/help/7/artifacts/collection)</span>. This allows you to take your binning results elsewhere (including into another Anvi'o project with the command <span class="artifact-n">[anvi-import-collection](/software/anvio/help/7/programs/anvi-import-collection)</span>). 

You can run this program on a <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span> or <span class="artifact-n">[pan-db](/software/anvio/help/7/artifacts/pan-db)</span> as follows: 

<div class="codeblock" markdown="1">
anvi&#45;export&#45;collection &#45;C my_favorite_collection \
                        &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7/artifacts/profile&#45;db)</span> 
</div>

This will give you a <span class="artifact-n">[collection-txt](/software/anvio/help/7/artifacts/collection-txt)</span> file that describes the collection `my_favorite_collection`. 

To list the collections available in this database, you can run 

<div class="codeblock" markdown="1">
anvi&#45;export&#45;collection &#45;p <span class="artifact&#45;n">[pan&#45;db](/software/anvio/help/7/artifacts/pan&#45;db)</span> \
                        &#45;&#45;list&#45;colllections
</div>

You can also add the flag `--include-unbinned` to have all unbinned contigs in the database show up at the end of your <span class="artifact-n">[collection-txt](/software/anvio/help/7/artifacts/collection-txt)</span> file in a bin titled `UNBINNED`. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-export-collection.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-export-collection) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

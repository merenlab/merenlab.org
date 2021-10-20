---
layout: page
title: anvi-export-functions [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-export-functions
image:
  featurerelative: ../../../images/header.png
  display: true
---

Export functions of genes from an anvi&#x27;o contigs database for a given annotation source.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[functions](../../artifacts/functions) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[functions-txt](../../artifacts/functions-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Usage


This program **takes in a <span class="artifact-n">[functions](/software/anvio/help/7.1/artifacts/functions)</span> artifact to create a <span class="artifact-n">[functions-txt](/software/anvio/help/7.1/artifacts/functions-txt)</span>.** Basically, if you want to take the information in your <span class="artifact-n">[functions](/software/anvio/help/7.1/artifacts/functions)</span> artifact out of anvi'o or give it to a fellow anvi'o user (for them to [import it](http://merenlab.org/software/anvio/help/programs/anvi-import-functions/) into their own project), you get that information using this command. 

Simply provide the <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> that has been annotated with <span class="artifact-n">[functions](/software/anvio/help/7.1/artifacts/functions)</span>: 

<div class="codeblock" markdown="1">
anvi&#45;export&#45;functions &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> 
</div>

You can also get annotations for only a specific list of sources. For example:

<div class="codeblock" markdown="1">
anvi&#45;export&#45;functions &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                      &#45;&#45;annotation&#45;sources source_1,source_2,source_3
</div>



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-export-functions.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-export-functions) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

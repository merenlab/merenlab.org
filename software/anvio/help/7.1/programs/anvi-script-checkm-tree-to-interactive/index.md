---
layout: page
title: anvi-script-checkm-tree-to-interactive [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-script-checkm-tree-to-interactive
image:
  featurerelative: ../../../images/header.png
  display: true
---

A helper script to convert CheckM trees into anvio interactive with taxonomy information.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/ozcan.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Özcan C. Esen</span><div class="page-author-social-box"><a href="http://blog.ozcanesen.com/" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:ozcanesen@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/ozcanesen" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/ozcan" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[phylogeny](../../artifacts/phylogeny) <img src="../../images/icons/NEWICK.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[interactive](../../artifacts/interactive) <img src="../../images/icons/DISPLAY.png" class="artifact-icon-mini" /></span></p>


## Usage


A helper script to process CheckM tree output to generate files compatible with <span class="artifact-n">[anvi-interactive](/software/anvio/help/7.1/programs/anvi-interactive)</span>.

An example use:

<div class="codeblock" markdown="1">
anvi&#45;script&#45;checkm&#45;tree&#45;to&#45;interactive &#45;t CheckM_concatenated.tree \
                                       &#45;o OUTPUT_PATH
cd OUTPUT_PATH/
anvi&#45;interactive &#45;p PROFILE.db \
                 &#45;t newick.tree \
                 &#45;d view_data.txt \
                 &#45;&#45;manual
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-script-checkm-tree-to-interactive.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-checkm-tree-to-interactive) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

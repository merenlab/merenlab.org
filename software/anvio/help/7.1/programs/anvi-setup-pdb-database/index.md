---
layout: page
title: anvi-setup-pdb-database [program]
categories: [anvio]
comments: false
redirect_from: /m/anvi-setup-pdb-database
image:
  featurerelative: ../../../images/header.png
  display: true
---

Setup or update an offline database of representative PDB structures clustered at 95%.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/ekiefl.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Evan Kiefl</span><div class="page-author-social-box"><a href="mailto:kiefl.evan@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/evankiefl" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/ekiefl" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


This program seems to know what its doing. It needs no input material from its user. Good program.


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[pdb-db](../../artifacts/pdb-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>


## Usage



## Basic usage 

This program creates a <span class="artifact-n">[pdb-db](/software/anvio/help/main/artifacts/pdb-db)</span> local database that holds PDB structures from [this sequence database](https://salilab.org/modeller/supplemental.html), which is hosted by the [Sali lab](https://salilab.org/).  Their database comprises all PDB RCSB sequences that have been clustered at 95% sequence similarity. They seem to update their database every couple of months (thank you guys!).


The purpose of <span class="artifact-n">[anvi-setup-pdb-database](/software/anvio/help/main/programs/anvi-setup-pdb-database)</span> to have a local copy of reference structures that can be used to, for example, get template structures for homology modelling when <span class="artifact-n">[anvi-gen-structure-database](/software/anvio/help/main/programs/anvi-gen-structure-database)</span> is ran.


Running this program is easy:

<div class="codeblock" markdown="1">
anvi&#45;setup&#45;pdb&#45;database &#45;&#45;just&#45;do&#45;it
</div>

If you already have a <span class="artifact-n">[pdb-db](/software/anvio/help/main/artifacts/pdb-db)</span> artifact and are trying to redownload this data, run 

<div class="codeblock" markdown="1">
anvi&#45;setup&#45;pdb&#45;database &#45;&#45;reset
</div>

Or if you just want to update your database, run 

<div class="codeblock" markdown="1">
anvi&#45;setup&#45;pdb&#45;database &#45;&#45;update
</div>

## Notes

The output <span class="artifact-n">[pdb-db](/software/anvio/help/main/artifacts/pdb-db)</span> database is ~20GB and its contents may take several hours to download.



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-setup-pdb-database.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-setup-pdb-database) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

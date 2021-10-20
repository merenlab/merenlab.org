---
layout: page
title: anvi-setup-interacdome [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-setup-interacdome
image:
  featurerelative: ../../../images/header.png
  display: true
---

Setup InteracDome data.

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


<p style="text-align: left" markdown="1"><span class="artifact-p">[interacdome-data](../../artifacts/interacdome-data) <img src="../../images/icons/DATA.png" class="artifact-icon-mini" /></span></p>


## Usage



This program (much like all of the other programs that begin with `anvi-setup`) sets up a local copy of the InteracDome database for <span class="artifact-n">[anvi-run-interacdome](/software/anvio/help/7.1/programs/anvi-run-interacdome)</span> as well as a local copy of Pfam v31.0, which is what InteracDome is defined for. Note that anvi'o only needs this program to be run once.


Specifically, this downloads [InteracDome](https://interacdome.princeton.edu/)’s [tab-separated files](https://interacdome.princeton.edu/#tab-6136-4) and the Pfam v31.0 HMM profiles for the Pfams in your InteracDome data. This data is stored in the <span class="artifact-n">[interacdome-data](/software/anvio/help/7.1/artifacts/interacdome-data)</span> artifact. 


It's easy as 1-2-3:

<div class="codeblock" markdown="1">
anvi&#45;setup&#45;interacdome
</div>

When running this program, you can provide a path to store your InteracDome data in. The default path is `anvio/data/misc/InteracDome`; if you use a custom path, you will have to provide it to <span class="artifact-n">[anvi-run-interacdome](/software/anvio/help/7.1/programs/anvi-run-interacdome)</span> with the same parameter. Here is an example run: 


<div class="codeblock" markdown="1">
anvi&#45;setup&#45;interacdome &#45;&#45;interacdome&#45;data&#45;dir path/to/directory 
</div>

If you want to overwrite any data that you have already downloaded (for example if you suspect something went wrong in the download), add the `--reset` flag: 

<div class="codeblock" markdown="1">
anvi&#45;setup&#45;interacdome  &#45;&#45;interacdome&#45;data&#45;dir path/to/directory \ 
                        &#45;&#45;reset
</div>



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-setup-interacdome.md) to update this information.


## Additional Resources


* [The setup step in the InteracDome technical blogpost](http://merenlab.org/2020/07/22/interacdome/#anvi-setup-interacdome)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-setup-interacdome) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

---
layout: page
title: anvi-delete-hmms [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-delete-hmms
image:
  featurerelative: ../../../images/header.png
  display: true
---

Remove HMM hits from an anvi&#x27;o contigs database.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[hmm-source](../../artifacts/hmm-source) <img src="../../images/icons/HMM.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[hmm-hits](../../artifacts/hmm-hits) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>


## Can provide


This program does not seem to provide any artifacts. Such programs usually print out some information for you to see or alter some anvi'o artifacts without producing any immediate outputs.


## Usage


This program, as implied by the name, is used to delete a <span class="artifact-n">[hmm-hits](/software/anvio/help/7.1/artifacts/hmm-hits)</span> from a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>. This way, you can repopulate the function annotations with a different source or program or just delete data that's clogging up the interface.

It is generally a good idea to export your information before deleting it, just in case. The HMM hits will show up in most displays, so if you've already run <span class="artifact-n">[anvi-summarize](/software/anvio/help/7.1/programs/anvi-summarize)</span>, you should be good. 

To list available <span class="artifact-n">[hmm-source](/software/anvio/help/7.1/artifacts/hmm-source)</span>s in a database, call 

<div class="codeblock" markdown="1">
anvi&#45;delete&#45;hmms &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                 &#45;&#45;list&#45;hmm&#45;sources
</div>

Then, you can easily delete <span class="artifact-n">[hmm-hits](/software/anvio/help/7.1/artifacts/hmm-hits)</span> from a specific source with the command

<div class="codeblock" markdown="1">
anvi&#45;delete&#45;hmms &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                 &#45;&#45;hmm&#45;source <span class="artifact&#45;n">[hmm&#45;source](/software/anvio/help/7.1/artifacts/hmm&#45;source)</span> 
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-delete-hmms.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-delete-hmms) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

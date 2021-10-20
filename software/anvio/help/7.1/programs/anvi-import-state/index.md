---
layout: page
title: anvi-import-state [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-import-state
image:
  featurerelative: ../../../images/header.png
  display: true
---

Import an anvi&#x27;o state into a profile database.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[pan-db](../../artifacts/pan-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[state-json](../../artifacts/state-json) <img src="../../images/icons/JSON.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[state](../../artifacts/state) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>


## Usage


This program allows you to import a <span class="artifact-n">[state](/software/anvio/help/7.1/artifacts/state)</span> from a <span class="artifact-n">[state-json](/software/anvio/help/7.1/artifacts/state-json)</span>.

You can run this program on a <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span> or <span class="artifact-n">[pan-db](/software/anvio/help/7.1/artifacts/pan-db)</span> like so: 

<div class="codeblock" markdown="1">
anvi&#45;import&#45;state &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7.1/artifacts/profile&#45;db)</span> \
                  &#45;s <span class="artifact&#45;n">[state&#45;json](/software/anvio/help/7.1/artifacts/state&#45;json)</span> \
                  &#45;n MY_STATE
</div>

This will import the state described in your <span class="artifact-n">[state-json](/software/anvio/help/7.1/artifacts/state-json)</span> into your <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span> with the name `MY_STATE`. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-import-state.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-import-state) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

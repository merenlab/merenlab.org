---
layout: page
title: anvi-script-filter-hmm-hits-table [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-script-filter-hmm-hits-table
image:
  featurerelative: ../../../images/header.png
  display: true
---

Filter weak HMM hits from a given contigs database using a domain hits table reported by `anvi-run-hmms`..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/mschecht.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Matthew Schechter</span><div class="page-author-social-box"><a href="mailto:mschechter@uchicago.edu" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/mschecht_bio" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/mschecht" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[hmm-source](../../artifacts/hmm-source) <img src="../../images/icons/HMM.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[hmm-hits](../../artifacts/hmm-hits) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[hmm-hits](../../artifacts/hmm-hits) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>


## Usage


This program is for filtering a <span class="artifact-n">[hmm-source](/software/anvio/help/7.1/artifacts/hmm-source)</span> from a <span class="artifact-n">[hmm-hits](/software/anvio/help/7.1/artifacts/hmm-hits)</span> in a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> using HMM alignment parameters such as query-coverage and target-coverage. Briefly, the program will remove all records from an <span class="artifact-n">[hmm-source](/software/anvio/help/7.1/artifacts/hmm-source)</span> in the <span class="artifact-n">[hmm-hits](/software/anvio/help/7.1/artifacts/hmm-hits)</span>, then import a new <span class="artifact-n">[hmm-hits](/software/anvio/help/7.1/artifacts/hmm-hits)</span> table into the <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> that was filtered to your specifications.

For this, you first need to ask <span class="artifact-n">[anvi-run-hmms](/software/anvio/help/7.1/programs/anvi-run-hmms)</span> to ask HMMER to report a domain hits table by including `--domain-hits-table` flag in your command:

<div class="codeblock" markdown="1">
anvi&#45;run&#45;hmms &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
              &#45;I Bacteria_71 \
              &#45;&#45;hmmer&#45;output&#45;dir path/to/DOMTABLE.txt
              &#45;&#45;domain&#45;hits&#45;table
</div>

At the end of this run, your HMM hits will be stored in your contigs database as usual. But with the availability of the domain hits table from this run, you can filter out hits from your contigs database using thresholds for query or target coverage of each hit.

For instance following the command above, the command below will remove HMM hits from your contigs database for genes that had less than 90% coverage of the target:

<div class="codeblock" markdown="1">
anvi&#45;script&#45;filter&#45;hmm&#45;hits&#45;table &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                                  &#45;&#45;hmm&#45;source Bacteria_71 \
                                  &#45;&#45;domain&#45;hits&#45;table path/to/DOMTABLE.txt \
                                  &#45;&#45;target&#45;coverage 0.9
</div>


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-script-filter-hmm-hits-table.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-filter-hmm-hits-table) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

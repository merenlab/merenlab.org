---
layout: page
title: anvi-script-process-genbank [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-script-process-genbank
image:
  featurerelative: ../../../images/header.png
  display: true
---

This script takes a GenBank file, and outputs a FASTA file, as well as two additional TAB-delimited output files for external gene calls and gene functions that can be used with the programs `anvi-gen-contigs-database` and `anvi-import-functions`.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/AstrobioMike.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Mike Lee</span><div class="page-author-social-box"><a href="https://astrobiomike.github.io" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:michael.lee0517@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/AstrobioMike" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/AstrobioMike" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[genbank-file](../../artifacts/genbank-file) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[contigs-fasta](../../artifacts/contigs-fasta) <img src="../../images/icons/FASTA.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[external-gene-calls](../../artifacts/external-gene-calls) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[functions-txt](../../artifacts/functions-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Usage


This program extracts the data from a <span class="artifact-n">[genbank-file](/software/anvio/help/7.1/artifacts/genbank-file)</span> and converts it into anvi'o friendly artifacts: namely, a <span class="artifact-n">[contigs-fasta](/software/anvio/help/7.1/artifacts/contigs-fasta)</span>, <span class="artifact-n">[external-gene-calls](/software/anvio/help/7.1/artifacts/external-gene-calls)</span> and a <span class="artifact-n">[functions-txt](/software/anvio/help/7.1/artifacts/functions-txt)</span>.

The <span class="artifact-n">[contigs-fasta](/software/anvio/help/7.1/artifacts/contigs-fasta)</span> and <span class="artifact-n">[external-gene-calls](/software/anvio/help/7.1/artifacts/external-gene-calls)</span> can be given to <span class="artifact-n">[anvi-gen-contigs-database](/software/anvio/help/7.1/programs/anvi-gen-contigs-database)</span> to create a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>, and then you can use <span class="artifact-n">[anvi-import-functions](/software/anvio/help/7.1/programs/anvi-import-functions)</span> to bring the function data (in the <span class="artifact-n">[functions-txt](/software/anvio/help/7.1/artifacts/functions-txt)</span>) into the database. Then you'll have all of the data in your <span class="artifact-n">[genbank-file](/software/anvio/help/7.1/artifacts/genbank-file)</span> converted into a single <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>, which you can use for a variety of anvi'o analyses.

The parameters of this program entirely deal with the outputs. Besides telling the program where to put them, you can also give the function annotation source (in the <span class="artifact-n">[functions-txt](/software/anvio/help/7.1/artifacts/functions-txt)</span>) a custom name. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-script-process-genbank.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-process-genbank) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

---
layout: page
title: anvi-run-workflow [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-run-workflow
image:
  featurerelative: ../../../images/header.png
  display: true
---

Execute, manage, parallelize, and troubleshoot entire &#x27;omics workflows and chain together anvi&#x27;o and third party programs.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/ShaiberAlon.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Alon Shaiber</span><div class="page-author-social-box"><a href="mailto:alon.shaiber@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/alon_shaiber" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/ShaiberAlon" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/semiller10.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Samuel Miller</span><div class="page-author-social-box"><a href="https://semiller10.github.io" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:samuelmiller@uchicago.edu" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/smiller_science" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/semiller10" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/mschecht.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Matthew Schechter</span><div class="page-author-social-box"><a href="mailto:mschechter@uchicago.edu" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/mschecht_bio" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/mschecht" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[samples-txt](../../artifacts/samples-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[fasta-txt](../../artifacts/fasta-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[workflow-config](../../artifacts/workflow-config) <img src="../../images/icons/JSON.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[contigs-workflow](../../artifacts/contigs-workflow) <img src="../../images/icons/WORKFLOW.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[metagenomics-workflow](../../artifacts/metagenomics-workflow) <img src="../../images/icons/WORKFLOW.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[pangenomics-workflow](../../artifacts/pangenomics-workflow) <img src="../../images/icons/WORKFLOW.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[phylogenomics-workflow](../../artifacts/phylogenomics-workflow) <img src="../../images/icons/WORKFLOW.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[trnaseq-workflow](../../artifacts/trnaseq-workflow) <img src="../../images/icons/WORKFLOW.png" class="artifact-icon-mini" /></span></p>


## Usage


This program allows you to run [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflows for common anvi'o processes. It is described fully in [this tutorial](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#contigs-workflow). 

Essentially, an anvi'o workflow will run several anvi'o programs for you in quick succession (based on a standard set of intiial steps that will allow you to quickly get to a point where you can ask novel questions). 

As of now, the available workflows are the <span class="artifact-n">[contigs-workflow](/software/anvio/help/7.1/artifacts/contigs-workflow)</span>, the <span class="artifact-n">[metagenomics-workflow](/software/anvio/help/7.1/artifacts/metagenomics-workflow)</span>, the <span class="artifact-n">[pangenomics-workflow](/software/anvio/help/7.1/artifacts/pangenomics-workflow)</span>, the <span class="artifact-n">[phylogenomics-workflow](/software/anvio/help/7.1/artifacts/phylogenomics-workflow)</span>, and the <span class="artifact-n">[trnaseq-workflow](/software/anvio/help/7.1/artifacts/trnaseq-workflow)</span>. 

### Before running the workflow

Each workflow requires a <span class="artifact-n">[workflow-config](/software/anvio/help/7.1/artifacts/workflow-config)</span>: the file that details all of the parameters for the workflow. To get the <span class="artifact-n">[workflow-config](/software/anvio/help/7.1/artifacts/workflow-config)</span> with the default parameters, just run 

<div class="codeblock" markdown="1">
anvi&#45;run&#45;workflow &#45;w WORKFLOW&#45;NAME \
                  &#45;&#45;get&#45;default&#45;config CONFIG.json
</div>

Before running a workflow, it is also a good idea to check the required dependencies by running 

<div class="codeblock" markdown="1">
anvi&#45;run&#45;workflow &#45;w WORKFLOW&#45;NAME \
                  &#45;&#45;list&#45;dependencies
</div>

### The main run 

The main run of the workflow should look like this: 

<div class="codeblock" markdown="1">
anvi&#45;run&#45;workflow &#45;w WORKFLOW&#45;NAME \
                  &#45;c CONFIG.json
                  &#45;&#45;save&#45;workflow&#45;graph
</div>

The flag `--save-workflow-graph` creates a visual representation of the anvio programs that the workflow you're running used. 

You can also use the `-A` flag at the end of the parameter list to change other [Snakemake](https://snakemake.readthedocs.io/en/stable/) parameters. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-run-workflow.md) to update this information.


## Additional Resources


* [Tutorial](http://merenlab.org/2018/07/09/anvio-snakemake-workflows/)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-run-workflow) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

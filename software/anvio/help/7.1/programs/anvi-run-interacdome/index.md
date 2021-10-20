---
layout: page
title: anvi-run-interacdome [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-run-interacdome
image:
  featurerelative: ../../../images/header.png
  display: true
---

Run InteracDome on a contigs database.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/ekiefl.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Evan Kiefl</span><div class="page-author-social-box"><a href="mailto:kiefl.evan@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/evankiefl" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/ekiefl" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[interacdome-data](../../artifacts/interacdome-data) <img src="../../images/icons/DATA.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[binding-frequencies-txt](../../artifacts/binding-frequencies-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[misc-data-amino-acids](../../artifacts/misc-data-amino-acids) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>


## Usage



This program predicts per-residue binding scores for genes in your <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> via the [InteracDome](https://interacdome.princeton.edu/) database.


The full process is detailed in [this blog post](https://merenlab.org/2020/07/22/interacdome/). In fact, ideally, all of that information should really be in this very document, but because the blogpost has preceded this document, it hasn't been translated over yet. So really, you should really be reading that blogpost if you want to get into the nitty gritty details. Otherwise, the quick reference herein should be sufficient.


In summary, this program runs an HMM search of the genes in your <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> to all the Pfam gene families that have been annotated with InteracDome binding frequencies. Then, it parses and filters results, associates binding frequencies of HMM match states to the user's genes of interest, and then stores the resulting per-residue binding frequencies for each gene into the <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> as <span class="artifact-n">[misc-data-amino-acids](/software/anvio/help/7.1/artifacts/misc-data-amino-acids)</span>.


Before running this program, you'll have to run <span class="artifact-n">[anvi-setup-interacdome](/software/anvio/help/7.1/programs/anvi-setup-interacdome)</span> to set up a local copy of [InteracDome's tab-separated files](https://interacdome.princeton.edu/#tab-6136-4).



## Basic Usage

A basic run of this program looks like this:

<div class="codeblock" markdown="1">
anvi&#45;run&#45;interacdome &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> &#45;T 4
</div>

In addition to storing per-residue binding frequencies as <span class="artifact-n">[misc-data-amino-acids](/software/anvio/help/7.1/artifacts/misc-data-amino-acids)</span> in your <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>, this also outputs additional files prefixed with `INTERACDOME` by default (the prefix can be changed with `-O`). These are provided as <span class="artifact-n">[binding-frequencies-txt](/software/anvio/help/7.1/artifacts/binding-frequencies-txt)</span> files named `INTERACDOME-match_state_contributors.txt` and `INTERACDOME-domain_hits.txt`. See <span class="artifact-n">[binding-frequencies-txt](/software/anvio/help/7.1/artifacts/binding-frequencies-txt)</span> for details.


## Parameters

[InteracDome](https://interacdome.princeton.edu/) offers two different binding frequency datasets that can be chosen with `--interacdome-dataset`.  Choose 'representable' to include Pfams that correspond to domain-ligand interactions that had nonredundant instances across three or more distinct PDB structures. InteracDome authors recommend using this collection to learn more about domain binding properties. Choose 'confident' to include Pfams that correspond to domain-ligand interactions that had nonredundant instances across three or more distinct PDB entries and achieved a cross-validated precision of at least 0.5. The default is 'representable', and you can change it like so:


<div class="codeblock" markdown="1">
anvi&#45;run&#45;interacdome &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                     &#45;&#45;interacdome&#45;dataset confident
</div>

This progarm is multi-threaded, so be sure to make use of it:

<div class="codeblock" markdown="1">
anvi&#45;run&#45;interacdome &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                     &#45;&#45;interacdome&#45;dataset confident \
                     &#45;T 8
</div>

Additionally, there are numerous thresholds that you can set: 

1. [`--min-binding-frequency` to ignore very low frequencies](https://merenlab.org/2020/07/22/interacdome/#filtering-low-binding-frequency-scores). The InteracDome scale is from 0 (most likely not involved in binding) to 1 (most likely involved in binding). The default cutoff is 0.200000. 
2. [`--min-hit-fraction` to remove poor quality HMM hits]((https://merenlab.org/2020/07/22/interacdome/#filtering-partial-hits)). The default value is 0.5, so at least half of a profile HMM's length must align to your gene, otherwise the hit will be discarded.
3. [`--information-content-cutoff` to ignore low-qulaity domain hits](https://merenlab.org/2020/07/22/interacdome/#filtering-bad-hits-with-information-content). The default value is 4, which means every amino acid of your gene must match the consensus amino acid of the match state for each mate state with [information content](https://en.wikipedia.org/wiki/Sequence_logo) greater than 4. Decreasing this cutoff yields an increasingly stringent filter.




{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-run-interacdome.md) to update this information.


## Additional Resources


* [Estimating per-residue binding frequencies with InteracDome](http://merenlab.org/2020/07/22/interacdome/)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-run-interacdome) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

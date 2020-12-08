---
layout: page
title: anvi-run-interacdome [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Run InteracDome on a contigs database.

See **[program help menu](../../../../vignette#anvi-run-interacdome)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[binding-frequencies-txt](../../artifacts/binding-frequencies-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[misc-data-nucleotides](../../artifacts/misc-data-nucleotides) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[misc-data-amino-acids](../../artifacts/misc-data-amino-acids) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[interacdome-data](../../artifacts/interacdome-data) <img src="../../images/icons/DATA.png" class="artifact-icon-mini" /></span></p>

## Usage


This program runs [InteracDome](https://interacdome.princeton.edu/) on your <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span>, which **finds the per-residue binding scores for all of your genes**. 

The full process that this program goes through is detailed in [this blog post by Evan Kiefl](https://merenlab.org/2020/07/22/interacdome/). In summary, this program runs the HMM search against all of the genes in your <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span>, parses and filters the results, and then stores the per-residue binding frequencies for each gene into the <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span>.  

Before running this program, you'll have to have run <span class="artifact-n">[anvi-setup-interacdome](/software/anvio/help/main/programs/anvi-setup-interacdome)</span> to set up a local copy of [InteracDome's tab-separated files](https://interacdome.princeton.edu/#tab-6136-4).

### Parameters

A basic run of this program looks like this: 

<div class="codeblock" markdown="1">
anvi&#45;run&#45;interacdome &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/main/artifacts/contigs&#45;db)</span> \
</div>

If you want to annotate potential ligand-binding positions in your sequences instead of domain-binding properties, you can choose to only use Pfams that correspond to domain-ligand interactions that had nonredundant instances across three or more distinct PDB entries and achieved a cross-validated precision of at least 0.5. 

<div class="codeblock" markdown="1">
anvi&#45;run&#45;interacdome &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/main/artifacts/contigs&#45;db)</span> \
                     &#45;&#45;interacdome&#45;dataset confident
</div>

Additionally, there are three thresholds that you can set: 

1. [`--min-binding-frequency` to ignore very low frequencies](https://merenlab.org/2020/07/22/interacdome/#filtering-low-binding-frequency-scores). The InteracDome scale is from 0 (most likely not involved in binding) to 1 (most likely involved in binding). The default cutoff is 0.200000. 
2. [`--min-hit-fraction` to remove results with low detection]((https://merenlab.org/2020/07/22/interacdome/#filtering-partial-hits)). The default value is 0.5, so HMMs that are less than twice as long as the total hit length will not be considered. 
3. [`--information-content-cutoff` to ignore low-qulaity domain hits](https://merenlab.org/2020/07/22/interacdome/#filtering-bad-hits-with-information-content). The default value is 4, so for an alignment to count, the HMM sequence must be very conserved (with more than 95 percent consensus). Setting this cutoff to a very high value will keep all sequences regardless of percent consensus. 

By default, this program does not produce an output, just puts the final results into the <span class="artifact-n">[misc-data-amino-acids](/software/anvio/help/main/artifacts/misc-data-amino-acids)</span> (for SAAVs) or <span class="artifact-n">[misc-data-nucleotides](/software/anvio/help/main/artifacts/misc-data-nucleotides)</span> (for SNVs and SCVs) of your <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span>. If you'd like a tab-delimited output, simply provide the `-O` flag and anvi'o will create a <span class="artifact-n">[binding-frequencies-txt](/software/anvio/help/main/artifacts/binding-frequencies-txt)</span> for your data.  


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-run-interacdome.md) to update this information.


## Additional Resources


* [Estimating per-residue binding frequencies with InteracDome](http://merenlab.org/2020/07/22/interacdome/)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-run-interacdome) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

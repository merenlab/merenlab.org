---
layout: page
title: anvi-run-trna-taxonomy [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

The purpose of this program is to affiliate tRNA gene sequences in an anvi&#x27;o contigs database with taxonomic names. A properly setup local tRNA taxonomy database is required for this program to perform properly. After its successful run, `anvi-estimate-trna-taxonomy` will be useful to estimate taxonomy at genome-, collection-, or metagenome-level)..

See **[program help menu](../../../../vignette#anvi-run-trna-taxonomy)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[trna-taxonomy](../../artifacts/trna-taxonomy) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[trna-taxonomy-db](../../artifacts/trna-taxonomy-db) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>

## Usage


This program associates the tRNA reads found in your <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> with taxonomy information. 

Once these associations are stored in your <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> (represented by a <span class="artifact-n">[trna-taxonomy](/software/anvio/help/7/artifacts/trna-taxonomy)</span> artifact), you'll be able to run <span class="artifact-n">[anvi-estimate-trna-taxonomy](/software/anvio/help/7/programs/anvi-estimate-trna-taxonomy)</span> to use the associations to estimate taxonomy on a larger scale (i.e. for a genome or metagenome). 

To run this program, you'll need to have set up two things: 
1. a <span class="artifact-n">[trna-taxonomy-db](/software/anvio/help/7/artifacts/trna-taxonomy-db)</span>, which you can set up by running <span class="artifact-n">[anvi-setup-trna-taxonomy](/software/anvio/help/7/programs/anvi-setup-trna-taxonomy)</span>.
2. the 'transfer-RNAs' HMM hits in your <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span>, which you can set up by running <span class="artifact-n">[anvi-scan-trnas](/software/anvio/help/7/programs/anvi-scan-trnas)</span>

This program will then go through the tRNA hits in your contigs database and search them against the sequences in the [GTDB](https://gtdb.ecogenomic.org/) databases that you downloaded to assign them taxonomy. 

### Basic run

The following is a basic run of this program: 

<div class="codeblock" markdown="1">
anvi&#45;run&#45;trna&#45;taxonomy &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span>
</div>

If you have set up the two requirements listed above, this should run smoothly. 

### Additional Parameters

When changing these parameters, it might be a good idea to run <span class="artifact-n">[anvi-estimate-trna-taxonomy](/software/anvio/help/7/programs/anvi-estimate-trna-taxonomy)</span> with the `--debug` flag so that you can see what your results look like under the hood. 

1. `--max-num-target-sequences`: the number of hits that this program considers for each tRNA sequence before making a final decision for the taxonomy association. The default is 100, but if you want to ensure that you have accurate data at the expense of some runtime, you can increase it. 
2. `--min-percent-identity`: the minimum percent alignment needed to consider something a hit.  The default is 90, but if you're not getting any hits on a specific sequence, you can decrease it at the risk of getting some nonsense results. 

Finally, this program does not usually have an output file, but if desired you can add the parameter `--all-hits-output-file` to store the list of hits that anvi'o looked at to determine the consensus hit for each sequence. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-run-trna-taxonomy.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-run-trna-taxonomy) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

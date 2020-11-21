---
layout: page
title: anvi-gen-structure-database [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Identifies genes in your contigs database that encode proteins that are homologous to proteins with solved structures. If sufficiently similar homologs are identified, they are used as structural templates to predict the 3D structure of proteins in your contigs database.

See **[program help menu](../../../vignette#anvi-gen-structure-database)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[structure-db](../../artifacts/structure-db)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db)</span> <span class="artifact-r">[pdb-db](../../artifacts/pdb-db)</span></p>

## Usage


This program attempts to solve for the 3D strucutres of proteins encoded by genes in your <span class="artifact-n">[contigs-db](/software/anvio/help/artifacts/contigs-db)</span> using DIAMOND and MODELLER. 

MODELLER first searches your sequence(s) against a database of proteins with a known structure (in Anvi'o, this is either your <span class="artifact-n">[pdb-db](/software/anvio/help/artifacts/pdb-db)</span> or the online copy of [the RCSB database](https://www.rcsb.org/) using [DIAMOND](http://www.diamondsearch.org/index.php). After sequence alignments, the program will select a base template based on the best hits. Then, the program creates a 3D alignment for your sequence and makes final adjustments to it based off of intermolecular interactions. For more information, see [here](http://merenlab.org/2018/09/04/getting-started-with-anvi-3dev/#how-modeller-works). 

The output of this is a <span class="artifact-n">[structure-db](/software/anvio/help/artifacts/structure-db)</span>, which can be used to run <span class="artifact-n">[anvi-3dev](/software/anvio/help/programs/anvi-3dev)</span> to visualize all of this information. You can also export your strucutres into external .pdb files (<span class="artifact-n">[anvi-export-structures](/software/anvio/help/programs/anvi-export-structures)</span>), generate the fixation index matrix (<span class="artifact-n">[anvi-gen-fixation-index-matrix](/software/anvio/help/programs/anvi-gen-fixation-index-matrix)</span>), or the variability profile (<span class="artifact-n">[anvi-gen-variability-profile](/software/anvio/help/programs/anvi-gen-variability-profile)</span>). 

### Basic run 

Here is a simple run:

<div class="codeblock" markdown="1">
anvi&#45;gen&#45;structure&#45;database &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/artifacts/contigs&#45;db)</span> \
                            &#45;&#45;gene&#45;caller&#45;ids 1,2,3 \
                            &#45;o STRUCTURE.db 
</div>

Following this, you will have the strucutures for genes 1, 2, and 3 stored in `STRUCTURE.db`. Alternatively, you can provide a file name with the gene caller-ids (one ID per line) with the flag `--genes-of-interest`. 

To indicate that you have already run <span class="artifact-n">[anvi-setup-pdb-database](/software/anvio/help/programs/anvi-setup-pdb-database)</span> to set up a local copy of representative PDB structures, provide the path to your <span class="artifact-n">[pdb-db](/software/anvio/help/artifacts/pdb-db)</span>:

<div class="codeblock" markdown="1">
anvi&#45;gen&#45;structure&#45;database &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/artifacts/contigs&#45;db)</span> \
                            &#45;&#45;gene&#45;caller&#45;ids 1,2,3 \
                            &#45;&#45;pdb&#45;database <span class="artifact&#45;n">[pdb&#45;db](/software/anvio/help/artifacts/pdb&#45;db)</span> \
                            &#45;o STRUCTURE.db 
</div>

To quickly get a very rough estimate for your structures, you can run with the flag `--very-fast`. 

### Advanced Parameters

Here, we will go through a brief overview of the MODELLER parameters that you are able to change. See [this page](http://merenlab.org/2018/09/04/getting-started-with-anvi-3dev/#description-of-all-modeller-parameters) for more information. 

- The number of models to be simultated. The default is 1. 
- The standard deviation of atomic perturbation of the initial strucutre (i.e. how much you change the position of the atoms before fine tuning with other analysis). The default is 4.
- The MODELLER database used. The default is `pdb_95`, which can be found [here](https://salilab.org/modeller/supplemental.html). This is the same database that is downloaded by <span class="artifact-n">[anvi-setup-pdb-database](/software/anvio/help/programs/anvi-setup-pdb-database)</span>. 
- The socring function used to compare potential models. The default is `DOPE_score`.
- The normalized percent identity cutoff for a template from the database to be further considered. 
- The maximum number of templates that the program will consider. The default is 5. 
- The MODELLER program to use. The default is `mod9.19`. 

For a case study on how some of these parameters matter, see [here](http://merenlab.org/2018/09/04/getting-started-with-anvi-3dev/#a-quick-case-study-on-the-importance-of-key-parameters). 

You also have the option to 

- Skip the use of DSSP, which predicts the locations of beta sheets, alpha helices, etc. 
- Turn on `offline mode`, which will prevent the system from looking up sequences it did not find a match to if it cannot find one in the downloaded database. 

Additionally, to output the raw data, just provide a path to the desired directory with the flag `--dump-dir`. 




{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-gen-structure-database.md) to update this information.


## Additional Resources


* [A conceptual tutorial on the structural biology capabilities of anvio](http://merenlab.org/2018/09/04/structural-biology-with-anvio/)

* [A practical tutorial section in the infant gut tutorial](http://merenlab.org/tutorials/infant-gut/#chapter-vii-linking-genomic-heterogeneity-to-protein-structures)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-gen-structure-database) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

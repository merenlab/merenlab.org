---
layout: page
title: anvi-gen-structure-database [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-gen-structure-database
image:
  featurerelative: ../../../images/header.png
  display: true
---

Creates a database of protein structures. Predict protein structures using template-based homology modelling of genes in your contigs database, or import pre-computed PDB structures you already have..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/ekiefl.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Evan Kiefl</span><div class="page-author-social-box"><a href="mailto:kiefl.evan@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/evankiefl" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/ekiefl" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[pdb-db](../../artifacts/pdb-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[structure-db](../../artifacts/structure-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>


## Usage



This program creates a <span class="artifact-n">[structure-db](/software/anvio/help/7.1/artifacts/structure-db)</span> either by (a) attempting to solve for the 3D structures of proteins encoded by genes in your <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> using DIAMOND and MODELLER, or (b) importing pre-existing structures provided by the user using an <span class="artifact-n">[external-structures](/software/anvio/help/7.1/artifacts/external-structures)</span> file.

### The basics of the pipeline

This section covers option (a), where the user is interested in having structures predicted for them.

DIAMOND first searches your sequence(s) against a database of proteins with a known structure.  This database is downloaded from the [Sali lab](https://salilab.org/modeller/supplemental.html), who created and maintain MODELLER, and contains all of the PDB sequences clustered at 95% identity.

If any good hits are found, they are selected as templates, and their structures are nabbed either from [the RCSB directly](https://www.rcsb.org/), or from a local <span class="artifact-n">[pdb-db](/software/anvio/help/7.1/artifacts/pdb-db)</span> database which you can create yourself with <span class="artifact-n">[anvi-setup-pdb-database](/software/anvio/help/7.1/programs/anvi-setup-pdb-database)</span>. Then, anvi'o passes control over to MODELLER, which creates a 3D alignment for your sequence to the template structures, and makes final adjustments to it based off of empirical distributions of bond angles. For more information, check [this blogpost](http://merenlab.org/2018/09/04/getting-started-with-anvio-structure/#how-modeller-works).

The output of this program is a <span class="artifact-n">[structure-db](/software/anvio/help/7.1/artifacts/structure-db)</span>, which contains all of the modelled structures. Currently, the primary use of the <span class="artifact-n">[structure-db](/software/anvio/help/7.1/artifacts/structure-db)</span> is for interactive exploration with <span class="artifact-n">[anvi-display-structure](/software/anvio/help/7.1/programs/anvi-display-structure)</span>. You can also export your structures into external .pdb files with <span class="artifact-n">[anvi-export-structures](/software/anvio/help/7.1/programs/anvi-export-structures)</span>, or incorporate structural information in the <span class="artifact-n">[variability-profile-txt](/software/anvio/help/7.1/artifacts/variability-profile-txt)</span> with <span class="artifact-n">[anvi-gen-variability-profile](/software/anvio/help/7.1/programs/anvi-gen-variability-profile)</span>.

### Basic standard run

Here is a simple run: 

<div class="codeblock" markdown="1">
anvi&#45;gen&#45;structure&#45;database &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                            &#45;&#45;gene&#45;caller&#45;ids 1,2,3 \
                            &#45;o STRUCTURE.db 
</div>

Following this, you will have the structures for genes 1, 2, and 3 stored in `STRUCTURE.db`, assuming reasonable templates were found. Alternatively, you can provide a file name with the gene caller IDs (one ID per line) with the flag `--genes-of-interest`.  

If you have already run <span class="artifact-n">[anvi-setup-pdb-database](/software/anvio/help/7.1/programs/anvi-setup-pdb-database)</span> and therefore have a local copy of representative PDB structures, make sure you use it by providing the `--offline` flag. If you put it in a non-default location, provide the path to your <span class="artifact-n">[pdb-db](/software/anvio/help/7.1/artifacts/pdb-db)</span>: 

<div class="codeblock" markdown="1">
anvi&#45;gen&#45;structure&#45;database &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                            &#45;&#45;gene&#45;caller&#45;ids 1,2,3 \
                            &#45;&#45;pdb&#45;database <span class="artifact&#45;n">[pdb&#45;db](/software/anvio/help/7.1/artifacts/pdb&#45;db)</span> \
                            &#45;o STRUCTURE.db 
</div>

To quickly get a very rough estimate for your structures, you can run with the flag `--very-fast`. 

### Basic import run

If you already possess structures and would like to create a <span class="artifact-n">[structure-db](/software/anvio/help/7.1/artifacts/structure-db)</span> for downstream anvi'o uses such as <span class="artifact-n">[anvi-display-structure](/software/anvio/help/7.1/programs/anvi-display-structure)</span>, you should create a <span class="artifact-n">[external-structures](/software/anvio/help/7.1/artifacts/external-structures)</span> file. Then, create the database as follows:

<div class="codeblock" markdown="1">
anvi&#45;gen&#45;structure&#45;database &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                            &#45;&#45;external&#45;structures <span class="artifact&#45;n">[external&#45;structures](/software/anvio/help/7.1/artifacts/external&#45;structures)</span> \
                            &#45;o STRUCTURE.db 
</div>

{:.notice}
Please avoid using any MODELLER-specific parameters when using this mode, as they will be silently ignored.


### Advanced Parameters

Here, we will go through a brief overview of the MODELLER parameters that you are able to change. See [this page](http://merenlab.org/2018/09/04/getting-started-with-anvio-structure/#description-of-all-modeller-parameters) for more information. 

- The number of models to be simulated. The default is 1. 
- The standard deviation of atomic perturbation of the initial structure (i.e. how much you change the position of the atoms before fine tuning with other analysis). The default is 4 angstroms.
- The MODELLER database used. The default is `pdb_95`, which can be found [here](https://salilab.org/modeller/supplemental.html). This is the same database that is downloaded by <span class="artifact-n">[anvi-setup-pdb-database](/software/anvio/help/7.1/programs/anvi-setup-pdb-database)</span>.
- The scoring function used to compare potential models. The default is `DOPE_score`.
- The minimum percent identity cutoff for a template to be further considered.
- The minimum alignment fraction that the sequence is covered by the template in order to be further considered.
- The maximum number of templates that the program will consider. The default is 5. 
- The MODELLER program to use. The default is `mod9.19`, but anvi'o is somewhat intelligent and will
  look for the most recent version it can find.

For a case study on how some of these parameters matter, see [here](http://merenlab.org/2018/09/04/getting-started-with-anvio-structure/#a-quick-case-study-on-the-importance-of-key-parameters). 

You also have the option to

- Skip the use of DSSP, which predicts beta sheets, alpha helices, certain bond angles, and relative
  solvent acessibility of residues.
- Output **all** the raw data, just provide a path to the desired directory with the flag `--dump-dir`.




{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-gen-structure-database.md) to update this information.


## Additional Resources


* [A conceptual tutorial on the structural biology capabilities of anvio](http://merenlab.org/2018/09/04/structural-biology-with-anvio/)

* [A practical tutorial section in the infant gut tutorial](http://merenlab.org/tutorials/infant-gut/#chapter-vii-linking-genomic-heterogeneity-to-protein-structures)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-gen-structure-database) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

---
layout: page
title: genome-taxonomy [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/CONCEPT.png" alt="CONCEPT" style="width:100px; border:none" />

A CONCEPT-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-estimate-scg-taxonomy](../../programs/anvi-estimate-scg-taxonomy)</span> <span class="artifact-p">[anvi-estimate-trna-taxonomy](../../programs/anvi-estimate-trna-taxonomy)</span></p>


## Required or used by


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

This artifact is the output tables that are displayed when you run <span class="artifact-n">[anvi-estimate-scg-taxonomy](/software/anvio/help/7/programs/anvi-estimate-scg-taxonomy)</span> or <span class="artifact-n">[anvi-estimate-trna-taxonomy](/software/anvio/help/7/programs/anvi-estimate-trna-taxonomy)</span>. 

By default, they won't be outputed anywhere, just displayed in the terminal for your viewing pleasure. If you want them in a tab-delimited file (as a <span class="artifact-n">[genome-taxonomy-txt](/software/anvio/help/7/artifacts/genome-taxonomy-txt)</span>), just provide the `-o` or the `-O` prefix and anvi'o will do that for you.

The content of these tables will depend on how you ran <span class="artifact-n">[anvi-estimate-trna-taxonomy](/software/anvio/help/7/programs/anvi-estimate-trna-taxonomy)</span> or <span class="artifact-n">[anvi-estimate-scg-taxonomy](/software/anvio/help/7/programs/anvi-estimate-scg-taxonomy)</span>. [This blog post](http://merenlab.org/2019/10/08/anvio-scg-taxonomy/#estimating-taxonomy-in-the-terminal) gives you examples of what this looks like for each of the input scenarios for anvi-estimate-scg-taxonomy. Anvi-estimate-scg-taxonomy's output is very similar, just with the results coming from different gene types. They will also be briefly described below. 

When you run <span class="artifact-n">[anvi-estimate-scg-taxonomy](/software/anvio/help/7/programs/anvi-estimate-scg-taxonomy)</span> or <span class="artifact-n">[anvi-estimate-scg-taxonomy](/software/anvio/help/7/programs/anvi-estimate-scg-taxonomy)</span> on 

- a single genome, this table will contain a single line telling you the taxonomy estimate for your genome. It will also show the number of single-copy core genes or tRNA genes that support this estimate. If you run the `--debug` flag, it will also display the hits for all of the single-copy core genes.  
- a single metagenome, this table will list all of the hits for the chosen single-copy core gene or anticodon (by default, the one with the most hits) and their taxonomy information.   
- a <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> and <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span> with the flag `--compute-scg-coverages`, additional columns will be added that describe the coverage values for your single-copy core gene or tRNA gene hits across your samples.   
- a <span class="artifact-n">[collection](/software/anvio/help/7/artifacts/collection)</span>, this table will show you each of your bins, and the best taxonomy estimate for each one, similarly to how it's displayed for a run on a single genome. 
- a <span class="artifact-n">[metagenomes](/software/anvio/help/7/artifacts/metagenomes)</span> artifact, this table will give a gene entry ID, its taxonomy, and its corresponding coverage in your metagenomes. This format is essentially identical to the output for a single metagenome. If you provide the flag `--matrix-format`, then it will list taxonomy information in each row, and tell you the coverage of each in each of your metagenomes.   

This may sound confusing, but it is easier to understand when looking at the functionality of <span class="artifact-n">[anvi-estimate-scg-taxonomy](/software/anvio/help/7/programs/anvi-estimate-scg-taxonomy)</span> and the comprehensive examples given on [this page](http://merenlab.org/2019/10/08/anvio-scg-taxonomy/#estimating-taxonomy-in-the-terminal).


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/genome-taxonomy.md) to update this information.


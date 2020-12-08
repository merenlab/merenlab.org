---
layout: page
title: anvi-compute-genome-similarity [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Export sequences from sequence sources and compute a similarity metric (e.g. ANI). If a Pan Database is given anvi&#39;o will write computed output to misc data tables of Pan Database.

See **[program help menu](../../../vignette#anvi-compute-genome-similarity)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[genome-similarity](../../artifacts/genome-similarity)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[external-genomes](../../artifacts/external-genomes)</span> <span class="artifact-r">[internal-genomes](../../artifacts/internal-genomes)</span> <span class="artifact-r">[pan-db](../../artifacts/pan-db)</span></p>

## Usage


This program uses the user's similarity metric of choice to calculate the distance between the input genomes.

The currently available distance metrics include:
- [PyANI](https://github.com/widdowquinn/pyani)) to calculate the average nucleotide identity (ANI) (i.e. what portion of orthologous gene pairs align)
- [fastANI](https://github.com/ParBLiSS/FastANI) also to calcualte the ANI but at a faster speed (at the drawback of a slight reduction in accuracy)
- [sourmash](https://sourmash.readthedocs.io/en/latest/) to calculate the mash distance between genomes

### Input/Output

The expected input is any combination of <span class="artifact-n">[external-genomes](/software/anvio/help/main/artifacts/external-genomes)</span>, <span class="artifact-n">[internal-genomes](/software/anvio/help/main/artifacts/internal-genomes)</span>, and text files that contains paths to <span class="artifact-n">[fasta](/software/anvio/help/main/artifacts/fasta)</span> files that describe each of your genomes. This is a tab-delimited file with two columns (`name` and `path` to the fasta files, each of which is assumed to be a single genome),

You also have the option to provide a <span class="artifact-n">[pan-db](/software/anvio/help/main/artifacts/pan-db)</span>, in which case the output data will be written to the database as <span class="artifact-n">[misc-data-layers](/software/anvio/help/main/artifacts/misc-data-layers)</span> and <span class="artifact-n">[misc-data-layer-orders](/software/anvio/help/main/artifacts/misc-data-layer-orders)</span> data. This was done in the [pangenomic tutorial](http://merenlab.org/2016/11/08/pangenomics-v2/#computing-the-average-nucleotide-identity-for-genomes-and-other-genome-similarity-metrics-too). Else, you'll need to provide an output path for the <span class="artifact-n">[genome-similarity](/software/anvio/help/main/artifacts/genome-similarity)</span> files.

Here is an example run with pyANI from an <span class="artifact-n">[external-genomes](/software/anvio/help/main/artifacts/external-genomes)</span> without any parameter changes:

<div class="codeblock" markdown="1">
anvi&#45;compute&#45;genome&#45;similarity &#45;e <span class="artifact&#45;n">[external&#45;genomes](/software/anvio/help/main/artifacts/external&#45;genomes)</span> \
                               &#45;o path/for/<span class="artifact&#45;n">[genome&#45;similarity](/software/anvio/help/main/artifacts/genome&#45;similarity)</span> \
                               &#45;&#45;program pyANI
</div>

### Genome similarity metrics: parameters

#### pyANI

You have the option to change any of the follow parameters:

- The method used for alignment. The options are:
    - `ANIb` (default): uses [BLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome)+ to align 1020 nt fragments of the inputs
    - `ANIm`: uses [MUMmer](http://mummer.sourceforge.net/) to align
    - `ANIblastall`: Uses legacy [BLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) to align 1020 nt fragments
    - `TETRA`: caclulates tetranucleotide frequencies for each input

- The minimum alignment fraction (all percent identity scores lower than this will be automatically set to 0). The default is 0. When using this, you can also use `--significant-alignment-length` to overwrite this discard for sequences with alignmens longer than a specific absolute length
-Similarly, you can discard all results less than some full percent identity (which accounts for what total percent of the two sequences of interest aligned).

#### fastANI

You can change any of the following fastANI parameters:

* The kmer size. The default is 16.

* The fragement length. The default is 30.

* The minimum number of fragments for a result to count. The default is 50.

#### sourmash

You have the option to change the `kmer-size`. This value should depend on the relationship between your samples. The default is 31 ([as recommended by sourmash for genus-level distances](https://sourmash.readthedocs.io/en/latest/using-sourmash-a-guide.html), but we found that 13 most closely parallels that results from an ANI alignment.

You can also set the compression ratio for your fasta files. Decreasing this from the default (1000) will decrease sensitivity.

### Other Parameters

If requested (i.e. a <span class="artifact-n">[pan-db](/software/anvio/help/main/artifacts/pan-db)</span> is not provided), this program outputs similarity matrix files, which can be clustered into a <span class="artifact-n">[dendrogram](/software/anvio/help/main/artifacts/dendrogram)</span>. You can choose to change the distance metric or linkage algorithm for hierarchical clustering.

If your getting a lot of debug/output messages, you can turn them off with `--just-do-it` or helpfully store them into a file with `--log-file`.



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-compute-genome-similarity.md) to update this information.


## Additional Resources


* [In action in the pangenomic workflow tutorial](http://merenlab.org/2016/11/08/pangenomics-v2/#computing-the-average-nucleotide-identity-for-genomes-and-other-genome-similarity-metrics-too)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-compute-genome-similarity) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

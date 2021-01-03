---
layout: page
title: anvi-script-process-genbank-metadata [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

This script takes the &#x27;metadata&#x27; output of the program `ncbi-genome-download` (see [https://github.com/kblin/ncbi-genome-download](https://github.com/kblin/ncbi-genome-download) for details), and processes each GenBank file found in the metadata file to generate a FASTA file, as well as genes and functions files for each entry. Plus, it autmatically generates a FASTA TXT file descriptor for anvi&#x27;o snakemake workflows. So it is a multi-talented program like that.

See **[program help menu](../../../../vignette#anvi-script-process-genbank-metadata)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[contigs-fasta](../../artifacts/contigs-fasta) <img src="../../images/icons/FASTA.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[functions-txt](../../artifacts/functions-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[external-gene-calls](../../artifacts/external-gene-calls) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"></p>

## Usage


Suppose you have downloaded some genomes from NCBI (using [this](https://github.com/kblin/ncbi-genome-download) incredibly useful program) and you have a metadata table describing those genomes. This program will convert that metadata table into some useful files, namely: a FASTA file of contig sequences, an external gene calls file, and an external functions file for each genome you have downloaded; as well as a single tab-delimited fasta-txt file (like the one shown [here](https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#fastatxt)) describing the path to each of these files for all downloaded genomes (that you can pass directly to a snakemake workflow if you need to). Yay.

### The metadata file

The prerequisite for running this program is to have a tab-delimited metadata file containing information about each of the genomes you downloaded from NCBI. Let's say your download command started like this: `ncbi-genome-download --metadata-table ncbi_metadata.txt -t ....` So for the purposes of this usage tutorial, your metadata file is called `ncbi_metadata.txt`.

In case you are wondering, that file should have a header that looks something like this:
```
assembly_accession	bioproject	biosample	wgs_master	excluded_from_refseq	refseq_category	relation_to_type_material	taxid	species_taxid	organism_name	infraspecific_name	isolate	version_status	assembly_level	release_type	genome_rep	seq_rel_dateasm_name	submitter	gbrs_paired_asm	paired_asm_comp	ftp_path	local_filename
```

### Basic usage

If you run this, all the output files will show up in your current working directory.

<div class="codeblock" markdown="1">
anvi&#45;script&#45;process&#45;genbank&#45;metadata &#45;m ncbi_metadata.txt
</div>

### Choosing an output directory

Alternatively, you can specify a directory in which to generate the output:

<div class="codeblock" markdown="1">
anvi&#45;script&#45;process&#45;genbank&#45;metadata &#45;m ncbi_metadata.txt &#45;o DOWNLOADED_GENOMES
</div>

### Picking a name for the fasta-txt file

The default name for the fasta-txt file is `fasta-input.txt`, but you can change that with the `--output-fasta-txt` parameter.

<div class="codeblock" markdown="1">
anvi&#45;script&#45;process&#45;genbank&#45;metadata &#45;m ncbi_metadata.txt &#45;&#45;output&#45;fasta&#45;txt ncbi_fasta.txt
</div>

### Make a fasta-txt without the gene calls and functions columns

The default columns in the fasta-txt file are:
```
name	path	external_gene_calls	gene_functional_annotation
```

But sometimes, you don't want your downstream snakemake workflow to use those external gene calls or functional annotations files. So to skip adding those columns into the fasta-txt file, you can use the `-E` flag:
<div class="codeblock" markdown="1">
anvi&#45;script&#45;process&#45;genbank&#45;metadata &#45;m ncbi_metadata.txt &#45;&#45;output&#45;fasta&#45;txt ncbi_fasta.txt &#45;E
</div>

Then the fasta-txt will only contain a `name` column and a `path` column.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-script-process-genbank-metadata.md) to update this information.


## Additional Resources


* [A tutorial on using this program to access NCBI genomes for &#x27;omics analyses in Anvi&#x27;o](http://merenlab.org/2019/03/14/ncbi-genome-download-magic/)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-process-genbank-metadata) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

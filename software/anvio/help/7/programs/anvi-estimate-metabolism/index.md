---
layout: page
title: anvi-estimate-metabolism [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Reconstructs metabolic pathways and estimates pathway completeness for a given set of contigs.

See **[program help menu](../../../../vignette#anvi-estimate-metabolism)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[kegg-metabolism](../../artifacts/kegg-metabolism) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[kegg-data](../../artifacts/kegg-data) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[kegg-functions](../../artifacts/kegg-functions) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[collection](../../artifacts/collection) <img src="../../images/icons/COLLECTION.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[bin](../../artifacts/bin) <img src="../../images/icons/BIN.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[external-genomes](../../artifacts/external-genomes) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[internal-genomes](../../artifacts/internal-genomes) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[metagenomes](../../artifacts/metagenomes) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>

## Usage


<span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/7/programs/anvi-estimate-metabolism)</span> predicts the metabolic capabilities of organisms based on their genetic content. It relies upon <span class="artifact-n">[kegg-functions](/software/anvio/help/7/artifacts/kegg-functions)</span> and metabolism information from the KEGG resource, which is stored in a <span class="artifact-n">[modules-db](/software/anvio/help/7/artifacts/modules-db)</span>.

The metabolic pathways that this program currently considers are those defined by KOs in the [KEGG MODULES resource](https://www.genome.jp/kegg/module.html). Each KO represents a gene function, and a KEGG module is a set of KOs that collectively carry out the steps in a metabolic pathway. Therefore, for this to work, you need to have annotated your <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> with hits to the KEGG KOfam database by running <span class="artifact-n">[anvi-run-kegg-kofams](/software/anvio/help/7/programs/anvi-run-kegg-kofams)</span> prior to using this program.

Given a properly annotated <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span>, this program determines which KOs are present and from those determines the completeness of each KEGG module. The results are described in a set of output text files, collectively referred to as <span class="artifact-n">[kegg-metabolism](/software/anvio/help/7/artifacts/kegg-metabolism)</span>.

## Running metabolism estimation on a single contigs database

There are several possible inputs to this program. For single genomes (isolate genomes or MAGs, for example) you can provide a <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span>. If your <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> describes a metagenome rather than a single genome, you can provide the flag `--metagenome-mode`. In metagenome mode, estimation is run on each contig individually - that is, only KOfam hits within the same contig are allowed to contribute to the completeness score of a given KEGG module. Alternatively, if you have binned your metagenome sequences into separate populations and would like metabolism estimation to be run separately on each bin, you can provide a <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span> and a <span class="artifact-n">[collection](/software/anvio/help/7/artifacts/collection)</span>.

This program always takes one or more contigs database(s) as input, but what is in those contigs dbs depends on the context (ie, genome, metagenome, bin). In the case of internal genomes (or bins), is possible to have multiple inputs but only one input contigs db. So for clarity's sake, we sometimes refer to the inputs as 'samples' in the descriptions below. If you are getting confused, just try to remember that a 'sample' can be a genome, a metagenome, or a bin.

### Estimation for a single genome

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db
</div>

### Estimation for a metagenome

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db &#45;&#45;metagenome&#45;mode
</div>

{: .notice}
In metagenome mode, this program will estimate metabolism for each contig in the metagenome separately. This will tend to underestimate module completeness because it is likely that some modules will be broken up across multiple contigs belonging to the same population. If you prefer to instead treat all KOs in the metagenome as belonging to one collective genome, you can do so by simply leaving out the `--metagenome-mode` flag (to effectively pretend that you are doing estimation for a single genome, although in your heart you will know that your contigs database really contains a metagenome). Please note that this will result in the opposite tendency to overestimate module completeness (as the KOs will in reality be coming from multiple different populations), and there will be a lot of redundancy. We are working on improving our estimation algorithm for metagenome mode. In the meantime, if you are worried about the misleading results from either of these situations, we suggest binning your metagenomes first and running estimation for the bins as described below.

### Estimation for bins in a metagenome

You can estimate metabolism for each bin in a <span class="artifact-n">[collection](/software/anvio/help/7/artifacts/collection)</span>:

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db &#45;p PROFILE.db &#45;C COLLECTION_NAME
</div>

You can also provide a specific <span class="artifact-n">[bin](/software/anvio/help/7/artifacts/bin)</span> in that <span class="artifact-n">[collection](/software/anvio/help/7/artifacts/collection)</span> to run on:

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db &#45;p PROFILE.db &#45;C COLLECTION_NAME &#45;b BIN_NAME
</div>

Or, you can provide a specific list of bins in a text file:

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db &#45;p PROFILE.db &#45;C COLLECTION_NAME &#45;B bin_ids.txt
</div>

Each line in the `bin_ids.txt` file should be a bin name from the <span class="artifact-n">[collection](/software/anvio/help/7/artifacts/collection)</span>.

## MULTI-MODE: Running metabolism estimation on multiple contigs databases

If you have a set of contigs databases of the same type (i.e., all of them are single genomes or all are binned metagenomes), you can analyze them all at once. What you need to do is put the relevant information for each <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> into a text file and pass that text file to <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/7/programs/anvi-estimate-metabolism)</span>. The program will then run estimation individually on each contigs database in the file. The estimation results for each database will be aggregated and printed to the same output file(s).

### Estimation for multiple single genomes

Multiple single genomes (also known as <span class="artifact-n">[external-genomes](/software/anvio/help/7/artifacts/external-genomes)</span>) can be analyzed with the same command by providing an external genomes file to <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/7/programs/anvi-estimate-metabolism)</span>. To see the required format for the external genomes file, see <span class="artifact-n">[external-genomes](/software/anvio/help/7/artifacts/external-genomes)</span>.

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;e external&#45;genomes.txt
</div>

### Estimation for multiple metagenomes

Multiple metagenomes can be analyzed with the same command by providing a metagenomes input file. Metagenome mode will be used to analyze each contigs database in the file. To see the required format for the metagenomes file, see <span class="artifact-n">[metagenomes](/software/anvio/help/7/artifacts/metagenomes)</span>.

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;M metagenomes.txt
</div>

### Estimation for multiple bins in different metagenomes

If you have multiple bins (also known as <span class="artifact-n">[internal-genomes](/software/anvio/help/7/artifacts/internal-genomes)</span>) across different collections or even different metagenomes, they can be analyzed with the same command by providing an internal genomes file to <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/7/programs/anvi-estimate-metabolism)</span>. To see the required format for the internal genomes file, see <span class="artifact-n">[internal-genomes](/software/anvio/help/7/artifacts/internal-genomes)</span>.

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;i internal&#45;genomes.txt
</div>

## Adjusting module completion threshold

KEGG module completeness is computed as the percentage of steps in the metabolic pathway that are 'present' based on the KOs found in the contigs database. If this completeness is larger than a certain percentage, then the entire module is considered to be 'present' in the genome or metagenome. By default, this module completion threshold is 0.75; that is, 75 percent of the KOs in a module must have a KOfam hit in the contigs database in order for the module to be considered 'complete' as a whole. This threshold can be adjusted.

### Changing the module completion threshold

In this example, we change the threshold to 50 percent.

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db &#45;&#45;module&#45;completion&#45;threshold 0.5
</div>

## Working with a non-default KEGG data directory
If you have previously annotated your contigs databases using a non-default KEGG data directory with `--kegg-data-dir` (see <span class="artifact-n">[anvi-run-kegg-kofams](/software/anvio/help/7/programs/anvi-run-kegg-kofams)</span>), or have moved the KEGG data directory that you wish to use to a non-default location, then you will need to specify where to find the KEGG data so that this program can use the right one. In that case, this is how you do it:

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db &#45;&#45;kegg&#45;data&#45;dir /path/to/directory/KEGG
</div>

## Controlling output

<span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/7/programs/anvi-estimate-metabolism)</span> can produce a variety of output files. All will be prefixed with the same string, which by default is "kegg-metabolism".

### Changing the output file prefix

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db &#45;O my&#45;cool&#45;prefix
</div>

### Including only complete modules in the output
Remember that module completion threshold? Well, you can use that to control which modules make it into your output files. If you provide the `--only-complete` flag, then any module-related output files will only include modules that have a completeness score at or above the module completion threshold. (This doesn't affect KO-related outputs, for obvious reasons.)

Here is an example of using this flag with long format output (which is the default, as described below, but we are asking for it explicitly here just to be clear):
<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db &#45;&#45;kegg&#45;output&#45;modes modules &#45;&#45;only&#45;complete
</div>

And here is an example of using this flag with matrix output. In this case, we are working with multiple input samples, and the behavior of this flag is slightly different: a module will be included in the matrix if it is at or above the module completion threshold in **at least one sample**. If there are any samples in which that module's completeness is below the threshold, its completeness in that sample will be **represented by a 0.0** in the matrix, regardless of its actual completeness score.
<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;i internal&#45;genomes.txt &#45;&#45;matrix&#45;format &#45;&#45;only&#45;complete
</div>


## Output options
This program has two major output options: long format (tab-delimited) output files and matrices.

### Long Format Output
Long format output has several preset "modes" as well as a "custom" mode in which the user can define the contents of the output file. Multiple modes can be used at once, and each requested "mode" will result in a separate output file. The default output mode is "modules" mode.

You can find more details on the output format by looking at <span class="artifact-n">[kegg-metabolism](/software/anvio/help/7/artifacts/kegg-metabolism)</span>.

**Viewing available output modes**

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db &#45;&#45;list&#45;available&#45;modes
</div>

**Using a non-default output mode**

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db &#45;&#45;kegg&#45;output&#45;modes kofam_hits
</div>

**Using multiple output modes**

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db &#45;&#45;kegg&#45;output&#45;modes kofam_hits,modules
</div>

**Viewing available output headers for 'custom' mode**

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db &#45;&#45;list&#45;available&#45;output&#45;headers
</div>

**Using custom output mode**

Here is an example of defining the modules output to contain columns with the module number, the module name, and the completeness score.

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db &#45;&#45;kegg&#45;output&#45;modes custom &#45;&#45;custom&#45;output&#45;headers kegg_module,module_name,module_is_complete
</div>

**Including modules with 0% completeness in long-format output**

By default, modules with a completeness score of 0 are left out of the output files to save on space. But you can explicitly include them by adding the `--include-zeros` flag.

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;c CONTIGS.db &#45;&#45;kegg&#45;output&#45;modes modules &#45;&#45;include&#45;zeros
</div>


### Matrix Output
Matrix format is only available when working with multiple contigs databases. Several output matrices will be generated, each of which describes one statistic such as module completion score, module presence/absence, or KO hit counts. Rows will describe modules or KOs, columns will describe your input samples (ie genomes, metagenomes, bins), and each cell of the matrix will be the corresponding statistic for a module in a sample. You can see examples of this output format by viewing <span class="artifact-n">[kegg-metabolism](/software/anvio/help/7/artifacts/kegg-metabolism)</span>.

**Obtaining matrix-formatted output**

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;i internal&#45;genomes.txt &#45;&#45;matrix&#45;format
</div>

**Including KEGG metadata in the matrix output**

By default, the matrix output is a matrix ready for use in other computational applications, like visualizing as a heatmap or performing clustering. But you may want to instead have a matrix that is annotated with more information, like the names and categories of each module or the functional annotations of each KO. To include this additional information in the matrix output (as columns that occur before the sample columns), use the `--include-metadata` flag.

<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;i internal&#45;genomes.txt &#45;&#45;matrix&#45;format &#45;&#45;include&#45;metadata
</div>

Note that this flag only works for matrix output because, well, the long-format output inherently includes KEGG metadata.

Note also that you can combine this flag with the `--only-complete` flag, like so:
<div class="codeblock" markdown="1">
anvi&#45;estimate&#45;metabolism &#45;i internal&#45;genomes.txt &#45;&#45;matrix&#45;format &#45;&#45;only&#45;complete &#45;&#45;include&#45;metadata
</div>


## Testing this program
You can see if this program is working by running the following suite of tests, which will check several common use-cases:

<div class="codeblock" markdown="1">
anvi&#45;self&#45;test &#45;&#45;suite metabolism
</div>


## Help! I'm getting version errors!
If you have gotten an error that looks something like this:

```
Config Error: The contigs DB that you are working with has been annotated with a different version of the MODULES.db than you are working with now.
```

This means that the MODULES.db used by <span class="artifact-n">[anvi-run-kegg-kofams](/software/anvio/help/7/programs/anvi-run-kegg-kofams)</span> has different contents (different KOs and/or different modules) than the one you are currently using to estimate metabolism, which would lead to mismatches if metabolism estimation were to continue. There are a few ways this can happen, which of course have different solutions:

1. You annotated your <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> with a former version of <span class="artifact-n">[kegg-data](/software/anvio/help/7/artifacts/kegg-data)</span>, and subsequently set up a new <span class="artifact-n">[anvi-setup-kegg-kofams](/software/anvio/help/7/programs/anvi-setup-kegg-kofams)</span> (probably with the `--kegg-archive` or `--download-from-kegg` options, which get you a non-default version of KEGG data). Then you tried to run <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/7/programs/anvi-estimate-metabolism)</span> with the new <span class="artifact-n">[kegg-data](/software/anvio/help/7/artifacts/kegg-data)</span> version. If this is you, and you have saved your former version of <span class="artifact-n">[kegg-data](/software/anvio/help/7/artifacts/kegg-data)</span> somewhere, then you are in luck - you can simply direct <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/7/programs/anvi-estimate-metabolism)</span> to use the old version of KEGG with `--kegg-data-dir`. If you didn't save it, then unfortunately you will most likely have to re-run <span class="artifact-n">[anvi-run-kegg-kofams](/software/anvio/help/7/programs/anvi-run-kegg-kofams)</span> on your <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> to re-annotate it with the new version before continuing with metabolism estimation.
2. You have multiple versions of <span class="artifact-n">[kegg-data](/software/anvio/help/7/artifacts/kegg-data)</span> on your computer in different locations, and you used different ones for <span class="artifact-n">[anvi-run-kegg-kofams](/software/anvio/help/7/programs/anvi-run-kegg-kofams)</span> and<span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/7/programs/anvi-estimate-metabolism)</span>. If this is what you did, then there is an easy fix - simply find the KEGG data directory containing the MODULES.db with the same content hash (you can use <span class="artifact-n">[anvi-db-info](/software/anvio/help/7/programs/anvi-db-info)</span> on the MODULES.db to find the hash value) as the one used by <span class="artifact-n">[anvi-run-kegg-kofams](/software/anvio/help/7/programs/anvi-run-kegg-kofams)</span> and submit that location with `--kegg-data-dir` to this program.
3. Your collaborator gave you some databases that they annotated with a different version of <span class="artifact-n">[kegg-data](/software/anvio/help/7/artifacts/kegg-data)</span> than you have on your computer. In this case, either you or they (or both) has probably set up a non-default version of <span class="artifact-n">[kegg-data](/software/anvio/help/7/artifacts/kegg-data)</span>. If they were using the default snapshot of KEGG data but you are not, then you'll need to get that version onto your computer using the default usage of <span class="artifact-n">[anvi-setup-kegg-kofams](/software/anvio/help/7/programs/anvi-setup-kegg-kofams)</span>. Otherwise, your collaborator will need to somehow share all or part of their KEGG data directory with you before you can work on their databases. See <span class="artifact-n">[anvi-setup-kegg-kofams](/software/anvio/help/7/programs/anvi-setup-kegg-kofams)</span> for details on how to share non-default setups of <span class="artifact-n">[kegg-data](/software/anvio/help/7/artifacts/kegg-data)</span>.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-estimate-metabolism.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-estimate-metabolism) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

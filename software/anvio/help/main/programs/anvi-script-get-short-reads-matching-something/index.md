---
layout: page
title: anvi-script-get-short-reads-matching-something [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

You give this program one or more FASTQ files and a short sequence, and it returns all short reads from the FASTQ file that matches to it.

See **[program help menu](../../../vignette#anvi-script-get-short-reads-matching-something)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Provides

<p style="text-align: left" markdown="1"><span class="artifact-p">[short-reads-fasta](../../artifacts/short-reads-fasta)</span></p>

## Requires or uses

<p style="text-align: left" markdown="1"><span class="artifact-r">[fasta](../../artifacts/fasta)</span></p>

## Usage


This script takes a FASTQ file and a short input sequence and finds all of the short reads in your fastq file that align to your short sequence. 

The purpose of this is to get back short reads that may be extending into hypervariable regions of genomes, resulting a decreased mappability of short reads in the metagenome given a reference. You often see those areas of genomes as significant dips in coverage, and in most cases with a large number of SNVs. When you provide the downstream conserved sequence, this program allows you to take a better look at those regions at the short read level without any mapping.

To instead get short reads mapping to a gene, use <span class="artifact-n">[anvi-get-short-reads-mapping-to-a-gene](/software/anvio/help/main/programs/anvi-get-short-reads-mapping-to-a-gene)</span>.

Here is an example run of this program with the default parameters, where the user is searching for alignments to `AAAAAAAAAAAA` in the sample named `example_sample` stored in the two fastq files `fastaq_one.fastq` and `fastq_two.fastq`: 

<div class="codeblock" markdown="1">
anvi&#45;script&#45;get&#45;short&#45;reads&#45;matching&#45;something &#45;&#45;match&#45;sequence AAAAAAAAAAAA \
                                               &#45;s example_sample \ 
                                               &#45;O example_sample_AAAAAAA_results
                                               fastaq_one.fastq fastq_two.fastq
</div>

This will output all of the matching sequences into three <span class="artifact-n">[fasta](/software/anvio/help/main/artifacts/fasta)</span> files in the directory `example_sample_AAAAAAA_results`. These <span class="artifact-n">[fasta](/software/anvio/help/main/artifacts/fasta)</span> files differ in their format: (1) raw sequences, (2) the same sequences trimmed to the shortest one, (3) gaps `-` added to eliminate length variation. The last two formats provide downstream possibilities with oligotyping to cluster the short reads from an hypervariable region and estimatate their relative proportion. 

Note that this will only report sequences where the length of the short read after the matching sequence is above a certain threshold. The default is 60. For example, if this dataset has the sequence `TTAAAAAAAAAAAAGGGGGGGGG`, this would not be included in the results, but if the sequence was followed by 60 `G` nucleotides, it would be because the length of the sequence after the match is longer than the threshold. You can change this threshold with the parameter `--min-remainder-length`.

You can also choose to stop the program after it finds a certain number of matches or report the raw sequences instead of trimming them to the relevant sections. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-script-get-short-reads-matching-something.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-get-short-reads-matching-something) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

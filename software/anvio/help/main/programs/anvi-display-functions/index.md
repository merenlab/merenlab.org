---
layout: page
title: anvi-display-functions [program]
categories: [anvio]
comments: false
redirect_from: /m/anvi-display-functions
image:
  featurerelative: ../../../images/header.png
  display: true
---

Start an anvi&#x27;o interactive display to see functions across genomes.

See **[program help menu](../../../../vignette#anvi-display-functions)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[interactive](../../artifacts/interactive) <img src="../../images/icons/DISPLAY.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[functions](../../artifacts/functions) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[genomes-storage-db](../../artifacts/genomes-storage-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[internal-genomes](../../artifacts/internal-genomes) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[external-genomes](../../artifacts/external-genomes) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>

## Usage


For a given annotation source for <span class="artifact-n">[functions](/software/anvio/help/main/artifacts/functions)</span>, this program will display distribution patterns of unique function names (or accession numbers) across genomes stored in anvi'o databases.

It is a powerful way to analyze differentially occurring functions for any source of annotation that is shared across all genomes.

Currently, <span class="artifact-n">[anvi-display-functions](/software/anvio/help/main/programs/anvi-display-functions)</span> can work with any combination of genomes from <span class="artifact-n">[external-genomes](/software/anvio/help/main/artifacts/external-genomes)</span>, <span class="artifact-n">[internal-genomes](/software/anvio/help/main/artifacts/internal-genomes)</span>, and <span class="artifact-n">[genomes-storage-db](/software/anvio/help/main/artifacts/genomes-storage-db)</span>.

### Basic Run 

The simplest way to run this program is as follows:

<div class="codeblock" markdown="1">
anvi&#45;display&#45;functions &#45;e <span class="artifact&#45;n">[external&#45;genomes](/software/anvio/help/main/artifacts/external&#45;genomes)</span> \
                       &#45;&#45;annotation&#45;source KOfam \
                       &#45;&#45;profile&#45;db KOFAM&#45;PROFILE.db
</div>


You can replace the annotation source based on what is available across your genomes. You can use the program <span class="artifact-n">[anvi-db-info](/software/anvio/help/main/programs/anvi-db-info)</span> to see all available function annotation sources in a given <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span> or <span class="artifact-n">[genomes-storage-db](/software/anvio/help/main/artifacts/genomes-storage-db)</span>. Please see <span class="artifact-n">[functions](/software/anvio/help/main/artifacts/functions)</span> for more information on functions and how to obtain them.


{:.notice}
Please note that a <span class="artifact-n">[profile-db](/software/anvio/help/main/artifacts/profile-db)</span> will be automatically generated for you. Once it is generated, the same profile database can be visualized over and over again using <span class="artifact-n">[anvi-interactive](/software/anvio/help/main/programs/anvi-interactive)</span> in manual mode, without having to retain any other files.

### Aggregating functions using accession IDs

Once it is run, this program essentially aggregates all function names that occur in one or more genomes in the set of genomes found in input sources. The user can ask the program to use accession IDs to aggregate functions rather than function names:

<div class="codeblock" markdown="1">
anvi&#45;display&#45;functions &#45;e <span class="artifact&#45;n">[external&#45;genomes](/software/anvio/help/main/artifacts/external&#45;genomes)</span> \
                       &#45;&#45;annotation&#45;source KOfam \
                       &#45;&#45;profile&#45;db KOFAM&#45;PROFILE.db \
                       &#45;&#45;aggregate&#45;based&#45;on&#45;accession
</div>

While the default setting which is to use function names will be appropriate for most applications, using accession IDs instead of function names may be important for specific applications. There may be an actual difference between using functions or accession to aggregate data since multiple accession IDs in various databases may correspond to the same function. This may lead to misleading enrichment analyses downstream as identical function annotations may be over-split into multiple groups. Thus, the default aggregation method uses function names.

### Aggregating functions using accession IDs

In some cases a gene may be annotated with multiple functions. This is a decision often made at the function annotation tool level. For instance <span class="artifact-n">[anvi-run-ncbi-cogs](/software/anvio/help/main/programs/anvi-run-ncbi-cogs)</span> may yield two COG annotations for a single gene because the significance score for both hits may exceed the default cutoff. While this can be useful in <span class="artifact-n">[anvi-summarize](/software/anvio/help/main/programs/anvi-summarize)</span> output where things should be most comprehensive, having some genes annotated with multiple functions and others with one function may over-split them (since in this scenario a gene with COGXXX and COGXXX;COGYYY would end up in different bins). Thus, <span class="artifact-n">[anvi-display-functions](/software/anvio/help/main/programs/anvi-display-functions)</span> will will use the best hit for any gene that has multiple hits. But this behavior can be turned off the following way:

<div class="codeblock" markdown="1">
anvi&#45;display&#45;functions &#45;e <span class="artifact&#45;n">[external&#45;genomes](/software/anvio/help/main/artifacts/external&#45;genomes)</span> \
                       &#45;&#45;annotation&#45;source KOfam \
                       &#45;&#45;profile&#45;db KOFAM&#45;PROFILE.db \
                       &#45;&#45;aggregate&#45;using&#45;all&#45;hits
</div>

---

The user also can keep only functions that occur in more than a minimum number of genomes:

<div class="codeblock" markdown="1">
anvi&#45;display&#45;functions &#45;e <span class="artifact&#45;n">[external&#45;genomes](/software/anvio/help/main/artifacts/external&#45;genomes)</span> \
                       &#45;&#45;annotation&#45;source KOfam \
                       &#45;&#45;profile&#45;db KOFAM&#45;PROFILE.db \
                       &#45;&#45;min&#45;occurrence 5
</div>

### Combining genomes from multiple sources

Alternatively, you can run the program by combining genomes from multiple sources:

<div class="codeblock" markdown="1">
anvi&#45;display&#45;functions &#45;e <span class="artifact&#45;n">[external&#45;genomes](/software/anvio/help/main/artifacts/external&#45;genomes)</span> \
                       &#45;i <span class="artifact&#45;n">[internal&#45;genomes](/software/anvio/help/main/artifacts/internal&#45;genomes)</span> \
                       &#45;g <span class="artifact&#45;n">[genomes&#45;storage&#45;db](/software/anvio/help/main/artifacts/genomes&#45;storage&#45;db)</span> \
                       &#45;&#45;annotation&#45;source KOfam \
                       &#45;&#45;profile&#45;db KOFAM&#45;PROFILE.db

</div>

### A real-world example

Assume we have a list of <span class="artifact-n">[external-genomes](/software/anvio/help/main/artifacts/external-genomes)</span> that include three different species of *Bifidobacterium*. Running the following command,

<div class="codeblock" markdown="1">
anvi&#45;display&#45;functions &#45;&#45;external&#45;genomes Bifidobacterium.txt \
                       &#45;&#45;annotation&#45;source COG20_FUNCTION \
                       &#45;&#45;profile&#45;db COG20&#45;PROFILE.db \
                       &#45;&#45;min&#45;occurrence 3
</div>

Would produce the following display by default, where each layer is one of the genomes described in the <span class="artifact-n">[external-genomes](/software/anvio/help/main/artifacts/external-genomes)</span> file, and each item is a unique function name that occur in `COG20_FUNCTION` (which was obtained by running <span class="artifact-n">[anvi-run-ncbi-cogs](/software/anvio/help/main/programs/anvi-run-ncbi-cogs)</span> on each <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span> in the external genomes file) that were found in more than three genomes:

[![Example output](../../images/anvi-display-functions-01.png){:.center-img .width-50}](../../images/anvi-display-functions-01.png)

The outermost layer shows the function names:

[![Example output](../../images/anvi-display-functions-02.png){:.center-img .width-50}](../../images/anvi-display-functions-02.png)

After a quick prettification through the <span class="artifact-n">[interactive](/software/anvio/help/main/artifacts/interactive)</span> interface, leads to a cleaner display of three distinct species in this group, and functions that are uniquely enriched in either of them: 

[![Example output](../../images/anvi-display-functions-03.png){:.center-img .width-80}](../../images/anvi-display-functions-03.png)

Now the resulting <span class="artifact-n">[profile-db](/software/anvio/help/main/artifacts/profile-db)</span> can be used by <span class="artifact-n">[anvi-interactive](/software/anvio/help/main/programs/anvi-interactive)</span> to re-visualize these data, or can be shared with the community without sharing the underlying contigs databases.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-display-functions.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-display-functions) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

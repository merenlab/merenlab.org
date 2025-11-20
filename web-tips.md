---
layout: page
title: Tips for web developers
modified: 2020-12-12
comments: false
---

The purpose of this page is to show some of the fancy things our web page can do. If you add new features that are specific to our web code, or discover uses of features that are not described here, please feel free to extend this list.

# Show thumbnails for blog posts

For a blog post to have its own thumbnail, you need to add the following to the frontmatter of the Markdown file:

```
thumbnail: /images/thumbnails/image.jpg
```

Since thumbnails will be showed for every post in the blog listing page, it is critical to make sure that they have a small file size. The best practice here is the following. Let's assume the markdown file for your blogpost is named as `2025-11-17-mmoff-2025.md`, and it contains an image you wish to use as a thumbnail located at,

```
images/miscellaneous/mmoff-2025/mmoff-2025-workshop-discussion-round.jpg
```

The best way to generate a thumbnail for this post would be to run the following command in your terminal (which will require you to install [ImageMagick](https://imagemagick.org/) if you don't have it already (`apt` or `brew` will work perfectly for this)):

```
magick images/miscellaneous/mmoff-2025/mmoff-2025-workshop-discussion-round.jpg \
       -thumbnail 600x600^ \
       -gravity center \
       -extent 600x600 \
       -quality 85 \
       images/thumbnails/2025-11-17-mmoff-2025.jpg
```

For instance, in this particular case, this command generates a 0.1 MB thumbnail for a 2.8 MB original image. Then you can add the thumnbnail path to the frontmatter of your blog post like this (don't forget the extra `/` at the beginning of the path):

```
thumbnail: /images/thumbnails/2025-11-17-mmoff-2025.jpg
```

After this, you should see your post with a thumbnail in the blog listing page.

# Show/hide content

If you want to show/hide content, you can use this notation in your markdown files.

**Raw**:

``` html
{% raw %}
<details markdown="1"><summary>Show/Hide SOME CONTENT</summary>

SOME CONTENT GOES HERE

</details>
{% endraw %}
```

**Rendered**:

<details markdown="1"><summary>Show/Hide SOME CONTENT</summary>

SOME CONTENT GOES HERE

</details>

# Show/hide other stuff

One can add to the header sections of Markdwon files (pages, posts, etc) served by the web page

```
noleftpanel: true
```

to remove the left panel with all the links;

```
nonavigation: true
```

to remove the navigation links shown at the top, and,

```
image:
  feature: https://merenlab.org/images/the/image/you/want/to/display.png
  display: true
  nologo: true
```

To display a header image, but not the logo of the lab.

[Here](https://merenlab.org/2022/01/03/giant-viruses/) is an example blog post where each of these parameters were used.


## Code file name

When we show code, we usually use the standard markdown syntax and it looks like this:

```sh
echo "Script Name: $0"
echo "First Argument: $1"
echo "Number of Arguments: $#"
echo "Exit Status of Last Command: $?"
```

You can use also indicate the code shown is the contents of a file.

**Raw**:

~~~
{% raw %}
```sh
echo "Script Name: $0"
echo "First Argument: $1"
echo "Number of Arguments: $#"
echo "Exit Status of Last Command: $?"
```
{% include CODEBLOCKFILENAME filename="variables-built-in.sh" %}
{% endraw %}
~~~

**Rendered**:

```sh
echo "Script Name: $0"
echo "First Argument: $1"
echo "Number of Arguments: $#"
echo "Exit Status of Last Command: $?"
```
{% include CODEBLOCKFILENAME filename="variables-built-in.sh" %}



## Information boxes

If you want to show summary sections with a different background color, you can use this notation in your markdown files.

**Raw**:

``` html
{% raw %}
<div class="extra-info" markdown="1">

<span class="extra-info-header">Smart title for extra info</span>

EXTRA INFO GOES HERE

</div>
{% endraw %}
```

**Rendered**:

<div class="extra-info" markdown="1">

<span class="extra-info-header">Smart title for extra info</span>

EXTRA INFO GOES HERE

</div>

## Warning and notice statements

You should feel free to use warning and notice statements.

**Raw**:

``` html
{% raw %}
{:.warning}
A warning messages goes here.

{:.notice}
A notice statement goes here.
{% endraw %}
```

**Rendered**:

{:.warning}
A warning messages goes here.

{:.notice}
A notice statement goes here.

## Turning TAB-delmited files to markdown notation

Tutorial writers often need to display contents of TAB-delimited files. One easy way to do that is to use the program {% include PROGRAM name="anvi-script-as-markdown" %}. See examples on the program page.


## Marking anvi'o artifacts and programs

You can use this notation to link anvi'o programs and artifacts to their documentation with fancy looking links. If you are a `vim` user, you can automatize this conversion. Please see the section in this file called "A useful shortcut for vim users".

**Raw**:

```
{% raw %}
In all fairness, an anvi'o {% include ARTIFACT name="contigs-db" text="contigs database" version="8" %} is just a FASTA file adjusted to 2020 and is generated by {% include PROGRAM name="anvi-gen-contigs-database" text="a fancy anvi'o program" version="8" %}.
{% endraw %}
```

**Rendered**:

In all fairness, an anvi'o {% include ARTIFACT name="contigs-db" text="contigs database" version="8" %} is just a FASTA file adjusted to 2020 and is generated by {% include PROGRAM name="anvi-gen-contigs-database" text="a fancy anvi'o program" version="8" %}.

**Raw**:

```
{% raw %}
An anvi'o {% include ARTIFACT name="pan-db" %} comes from {% include PROGRAM name="anvi-pan-genome" %}.
{% endraw %}
```

**Rendered**:

An anvi'o {% include ARTIFACT name="pan-db" %} comes from {% include PROGRAM name="anvi-pan-genome" %}.

## Fancy quotes

**Raw**:

``` html
<blockquote>
An interesting quote.
<div class="blockquote-author">Someone Interesting</div>
</blockquote>
```

**Rendered**:

<blockquote>
An interesting quote.
<div class="blockquote-author">Someone Interesting</div>
</blockquote>

## Images

Never load your images directly, and always use this notation, which will enable **lazy loading** and will standardize all images:

``` html
{% raw %}
{% include IMAGE path="/path/to/image.png" width=50 %}
{% endraw %}
```

Here are some examples:

``` html
{% raw %}
{% include IMAGE path="https://merenlab.org/logo.png" %}
{% endraw %}
```

``` html
{% raw %}
{% include IMAGE path="some_folder/store_c.gif" width=50 %}
{% endraw %}
```

You can add **title text** for your figures, which will appear as a tiny text when people hover their mouse pointers:

``` html
{% raw %}
{% include IMAGE path="/images/ref_bg.gif" title="Changing the background color of a layer" width=90 %}
{% endraw %}
```

You can also **add captions** to your figures, which will appear as italicized text closer to your figure:

``` html
{% raw %}
{% include IMAGE path="/images/some_image.gif" caption="A descriptive short paragraph to explain what this figure is showing" %}
{% endraw %}
```

If you don't want to have a **border around your figure**, you can also turn it off using the `noborder` parameter:

``` html
{% raw %}
{% include IMAGE path="/images/some_image.gif" noborder="1" %}
{% endraw %}
```



## A useful shortcut for vim users

If you are using vim as a text editor, you can automatically turn anvi'o program and artifact names into nice markings shown below instantenously with no effort as originally suggested by Evan Kiefl. With this solution, when you type the following in any document in `insert` mode (notice the extra `i` at the end of program and artifact names),

```
An anvi'o pan-dbi comes from anvi-pan-genomei.
```

vim will automatically turn it into the following as you write:

```
{% raw %}
An anvi'o {% include ARTIFACT name="pan-db" %} comes from {% include PROGRAM name="anvi-pan-genome" %}.
{% endraw %}
```

Which will then be rendered as "An anvi'o {% include ARTIFACT name="pan-db" %} comes from {% include PROGRAM name="anvi-pan-genome" %}", as you know.

For this to improve your tutorial writing, all you need to do is to add the following lines somewhere in your `~/.vimrc` file:

``` vim
{% raw %}
" Anvi'o artifacts
inoreabbrev pan-dbi {% include ARTIFACT name="pan-db" %}
inoreabbrev contigs-dbi {% include ARTIFACT name="contigs-db" %}
inoreabbrev trnaseq-dbi {% include ARTIFACT name="trnaseq-db" %}
inoreabbrev trnaseq-contigs-dbi {% include ARTIFACT name="trnaseq-contigs-db" %}
inoreabbrev trnaseq-profile-dbi {% include ARTIFACT name="trnaseq-profile-db" %}
inoreabbrev modules-dbi {% include ARTIFACT name="modules-db" %}
inoreabbrev fastai {% include ARTIFACT name="fasta" %}
inoreabbrev contigs-fastai {% include ARTIFACT name="contigs-fasta" %}
inoreabbrev dna-sequencei {% include ARTIFACT name="dna-sequence" %}
inoreabbrev trnaseq-fastai {% include ARTIFACT name="trnaseq-fasta" %}
inoreabbrev configuration-inii {% include ARTIFACT name="configuration-ini" %}
inoreabbrev external-gene-callsi {% include ARTIFACT name="external-gene-calls" %}
inoreabbrev external-structuresi {% include ARTIFACT name="external-structures" %}
inoreabbrev concatenated-gene-alignment-fastai {% include ARTIFACT name="concatenated-gene-alignment-fasta" %}
inoreabbrev short-reads-fastai {% include ARTIFACT name="short-reads-fasta" %}
inoreabbrev genes-fastai {% include ARTIFACT name="genes-fasta" %}
inoreabbrev bam-filei {% include ARTIFACT name="bam-file" %}
inoreabbrev bam-stats-txti {% include ARTIFACT name="bam-stats-txt" %}
inoreabbrev protein-structure-txti {% include ARTIFACT name="protein-structure-txt" %}
inoreabbrev samples-txti {% include ARTIFACT name="samples-txt" %}
inoreabbrev primers-txti {% include ARTIFACT name="primers-txt" %}
inoreabbrev fasta-txti {% include ARTIFACT name="fasta-txt" %}
inoreabbrev raw-bam-filei {% include ARTIFACT name="raw-bam-file" %}
inoreabbrev locus-fastai {% include ARTIFACT name="locus-fasta" %}
inoreabbrev structure-dbi {% include ARTIFACT name="structure-db" %}
inoreabbrev pdb-dbi {% include ARTIFACT name="pdb-db" %}
inoreabbrev kegg-datai {% include ARTIFACT name="kegg-data" %}
inoreabbrev single-profile-dbi {% include ARTIFACT name="single-profile-db" %}
inoreabbrev profile-dbi {% include ARTIFACT name="profile-db" %}
inoreabbrev genes-dbi {% include ARTIFACT name="genes-db" %}
inoreabbrev genomes-storage-dbi {% include ARTIFACT name="genomes-storage-db" %}
inoreabbrev contigs-statsi {% include ARTIFACT name="contigs-stats" %}
inoreabbrev svgi {% include ARTIFACT name="svg" %}
inoreabbrev bini {% include ARTIFACT name="bin" %}
inoreabbrev collectioni {% include ARTIFACT name="collection" %}
inoreabbrev collection-txti {% include ARTIFACT name="collection-txt" %}
inoreabbrev hmm-sourcei {% include ARTIFACT name="hmm-source" %}
inoreabbrev hmm-hitsi {% include ARTIFACT name="hmm-hits" %}
inoreabbrev completioni {% include ARTIFACT name="completion" %}
inoreabbrev cogs-datai {% include ARTIFACT name="cogs-data" %}
inoreabbrev pfams-datai {% include ARTIFACT name="pfams-data" %}
inoreabbrev misc-data-items-txti {% include ARTIFACT name="misc-data-items-txt" %}
inoreabbrev misc-data-itemsi {% include ARTIFACT name="misc-data-items" %}
inoreabbrev misc-data-layers-txti {% include ARTIFACT name="misc-data-layers-txt" %}
inoreabbrev misc-data-layersi {% include ARTIFACT name="misc-data-layers" %}
inoreabbrev misc-data-nucleotides-txti {% include ARTIFACT name="misc-data-nucleotides-txt" %}
inoreabbrev misc-data-nucleotidesi {% include ARTIFACT name="misc-data-nucleotides" %}
inoreabbrev misc-data-amino-acids-txti {% include ARTIFACT name="misc-data-amino-acids-txt" %}
inoreabbrev misc-data-amino-acidsi {% include ARTIFACT name="misc-data-amino-acids" %}
inoreabbrev genome-similarityi {% include ARTIFACT name="genome-similarity" %}
inoreabbrev misc-data-layer-orders-txti {% include ARTIFACT name="misc-data-layer-orders-txt" %}
inoreabbrev misc-data-layer-ordersi {% include ARTIFACT name="misc-data-layer-orders" %}
inoreabbrev misc-data-items-order-txti {% include ARTIFACT name="misc-data-items-order-txt" %}
inoreabbrev misc-data-items-orderi {% include ARTIFACT name="misc-data-items-order" %}
inoreabbrev dendrogrami {% include ARTIFACT name="dendrogram" %}
inoreabbrev metapangenomei {% include ARTIFACT name="metapangenome" %}
inoreabbrev oligotypesi {% include ARTIFACT name="oligotypes" %}
inoreabbrev linkmers-txti {% include ARTIFACT name="linkmers-txt" %}
inoreabbrev palindromes-txti {% include ARTIFACT name="palindromes-txt" %}
inoreabbrev inversionsi {% include ARTIFACT name="inversions" %}
inoreabbrev phylogenyi {% include ARTIFACT name="phylogeny" %}
inoreabbrev gene-calls-txti {% include ARTIFACT name="gene-calls-txt" %}
inoreabbrev binding-frequencies-txti {% include ARTIFACT name="binding-frequencies-txt" %}
inoreabbrev interacdome-datai {% include ARTIFACT name="interacdome-data" %}
inoreabbrev functionsi {% include ARTIFACT name="functions" %}
inoreabbrev functions-txti {% include ARTIFACT name="functions-txt" %}
inoreabbrev functional-enrichment-txti {% include ARTIFACT name="functional-enrichment-txt" %}
inoreabbrev kegg-functionsi {% include ARTIFACT name="kegg-functions" %}
inoreabbrev interactivei {% include ARTIFACT name="interactive" %}
inoreabbrev view-datai {% include ARTIFACT name="view-data" %}
inoreabbrev layer-taxonomyi {% include ARTIFACT name="layer-taxonomy" %}
inoreabbrev layer-taxonomy-txti {% include ARTIFACT name="layer-taxonomy-txt" %}
inoreabbrev gene-taxonomyi {% include ARTIFACT name="gene-taxonomy" %}
inoreabbrev gene-taxonomy-txti {% include ARTIFACT name="gene-taxonomy-txt" %}
inoreabbrev genome-taxonomyi {% include ARTIFACT name="genome-taxonomy" %}
inoreabbrev genome-taxonomy-txti {% include ARTIFACT name="genome-taxonomy-txt" %}
inoreabbrev scgs-taxonomy-dbi {% include ARTIFACT name="scgs-taxonomy-db" %}
inoreabbrev scgs-taxonomyi {% include ARTIFACT name="scgs-taxonomy" %}
inoreabbrev trna-taxonomy-dbi {% include ARTIFACT name="trna-taxonomy-db" %}
inoreabbrev trna-taxonomyi {% include ARTIFACT name="trna-taxonomy" %}
inoreabbrev external-genomesi {% include ARTIFACT name="external-genomes" %}
inoreabbrev internal-genomesi {% include ARTIFACT name="internal-genomes" %}
inoreabbrev metagenomesi {% include ARTIFACT name="metagenomes" %}
inoreabbrev coverages-txti {% include ARTIFACT name="coverages-txt" %}
inoreabbrev detection-txti {% include ARTIFACT name="detection-txt" %}
inoreabbrev variability-profilei {% include ARTIFACT name="variability-profile" %}
inoreabbrev variability-profile-txti {% include ARTIFACT name="variability-profile-txt" %}
inoreabbrev codon-frequencies-txti {% include ARTIFACT name="codon-frequencies-txt" %}
inoreabbrev aa-frequencies-txti {% include ARTIFACT name="aa-frequencies-txt" %}
inoreabbrev fixation-index-matrixi {% include ARTIFACT name="fixation-index-matrix" %}
inoreabbrev trnaseq-seed-txti {% include ARTIFACT name="trnaseq-seed-txt" %}
inoreabbrev seeds-specific-txti {% include ARTIFACT name="seeds-specific-txt" %}
inoreabbrev seeds-non-specific-txti {% include ARTIFACT name="seeds-non-specific-txt" %}
inoreabbrev modifications-txti {% include ARTIFACT name="modifications-txt" %}
inoreabbrev trnaseq-ploti {% include ARTIFACT name="trnaseq-plot" %}
inoreabbrev summaryi {% include ARTIFACT name="summary" %}
inoreabbrev quick-summaryi {% include ARTIFACT name="quick-summary" %}
inoreabbrev split-binsi {% include ARTIFACT name="split-bins" %}
inoreabbrev statei {% include ARTIFACT name="state" %}
inoreabbrev ngramsi {% include ARTIFACT name="ngrams" %}
inoreabbrev state-jsoni {% include ARTIFACT name="state-json" %}
inoreabbrev kegg-metabolismi {% include ARTIFACT name="kegg-metabolism" %}
inoreabbrev augustus-gene-callsi {% include ARTIFACT name="augustus-gene-calls" %}
inoreabbrev genes-statsi {% include ARTIFACT name="genes-stats" %}
inoreabbrev vcfi {% include ARTIFACT name="vcf" %}
inoreabbrev blast-tablei {% include ARTIFACT name="blast-table" %}
inoreabbrev splits-txti {% include ARTIFACT name="splits-txt" %}
inoreabbrev genbank-filei {% include ARTIFACT name="genbank-file" %}
inoreabbrev groups-txti {% include ARTIFACT name="groups-txt" %}
inoreabbrev splits-taxonomy-txti {% include ARTIFACT name="splits-taxonomy-txt" %}
inoreabbrev hmm-hits-matrix-txti {% include ARTIFACT name="hmm-hits-matrix-txt" %}
inoreabbrev pn-ps-datai {% include ARTIFACT name="pn-ps-data" %}
inoreabbrev clustering-configurationi {% include ARTIFACT name="clustering-configuration" %}
inoreabbrev workflow-configi {% include ARTIFACT name="workflow-config" %}
inoreabbrev contigs-workflowi {% include ARTIFACT name="contigs-workflow" %}
inoreabbrev metagenomics-workflowi {% include ARTIFACT name="metagenomics-workflow" %}
inoreabbrev pangenomics-workflowi {% include ARTIFACT name="pangenomics-workflow" %}
inoreabbrev phylogenomics-workflowi {% include ARTIFACT name="phylogenomics-workflow" %}
inoreabbrev trnaseq-workflowi {% include ARTIFACT name="trnaseq-workflow" %}
inoreabbrev contig-inspectioni {% include ARTIFACT name="contig-inspection" %}
inoreabbrev gene-cluster-inspectioni {% include ARTIFACT name="gene-cluster-inspection" %}

" Anvi'o programs
inoreabbrev anvi-analyze-syntenyi {% include PROGRAM name="anvi-analyze-synteny" %}
inoreabbrev anvi-cluster-contigsi {% include PROGRAM name="anvi-cluster-contigs" %}
inoreabbrev anvi-compute-anii {% include PROGRAM name="anvi-compute-ani" %}
inoreabbrev anvi-compute-completenessi {% include PROGRAM name="anvi-compute-completeness" %}
inoreabbrev anvi-compute-functional-enrichmenti {% include PROGRAM name="anvi-compute-functional-enrichment" %}
inoreabbrev anvi-compute-functional-enrichment-across-genomesi {% include PROGRAM name="anvi-compute-functional-enrichment-across-genomes" %}
inoreabbrev anvi-compute-functional-enrichment-in-pani {% include PROGRAM name="anvi-compute-functional-enrichment-in-pan" %}
inoreabbrev anvi-compute-gene-cluster-homogeneityi {% include PROGRAM name="anvi-compute-gene-cluster-homogeneity" %}
inoreabbrev anvi-compute-genome-similarityi {% include PROGRAM name="anvi-compute-genome-similarity" %}
inoreabbrev anvi-compute-metabolic-enrichmenti {% include PROGRAM name="anvi-compute-metabolic-enrichment" %}
inoreabbrev anvi-db-infoi {% include PROGRAM name="anvi-db-info" %}
inoreabbrev anvi-delete-collectioni {% include PROGRAM name="anvi-delete-collection" %}
inoreabbrev anvi-delete-functionsi {% include PROGRAM name="anvi-delete-functions" %}
inoreabbrev anvi-delete-hmmsi {% include PROGRAM name="anvi-delete-hmms" %}
inoreabbrev anvi-delete-misc-datai {% include PROGRAM name="anvi-delete-misc-data" %}
inoreabbrev anvi-delete-statei {% include PROGRAM name="anvi-delete-state" %}
inoreabbrev anvi-dereplicate-genomesi {% include PROGRAM name="anvi-dereplicate-genomes" %}
inoreabbrev anvi-display-contigs-statsi {% include PROGRAM name="anvi-display-contigs-stats" %}
inoreabbrev anvi-display-functionsi {% include PROGRAM name="anvi-display-functions" %}
inoreabbrev anvi-display-metabolismi {% include PROGRAM name="anvi-display-metabolism" %}
inoreabbrev anvi-display-pani {% include PROGRAM name="anvi-display-pan" %}
inoreabbrev anvi-display-structurei {% include PROGRAM name="anvi-display-structure" %}
inoreabbrev anvi-estimate-genome-completenessi {% include PROGRAM name="anvi-estimate-genome-completeness" %}
inoreabbrev anvi-estimate-genome-taxonomyi {% include PROGRAM name="anvi-estimate-genome-taxonomy" %}
inoreabbrev anvi-estimate-metabolismi {% include PROGRAM name="anvi-estimate-metabolism" %}
inoreabbrev anvi-estimate-scg-taxonomyi {% include PROGRAM name="anvi-estimate-scg-taxonomy" %}
inoreabbrev anvi-estimate-trna-taxonomyi {% include PROGRAM name="anvi-estimate-trna-taxonomy" %}
inoreabbrev anvi-experimental-organizationi {% include PROGRAM name="anvi-experimental-organization" %}
inoreabbrev anvi-export-collectioni {% include PROGRAM name="anvi-export-collection" %}
inoreabbrev anvi-export-contigsi {% include PROGRAM name="anvi-export-contigs" %}
inoreabbrev anvi-export-functionsi {% include PROGRAM name="anvi-export-functions" %}
inoreabbrev anvi-export-gene-callsi {% include PROGRAM name="anvi-export-gene-calls" %}
inoreabbrev anvi-export-gene-coverage-and-detectioni {% include PROGRAM name="anvi-export-gene-coverage-and-detection" %}
inoreabbrev anvi-export-items-orderi {% include PROGRAM name="anvi-export-items-order" %}
inoreabbrev anvi-export-locusi {% include PROGRAM name="anvi-export-locus" %}
inoreabbrev anvi-export-misc-datai {% include PROGRAM name="anvi-export-misc-data" %}
inoreabbrev anvi-export-splits-and-coveragesi {% include PROGRAM name="anvi-export-splits-and-coverages" %}
inoreabbrev anvi-export-splits-taxonomyi {% include PROGRAM name="anvi-export-splits-taxonomy" %}
inoreabbrev anvi-export-statei {% include PROGRAM name="anvi-export-state" %}
inoreabbrev anvi-export-structuresi {% include PROGRAM name="anvi-export-structures" %}
inoreabbrev anvi-export-tablei {% include PROGRAM name="anvi-export-table" %}
inoreabbrev anvi-gen-contigs-databasei {% include PROGRAM name="anvi-gen-contigs-database" %}
inoreabbrev anvi-gen-fixation-index-matrixi {% include PROGRAM name="anvi-gen-fixation-index-matrix" %}
inoreabbrev anvi-gen-gene-consensus-sequencesi {% include PROGRAM name="anvi-gen-gene-consensus-sequences" %}
inoreabbrev anvi-gen-gene-level-stats-databasesi {% include PROGRAM name="anvi-gen-gene-level-stats-databases" %}
inoreabbrev anvi-gen-genomes-storagei {% include PROGRAM name="anvi-gen-genomes-storage" %}
inoreabbrev anvi-gen-networki {% include PROGRAM name="anvi-gen-network" %}
inoreabbrev anvi-gen-phylogenomic-treei {% include PROGRAM name="anvi-gen-phylogenomic-tree" %}
inoreabbrev anvi-gen-structure-databasei {% include PROGRAM name="anvi-gen-structure-database" %}
inoreabbrev anvi-gen-variability-matrixi {% include PROGRAM name="anvi-gen-variability-matrix" %}
inoreabbrev anvi-gen-variability-networki {% include PROGRAM name="anvi-gen-variability-network" %}
inoreabbrev anvi-gen-variability-profilei {% include PROGRAM name="anvi-gen-variability-profile" %}
inoreabbrev anvi-get-aa-countsi {% include PROGRAM name="anvi-get-aa-counts" %}
inoreabbrev anvi-get-codon-frequenciesi {% include PROGRAM name="anvi-get-codon-frequencies" %}
inoreabbrev anvi-get-pn-ps-ratioi {% include PROGRAM name="anvi-get-pn-ps-ratio" %}
inoreabbrev anvi-get-sequences-for-gene-callsi {% include PROGRAM name="anvi-get-sequences-for-gene-calls" %}
inoreabbrev anvi-get-sequences-for-gene-clustersi {% include PROGRAM name="anvi-get-sequences-for-gene-clusters" %}
inoreabbrev anvi-get-sequences-for-hmm-hitsi {% include PROGRAM name="anvi-get-sequences-for-hmm-hits" %}
inoreabbrev anvi-get-short-reads-from-bami {% include PROGRAM name="anvi-get-short-reads-from-bam" %}
inoreabbrev anvi-get-short-reads-mapping-to-a-genei {% include PROGRAM name="anvi-get-short-reads-mapping-to-a-gene" %}
inoreabbrev anvi-get-split-coveragesi {% include PROGRAM name="anvi-get-split-coverages" %}
inoreabbrev anvi-get-tlen-dist-from-bami {% include PROGRAM name="anvi-get-tlen-dist-from-bam" %}
inoreabbrev anvi-helpi {% include PROGRAM name="anvi-help" %}
inoreabbrev anvi-import-collectioni {% include PROGRAM name="anvi-import-collection" %}
inoreabbrev anvi-import-functionsi {% include PROGRAM name="anvi-import-functions" %}
inoreabbrev anvi-import-items-orderi {% include PROGRAM name="anvi-import-items-order" %}
inoreabbrev anvi-import-misc-datai {% include PROGRAM name="anvi-import-misc-data" %}
inoreabbrev anvi-import-statei {% include PROGRAM name="anvi-import-state" %}
inoreabbrev anvi-import-taxonomy-for-genesi {% include PROGRAM name="anvi-import-taxonomy-for-genes" %}
inoreabbrev anvi-import-taxonomy-for-layersi {% include PROGRAM name="anvi-import-taxonomy-for-layers" %}
inoreabbrev anvi-init-bami {% include PROGRAM name="anvi-init-bam" %}
inoreabbrev anvi-inspecti {% include PROGRAM name="anvi-inspect" %}
inoreabbrev anvi-interactivei {% include PROGRAM name="anvi-interactive" %}
inoreabbrev anvi-matrix-to-newicki {% include PROGRAM name="anvi-matrix-to-newick" %}
inoreabbrev anvi-mcg-classifieri {% include PROGRAM name="anvi-mcg-classifier" %}
inoreabbrev anvi-mergei {% include PROGRAM name="anvi-merge" %}
inoreabbrev anvi-merge-binsi {% include PROGRAM name="anvi-merge-bins" %}
inoreabbrev anvi-merge-trnaseqi {% include PROGRAM name="anvi-merge-trnaseq" %}
inoreabbrev anvi-meta-pan-genomei {% include PROGRAM name="anvi-meta-pan-genome" %}
inoreabbrev anvi-migratei {% include PROGRAM name="anvi-migrate" %}
inoreabbrev anvi-oligotype-linkmersi {% include PROGRAM name="anvi-oligotype-linkmers" %}
inoreabbrev anvi-pan-genomei {% include PROGRAM name="anvi-pan-genome" %}
inoreabbrev anvi-plot-trnaseqi {% include PROGRAM name="anvi-plot-trnaseq" %}
inoreabbrev anvi-profilei {% include PROGRAM name="anvi-profile" %}
inoreabbrev anvi-profile-blitzi {% include PROGRAM name="anvi-profile-blitz" %}
inoreabbrev anvi-pushi {% include PROGRAM name="anvi-push" %}
inoreabbrev anvi-refinei {% include PROGRAM name="anvi-refine" %}
inoreabbrev anvi-rename-binsi {% include PROGRAM name="anvi-rename-bins" %}
inoreabbrev anvi-report-inversionsi {% include PROGRAM name="anvi-report-inversions" %}
inoreabbrev anvi-report-linkmersi {% include PROGRAM name="anvi-report-linkmers" %}
inoreabbrev anvi-run-hmmsi {% include PROGRAM name="anvi-run-hmms" %}
inoreabbrev anvi-run-interacdomei {% include PROGRAM name="anvi-run-interacdome" %}
inoreabbrev anvi-run-kegg-kofamsi {% include PROGRAM name="anvi-run-kegg-kofams" %}
inoreabbrev anvi-run-ncbi-cogsi {% include PROGRAM name="anvi-run-ncbi-cogs" %}
inoreabbrev anvi-run-pfamsi {% include PROGRAM name="anvi-run-pfams" %}
inoreabbrev anvi-run-scg-taxonomyi {% include PROGRAM name="anvi-run-scg-taxonomy" %}
inoreabbrev anvi-run-trna-taxonomyi {% include PROGRAM name="anvi-run-trna-taxonomy" %}
inoreabbrev anvi-run-workflowi {% include PROGRAM name="anvi-run-workflow" %}
inoreabbrev anvi-scan-trnasi {% include PROGRAM name="anvi-scan-trnas" %}
inoreabbrev anvi-search-functionsi {% include PROGRAM name="anvi-search-functions" %}
inoreabbrev anvi-search-palindromesi {% include PROGRAM name="anvi-search-palindromes" %}
inoreabbrev anvi-search-sequence-motifsi {% include PROGRAM name="anvi-search-sequence-motifs" %}
inoreabbrev anvi-self-testi {% include PROGRAM name="anvi-self-test" %}
inoreabbrev anvi-setup-interacdomei {% include PROGRAM name="anvi-setup-interacdome" %}
inoreabbrev anvi-setup-kegg-kofamsi {% include PROGRAM name="anvi-setup-kegg-kofams" %}
inoreabbrev anvi-setup-ncbi-cogsi {% include PROGRAM name="anvi-setup-ncbi-cogs" %}
inoreabbrev anvi-setup-pdb-databasei {% include PROGRAM name="anvi-setup-pdb-database" %}
inoreabbrev anvi-setup-pfamsi {% include PROGRAM name="anvi-setup-pfams" %}
inoreabbrev anvi-setup-scg-taxonomyi {% include PROGRAM name="anvi-setup-scg-taxonomy" %}
inoreabbrev anvi-setup-trna-taxonomyi {% include PROGRAM name="anvi-setup-trna-taxonomy" %}
inoreabbrev anvi-show-collections-and-binsi {% include PROGRAM name="anvi-show-collections-and-bins" %}
inoreabbrev anvi-show-misc-datai {% include PROGRAM name="anvi-show-misc-data" %}
inoreabbrev anvi-spliti {% include PROGRAM name="anvi-split" %}
inoreabbrev anvi-summarizei {% include PROGRAM name="anvi-summarize" %}
inoreabbrev anvi-summarize-blitzi {% include PROGRAM name="anvi-summarize-blitz" %}
inoreabbrev anvi-tabulate-trnaseqi {% include PROGRAM name="anvi-tabulate-trnaseq" %}
inoreabbrev anvi-trnaseqi {% include PROGRAM name="anvi-trnaseq" %}
inoreabbrev anvi-update-db-descriptioni {% include PROGRAM name="anvi-update-db-description" %}
inoreabbrev anvi-update-structure-databasei {% include PROGRAM name="anvi-update-structure-database" %}
inoreabbrev anvi-upgradei {% include PROGRAM name="anvi-upgrade" %}
inoreabbrev anvi-script-add-default-collectioni {% include PROGRAM name="anvi-script-add-default-collection" %}
inoreabbrev anvi-script-as-markdowni {% include PROGRAM name="anvi-script-as-markdown" %}
inoreabbrev anvi-script-checkm-tree-to-interactivei {% include PROGRAM name="anvi-script-checkm-tree-to-interactive" %}
inoreabbrev anvi-script-compute-ani-for-fastai {% include PROGRAM name="anvi-script-compute-ani-for-fasta" %}
inoreabbrev anvi-script-compute-bayesian-pan-corei {% include PROGRAM name="anvi-script-compute-bayesian-pan-core" %}
inoreabbrev anvi-script-enrichment-statsi {% include PROGRAM name="anvi-script-enrichment-stats" %}
inoreabbrev anvi-script-estimate-genome-sizei {% include PROGRAM name="anvi-script-estimate-genome-size" %}
inoreabbrev anvi-script-filter-fasta-by-blasti {% include PROGRAM name="anvi-script-filter-fasta-by-blast" %}
inoreabbrev anvi-script-filter-hmm-hits-tablei {% include PROGRAM name="anvi-script-filter-hmm-hits-table" %}
inoreabbrev anvi-script-fix-homopolymer-indelsi {% include PROGRAM name="anvi-script-fix-homopolymer-indels" %}
inoreabbrev anvi-script-gen-CPR-classifieri {% include PROGRAM name="anvi-script-gen-CPR-classifier" %}
inoreabbrev anvi-script-gen-distribution-of-genes-in-a-bini {% include PROGRAM name="anvi-script-gen-distribution-of-genes-in-a-bin" %}
inoreabbrev anvi-script-gen-functions-per-group-stats-outputi {% include PROGRAM name="anvi-script-gen-functions-per-group-stats-output" %}
inoreabbrev anvi-script-gen-genomes-filei {% include PROGRAM name="anvi-script-gen-genomes-file" %}
inoreabbrev anvi-script-gen-help-pagesi {% include PROGRAM name="anvi-script-gen-help-pages" %}
inoreabbrev anvi-script-gen-hmm-hits-matrix-across-genomesi {% include PROGRAM name="anvi-script-gen-hmm-hits-matrix-across-genomes" %}
inoreabbrev anvi-script-gen-programs-networki {% include PROGRAM name="anvi-script-gen-programs-network" %}
inoreabbrev anvi-script-gen-programs-vignettei {% include PROGRAM name="anvi-script-gen-programs-vignette" %}
inoreabbrev anvi-script-gen-pseudo-paired-reads-from-fastqi {% include PROGRAM name="anvi-script-gen-pseudo-paired-reads-from-fastq" %}
inoreabbrev anvi-script-gen-scg-domain-classifieri {% include PROGRAM name="anvi-script-gen-scg-domain-classifier" %}
inoreabbrev anvi-script-gen-short-readsi {% include PROGRAM name="anvi-script-gen-short-reads" %}
inoreabbrev anvi-script-get-collection-infoi {% include PROGRAM name="anvi-script-get-collection-info" %}
inoreabbrev anvi-script-get-coverage-from-bami {% include PROGRAM name="anvi-script-get-coverage-from-bam" %}
inoreabbrev anvi-script-get-hmm-hits-per-gene-calli {% include PROGRAM name="anvi-script-get-hmm-hits-per-gene-call" %}
inoreabbrev anvi-script-get-primer-matchesi {% include PROGRAM name="anvi-script-get-primer-matches" %}
inoreabbrev anvi-script-merge-collectionsi {% include PROGRAM name="anvi-script-merge-collections" %}
inoreabbrev anvi-script-permute-trnaseq-seedsi {% include PROGRAM name="anvi-script-permute-trnaseq-seeds" %}
inoreabbrev anvi-script-pfam-accessions-to-hmms-directoryi {% include PROGRAM name="anvi-script-pfam-accessions-to-hmms-directory" %}
inoreabbrev anvi-script-predict-CPR-genomesi {% include PROGRAM name="anvi-script-predict-CPR-genomes" %}
inoreabbrev anvi-script-process-genbanki {% include PROGRAM name="anvi-script-process-genbank" %}
inoreabbrev anvi-script-process-genbank-metadatai {% include PROGRAM name="anvi-script-process-genbank-metadata" %}
inoreabbrev anvi-script-reformat-fastai {% include PROGRAM name="anvi-script-reformat-fasta" %}
inoreabbrev anvi-script-run-eggnog-mapperi {% include PROGRAM name="anvi-script-run-eggnog-mapper" %}
inoreabbrev anvi-script-snvs-to-interactivei {% include PROGRAM name="anvi-script-snvs-to-interactive" %}
inoreabbrev anvi-script-tabulatei {% include PROGRAM name="anvi-script-tabulate" %}
inoreabbrev anvi-script-transpose-matrixi {% include PROGRAM name="anvi-script-transpose-matrix" %}
inoreabbrev anvi-script-variability-to-vcfi {% include PROGRAM name="anvi-script-variability-to-vcf" %}
inoreabbrev anvi-script-visualize-split-coveragesi {% include PROGRAM name="anvi-script-visualize-split-coverages" %}
{% endraw %}
```

Just for posterity, this is how I generated the data above if you ever wish to update it.

For artifacts, I run this (assuming that you have the anvi'o codebase at `~/github/anvio/`):

``` bash
{% raw %}
for i in `grep '{' ~/github/anvio/anvio/docs/__init__.py | grep '"' | awk 'BEGIN{FS="\""}{print $2}'`
do
    echo inoreabbrev ${i}i {% include ARTIFACT name=\"$i\" %}
done
{% endraw%}
```

For programs, I run these two:

``` bash
{% raw %}
for i in `ls ~/github/anvio/bin/anvi-* ~/github/anvio/sandbox/anvi-* | grep -v '_' | grep -v "augustus"`
do
    echo inoreabbrev $(basename $i)i {% include PROGRAM name=\"$(basename $i)\" %}
done
{% endraw%}
```

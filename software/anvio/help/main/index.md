---
layout: page
title: Help pages for anvi'o programs and artifacts
categories: [anvio]
comments: false
image:
  featurerelative: images/header.png
  display: true
redirect_from:
    - /help
---

Here you will find a list of all anvi'o programs and artifacts that enable constructing workflows for integrated multi 'omics investigations.

If you need an introduction to the terminology used in 'omics research or in anvi'o, please take a look at <a href="http://merenlab.org/vocabulary/">our vocabulary page</a>. The anvi'o community is with you! If you have practical, technical, or science questions <a href="http://merenlab.org/2019/10/07/getting-help/">this page</a> to learn about resources available to you. If you are feeling overwhelmed, you can always scream towards the anvi'o Slack channel.


{% include _join-anvio-slack.html %}


<a href="/software/anvio/network/" target="_blank"><img src="/images/anvio-network.png" width="100%" /></a>

{:.notice}
The help contents were last updated on **12 Apr 21 20:21:26** for anvi'o version **7-dev (hope)**.


{% include _project-anvio-version.html %}
{% include _toc.html %}



## Anvi'o artifacts

Anvi'o artifacts represent **concepts, file types, or data types** anvi'o programs can work with. A given anvi'o artifact can be provided by the user (such as a FASTA file), produced by anvi'o (such as a profile database), or both (such as phylogenomic trees). Anvi'o artifacts link anvi'o programs to each other to build novel workflows.

Listed below **a total of 105 artifacts**.

<table class="artifacts-table">

<tr style="border:none;">
<td><div class="artifact-icon-div"><img src="images/icons/DB.png" class="artifact-icon" /></div></td><td class="artifact-list-td"><span class="artifact-n" markdown="1">[pan-db](artifacts/pan-db)</span> <span class="artifact-n" markdown="1">[contigs-db](artifacts/contigs-db)</span> <span class="artifact-n" markdown="1">[trnaseq-db](artifacts/trnaseq-db)</span> <span class="artifact-n" markdown="1">[modules-db](artifacts/modules-db)</span> <span class="artifact-n" markdown="1">[structure-db](artifacts/structure-db)</span> <span class="artifact-n" markdown="1">[pdb-db](artifacts/pdb-db)</span> <span class="artifact-n" markdown="1">[kegg-data](artifacts/kegg-data)</span> <span class="artifact-n" markdown="1">[single-profile-db](artifacts/single-profile-db)</span> <span class="artifact-n" markdown="1">[profile-db](artifacts/profile-db)</span> <span class="artifact-n" markdown="1">[genes-db](artifacts/genes-db)</span> <span class="artifact-n" markdown="1">[genomes-storage-db](artifacts/genomes-storage-db)</span></td>
</tr>

<tr>
<td><div class="artifact-icon-div"><img src="images/icons/FASTA.png" class="artifact-icon" /></div></td><td class="artifact-list-td"><span class="artifact-n" markdown="1">[fasta](artifacts/fasta)</span> <span class="artifact-n" markdown="1">[contigs-fasta](artifacts/contigs-fasta)</span> <span class="artifact-n" markdown="1">[trnaseq-fasta](artifacts/trnaseq-fasta)</span> <span class="artifact-n" markdown="1">[concatenated-gene-alignment-fasta](artifacts/concatenated-gene-alignment-fasta)</span> <span class="artifact-n" markdown="1">[short-reads-fasta](artifacts/short-reads-fasta)</span> <span class="artifact-n" markdown="1">[genes-fasta](artifacts/genes-fasta)</span> <span class="artifact-n" markdown="1">[locus-fasta](artifacts/locus-fasta)</span></td>
</tr>

<tr>
<td><div class="artifact-icon-div"><img src="images/icons/TXT.png" class="artifact-icon" /></div></td><td class="artifact-list-td"><span class="artifact-n" markdown="1">[configuration-ini](artifacts/configuration-ini)</span> <span class="artifact-n" markdown="1">[external-gene-calls](artifacts/external-gene-calls)</span> <span class="artifact-n" markdown="1">[protein-structure-txt](artifacts/protein-structure-txt)</span> <span class="artifact-n" markdown="1">[samples-txt](artifacts/samples-txt)</span> <span class="artifact-n" markdown="1">[fasta-txt](artifacts/fasta-txt)</span> <span class="artifact-n" markdown="1">[collection-txt](artifacts/collection-txt)</span> <span class="artifact-n" markdown="1">[misc-data-items-txt](artifacts/misc-data-items-txt)</span> <span class="artifact-n" markdown="1">[misc-data-layers-txt](artifacts/misc-data-layers-txt)</span> <span class="artifact-n" markdown="1">[misc-data-nucleotides-txt](artifacts/misc-data-nucleotides-txt)</span> <span class="artifact-n" markdown="1">[misc-data-amino-acids-txt](artifacts/misc-data-amino-acids-txt)</span> <span class="artifact-n" markdown="1">[misc-data-layer-orders-txt](artifacts/misc-data-layer-orders-txt)</span> <span class="artifact-n" markdown="1">[misc-data-items-order-txt](artifacts/misc-data-items-order-txt)</span> <span class="artifact-n" markdown="1">[linkmers-txt](artifacts/linkmers-txt)</span> <span class="artifact-n" markdown="1">[gene-calls-txt](artifacts/gene-calls-txt)</span> <span class="artifact-n" markdown="1">[binding-frequencies-txt](artifacts/binding-frequencies-txt)</span> <span class="artifact-n" markdown="1">[functions-txt](artifacts/functions-txt)</span> <span class="artifact-n" markdown="1">[functional-enrichment-txt](artifacts/functional-enrichment-txt)</span> <span class="artifact-n" markdown="1">[view-data](artifacts/view-data)</span> <span class="artifact-n" markdown="1">[layer-taxonomy-txt](artifacts/layer-taxonomy-txt)</span> <span class="artifact-n" markdown="1">[gene-taxonomy-txt](artifacts/gene-taxonomy-txt)</span> <span class="artifact-n" markdown="1">[genome-taxonomy-txt](artifacts/genome-taxonomy-txt)</span> <span class="artifact-n" markdown="1">[external-genomes](artifacts/external-genomes)</span> <span class="artifact-n" markdown="1">[internal-genomes](artifacts/internal-genomes)</span> <span class="artifact-n" markdown="1">[metagenomes](artifacts/metagenomes)</span> <span class="artifact-n" markdown="1">[coverages-txt](artifacts/coverages-txt)</span> <span class="artifact-n" markdown="1">[detection-txt](artifacts/detection-txt)</span> <span class="artifact-n" markdown="1">[variability-profile-txt](artifacts/variability-profile-txt)</span> <span class="artifact-n" markdown="1">[codon-frequencies-txt](artifacts/codon-frequencies-txt)</span> <span class="artifact-n" markdown="1">[aa-frequencies-txt](artifacts/aa-frequencies-txt)</span> <span class="artifact-n" markdown="1">[fixation-index-matrix](artifacts/fixation-index-matrix)</span> <span class="artifact-n" markdown="1">[kegg-metabolism](artifacts/kegg-metabolism)</span> <span class="artifact-n" markdown="1">[augustus-gene-calls](artifacts/augustus-gene-calls)</span> <span class="artifact-n" markdown="1">[vcf](artifacts/vcf)</span> <span class="artifact-n" markdown="1">[blast-table](artifacts/blast-table)</span> <span class="artifact-n" markdown="1">[splits-txt](artifacts/splits-txt)</span> <span class="artifact-n" markdown="1">[genbank-file](artifacts/genbank-file)</span> <span class="artifact-n" markdown="1">[groups-txt](artifacts/groups-txt)</span> <span class="artifact-n" markdown="1">[splits-taxonomy-txt](artifacts/splits-taxonomy-txt)</span> <span class="artifact-n" markdown="1">[hmm-hits-matrix-txt](artifacts/hmm-hits-matrix-txt)</span> <span class="artifact-n" markdown="1">[clustering-configuration](artifacts/clustering-configuration)</span></td>
</tr>

<tr>
<td><div class="artifact-icon-div"><img src="images/icons/BAM.png" class="artifact-icon" /></div></td><td class="artifact-list-td"><span class="artifact-n" markdown="1">[bam-file](artifacts/bam-file)</span> <span class="artifact-n" markdown="1">[raw-bam-file](artifacts/raw-bam-file)</span></td>
</tr>

<tr>
<td><div class="artifact-icon-div"><img src="images/icons/STATS.png" class="artifact-icon" /></div></td><td class="artifact-list-td"><span class="artifact-n" markdown="1">[contigs-stats](artifacts/contigs-stats)</span> <span class="artifact-n" markdown="1">[genes-stats](artifacts/genes-stats)</span></td>
</tr>

<tr>
<td><div class="artifact-icon-div"><img src="images/icons/SVG.png" class="artifact-icon" /></div></td><td class="artifact-list-td"><span class="artifact-n" markdown="1">[svg](artifacts/svg)</span></td>
</tr>

<tr>
<td><div class="artifact-icon-div"><img src="images/icons/BIN.png" class="artifact-icon" /></div></td><td class="artifact-list-td"><span class="artifact-n" markdown="1">[bin](artifacts/bin)</span></td>
</tr>

<tr>
<td><div class="artifact-icon-div"><img src="images/icons/COLLECTION.png" class="artifact-icon" /></div></td><td class="artifact-list-td"><span class="artifact-n" markdown="1">[collection](artifacts/collection)</span></td>
</tr>

<tr>
<td><div class="artifact-icon-div"><img src="images/icons/HMM.png" class="artifact-icon" /></div></td><td class="artifact-list-td"><span class="artifact-n" markdown="1">[hmm-source](artifacts/hmm-source)</span></td>
</tr>

<tr>
<td><div class="artifact-icon-div"><img src="images/icons/CONCEPT.png" class="artifact-icon" /></div></td><td class="artifact-list-td"><span class="artifact-n" markdown="1">[hmm-hits](artifacts/hmm-hits)</span> <span class="artifact-n" markdown="1">[completion](artifacts/completion)</span> <span class="artifact-n" markdown="1">[misc-data-items](artifacts/misc-data-items)</span> <span class="artifact-n" markdown="1">[misc-data-layers](artifacts/misc-data-layers)</span> <span class="artifact-n" markdown="1">[misc-data-nucleotides](artifacts/misc-data-nucleotides)</span> <span class="artifact-n" markdown="1">[misc-data-amino-acids](artifacts/misc-data-amino-acids)</span> <span class="artifact-n" markdown="1">[genome-similarity](artifacts/genome-similarity)</span> <span class="artifact-n" markdown="1">[misc-data-layer-orders](artifacts/misc-data-layer-orders)</span> <span class="artifact-n" markdown="1">[misc-data-items-order](artifacts/misc-data-items-order)</span> <span class="artifact-n" markdown="1">[metapangenome](artifacts/metapangenome)</span> <span class="artifact-n" markdown="1">[oligotypes](artifacts/oligotypes)</span> <span class="artifact-n" markdown="1">[functions](artifacts/functions)</span> <span class="artifact-n" markdown="1">[kegg-functions](artifacts/kegg-functions)</span> <span class="artifact-n" markdown="1">[layer-taxonomy](artifacts/layer-taxonomy)</span> <span class="artifact-n" markdown="1">[gene-taxonomy](artifacts/gene-taxonomy)</span> <span class="artifact-n" markdown="1">[genome-taxonomy](artifacts/genome-taxonomy)</span> <span class="artifact-n" markdown="1">[scgs-taxonomy-db](artifacts/scgs-taxonomy-db)</span> <span class="artifact-n" markdown="1">[scgs-taxonomy](artifacts/scgs-taxonomy)</span> <span class="artifact-n" markdown="1">[trna-taxonomy-db](artifacts/trna-taxonomy-db)</span> <span class="artifact-n" markdown="1">[trna-taxonomy](artifacts/trna-taxonomy)</span> <span class="artifact-n" markdown="1">[variability-profile](artifacts/variability-profile)</span> <span class="artifact-n" markdown="1">[split-bins](artifacts/split-bins)</span> <span class="artifact-n" markdown="1">[state](artifacts/state)</span> <span class="artifact-n" markdown="1">[ngrams](artifacts/ngrams)</span> <span class="artifact-n" markdown="1">[pn-ps-data](artifacts/pn-ps-data)</span></td>
</tr>

<tr>
<td><div class="artifact-icon-div"><img src="images/icons/DATA.png" class="artifact-icon" /></div></td><td class="artifact-list-td"><span class="artifact-n" markdown="1">[cogs-data](artifacts/cogs-data)</span> <span class="artifact-n" markdown="1">[pfams-data](artifacts/pfams-data)</span> <span class="artifact-n" markdown="1">[interacdome-data](artifacts/interacdome-data)</span></td>
</tr>

<tr>
<td><div class="artifact-icon-div"><img src="images/icons/NEWICK.png" class="artifact-icon" /></div></td><td class="artifact-list-td"><span class="artifact-n" markdown="1">[dendrogram](artifacts/dendrogram)</span> <span class="artifact-n" markdown="1">[phylogeny](artifacts/phylogeny)</span></td>
</tr>

<tr>
<td><div class="artifact-icon-div"><img src="images/icons/DISPLAY.png" class="artifact-icon" /></div></td><td class="artifact-list-td"><span class="artifact-n" markdown="1">[interactive](artifacts/interactive)</span></td>
</tr>

<tr>
<td><div class="artifact-icon-div"><img src="images/icons/SUMMARY.png" class="artifact-icon" /></div></td><td class="artifact-list-td"><span class="artifact-n" markdown="1">[summary](artifacts/summary)</span></td>
</tr>

<tr>
<td><div class="artifact-icon-div"><img src="images/icons/JSON.png" class="artifact-icon" /></div></td><td class="artifact-list-td"><span class="artifact-n" markdown="1">[state-json](artifacts/state-json)</span> <span class="artifact-n" markdown="1">[workflow-config](artifacts/workflow-config)</span></td>
</tr>

<tr>
<td><div class="artifact-icon-div"><img src="images/icons/WORKFLOW.png" class="artifact-icon" /></div></td><td class="artifact-list-td"><span class="artifact-n" markdown="1">[contigs-workflow](artifacts/contigs-workflow)</span> <span class="artifact-n" markdown="1">[metagenomics-workflow](artifacts/metagenomics-workflow)</span> <span class="artifact-n" markdown="1">[pangenomics-workflow](artifacts/pangenomics-workflow)</span> <span class="artifact-n" markdown="1">[phylogenomics-workflow](artifacts/phylogenomics-workflow)</span> <span class="artifact-n" markdown="1">[trnaseq-workflow](artifacts/trnaseq-workflow)</span></td>
</tr>

</table>

## Anvi'o programs

Anvi'o programs perform atomic tasks that can be weaved together to implement complete 'omics workflows. Please note that there may be programs that are not listed on this page. You can type 'anvi-' in your terminal, and press the TAB key twice to see the full list of programs available to you on your system, and type `anvi-program-name --help` to read the full list of command line options.

Listed below **a total of 125 programs**.


<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-analyze-synteny](programs/anvi-analyze-synteny)**</span>. <span markdown="1">Extract ngrams, as in &#x27;co-occurring genes in synteny&#x27;, from genomes</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[genomes-storage-db](artifacts/genomes-storage-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[functions](artifacts/functions) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[pan-db](artifacts/pan-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[ngrams](artifacts/ngrams) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-cluster-contigs](programs/anvi-cluster-contigs)**</span>. <span markdown="1">A program to cluster items in a merged anvi&#x27;o profile using automatic binning algorithms</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[collection](artifacts/collection) <img src="images/icons/COLLECTION.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[collection](artifacts/collection) <img src="images/icons/COLLECTION.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[bin](artifacts/bin) <img src="images/icons/BIN.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-compute-completeness](programs/anvi-compute-completeness)**</span>. <span markdown="1">A script to generate completeness info for a given list of _splits_</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[splits-txt](artifacts/splits-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[hmm-source](artifacts/hmm-source) <img src="images/icons/HMM.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> 
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-compute-functional-enrichment](programs/anvi-compute-functional-enrichment)**</span>. <span markdown="1">This is a driver program for `anvi-script-enrichment-stats`, a script that computes enrichment scores and group associations for annotated entities (ie, functions, KEGG Modules) across groups of genomes or samples.</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[kegg-metabolism](artifacts/kegg-metabolism) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[groups-txt](artifacts/groups-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[misc-data-layers](artifacts/misc-data-layers) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[pan-db](artifacts/pan-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[genomes-storage-db](artifacts/genomes-storage-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[external-genomes](artifacts/external-genomes) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[internal-genomes](artifacts/internal-genomes) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[functions](artifacts/functions) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[functional-enrichment-txt](artifacts/functional-enrichment-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-compute-gene-cluster-homogeneity](programs/anvi-compute-gene-cluster-homogeneity)**</span>. <span markdown="1">Compute homogeneity for gene clusters</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[pan-db](artifacts/pan-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[genomes-storage-db](artifacts/genomes-storage-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> 
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-compute-genome-similarity](programs/anvi-compute-genome-similarity)**</span>. <span markdown="1">Export sequences from sequence sources and compute a similarity metric (e.g. ANI). If a Pan Database is given anvi&#x27;o will write computed output to misc data tables of Pan Database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[external-genomes](artifacts/external-genomes) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[internal-genomes](artifacts/internal-genomes) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[pan-db](artifacts/pan-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[genome-similarity](artifacts/genome-similarity) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-convert-trnaseq-database](programs/anvi-convert-trnaseq-database)**</span>. <span markdown="1">A program that processes one or more anvio&#x27; tRNA-seq databases generated by `anvi-trnaseq` to generate anvi&#x27;o contigs and merged profile databases that are accessible to the rest of the tools in anvi&#x27;o software ecosystem. Briefly, this program will determine final seed sequences from input tRNA-seq databases, determine their coverages across samples, identify tRNA modification sites, and INDELs associated with transcripts in each sample against the seed sequences and store all these data into resulting databases for interactive visualization of the data or in-depth analysis using other anvi&#x27;o frameworks</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[trnaseq-db](artifacts/trnaseq-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-db-info](programs/anvi-db-info)**</span>. <span markdown="1">Access self tables, display values, or set new ones totally on your own risk</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[pan-db](artifacts/pan-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[genomes-storage-db](artifacts/genomes-storage-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[structure-db](artifacts/structure-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[genes-db](artifacts/genes-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> 
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-delete-collection](programs/anvi-delete-collection)**</span>. <span markdown="1">Remove a collection from a given profile database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[collection](artifacts/collection) <img src="images/icons/COLLECTION.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> 
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-delete-functions](programs/anvi-delete-functions)**</span>. <span markdown="1">Remove functional annotation sources from an anvi&#x27;o contigs database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[functions](artifacts/functions) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> 
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-delete-hmms](programs/anvi-delete-hmms)**</span>. <span markdown="1">Remove HMM hits from an anvi&#x27;o contigs database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[hmm-source](artifacts/hmm-source) <img src="images/icons/HMM.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[hmm-hits](artifacts/hmm-hits) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> 
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-delete-misc-data](programs/anvi-delete-misc-data)**</span>. <span markdown="1">Remove stuff from &#x27;additional data&#x27; or &#x27;order&#x27; tables for either items or layers in either pan or profile databases. OR, remove stuff from the &#x27;additional data&#x27; tables for nucleotides or amino acids in contigs databases</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[pan-db](artifacts/pan-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[misc-data-items](artifacts/misc-data-items) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[misc-data-layers](artifacts/misc-data-layers) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[misc-data-layer-orders](artifacts/misc-data-layer-orders) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[misc-data-nucleotides](artifacts/misc-data-nucleotides) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[misc-data-amino-acids](artifacts/misc-data-amino-acids) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> 
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-delete-state](programs/anvi-delete-state)**</span>. <span markdown="1">Delete an anvi&#x27;o state from a pan or profile database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[pan-db](artifacts/pan-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[state](artifacts/state) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> 
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-dereplicate-genomes](programs/anvi-dereplicate-genomes)**</span>. <span markdown="1">Identify redundant (highly similar) genomes</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[external-genomes](artifacts/external-genomes) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[internal-genomes](artifacts/internal-genomes) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[fasta](artifacts/fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[genome-similarity](artifacts/genome-similarity) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[fasta](artifacts/fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-display-contigs-stats](programs/anvi-display-contigs-stats)**</span>. <span markdown="1">Start the anvi&#x27;o interactive interactive for viewing or comparing contigs statistics</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[contigs-stats](artifacts/contigs-stats) <img src="images/icons/STATS.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[interactive](artifacts/interactive) <img src="images/icons/DISPLAY.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[svg](artifacts/svg) <img src="images/icons/SVG.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-display-functions](programs/anvi-display-functions)**</span>. <span markdown="1">Start an anvi&#x27;o interactive display to see functions across genomes</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[functions](artifacts/functions) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[genomes-storage-db](artifacts/genomes-storage-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[internal-genomes](artifacts/internal-genomes) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[external-genomes](artifacts/external-genomes) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[groups-txt](artifacts/groups-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[interactive](artifacts/interactive) <img src="images/icons/DISPLAY.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[functional-enrichment-txt](artifacts/functional-enrichment-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-display-metabolism](programs/anvi-display-metabolism)**</span>. <span markdown="1">Start the anvi&#x27;o interactive interactive for viewing KEGG metabolism data</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[kegg-data](artifacts/kegg-data) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[kegg-functions](artifacts/kegg-functions) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[collection](artifacts/collection) <img src="images/icons/COLLECTION.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[bin](artifacts/bin) <img src="images/icons/BIN.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[interactive](artifacts/interactive) <img src="images/icons/DISPLAY.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-display-pan](programs/anvi-display-pan)**</span>. <span markdown="1">Start an anvi&#x27;o server to display a pan-genome</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[pan-db](artifacts/pan-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[genomes-storage-db](artifacts/genomes-storage-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[collection](artifacts/collection) <img src="images/icons/COLLECTION.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[bin](artifacts/bin) <img src="images/icons/BIN.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[interactive](artifacts/interactive) <img src="images/icons/DISPLAY.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[svg](artifacts/svg) <img src="images/icons/SVG.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-display-structure](programs/anvi-display-structure)**</span>. <span markdown="1">Interactively visualize sequence variants on protein structures</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[structure-db](artifacts/structure-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[variability-profile-txt](artifacts/variability-profile-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[splits-txt](artifacts/splits-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[interactive](artifacts/interactive) <img src="images/icons/DISPLAY.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-estimate-genome-completeness](programs/anvi-estimate-genome-completeness)**</span>. <span markdown="1">Estimate completion and redundancy using domain-specific single-copy core genes</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[external-genomes](artifacts/external-genomes) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[collection](artifacts/collection) <img src="images/icons/COLLECTION.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[completion](artifacts/completion) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-estimate-metabolism](programs/anvi-estimate-metabolism)**</span>. <span markdown="1">Reconstructs metabolic pathways and estimates pathway completeness for a given set of contigs</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[kegg-data](artifacts/kegg-data) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[kegg-functions](artifacts/kegg-functions) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[collection](artifacts/collection) <img src="images/icons/COLLECTION.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[bin](artifacts/bin) <img src="images/icons/BIN.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[external-genomes](artifacts/external-genomes) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[internal-genomes](artifacts/internal-genomes) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[metagenomes](artifacts/metagenomes) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[kegg-metabolism](artifacts/kegg-metabolism) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-estimate-scg-taxonomy](programs/anvi-estimate-scg-taxonomy)**</span>. <span markdown="1">Estimates taxonomy at genome and metagenome level. This program is the entry point to estimate taxonomy for a given set of contigs (i.e., all contigs in a contigs database, or contigs described in collections as bins). For this, it uses single-copy core gene sequences and the GTDB database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[scgs-taxonomy](artifacts/scgs-taxonomy) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[collection](artifacts/collection) <img src="images/icons/COLLECTION.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[bin](artifacts/bin) <img src="images/icons/BIN.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[metagenomes](artifacts/metagenomes) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[genome-taxonomy](artifacts/genome-taxonomy) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[genome-taxonomy-txt](artifacts/genome-taxonomy-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-estimate-trna-taxonomy](programs/anvi-estimate-trna-taxonomy)**</span>. <span markdown="1">Estimates taxonomy at genome and metagenome level using tRNA sequences.</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[trna-taxonomy](artifacts/trna-taxonomy) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[collection](artifacts/collection) <img src="images/icons/COLLECTION.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[bin](artifacts/bin) <img src="images/icons/BIN.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[metagenomes](artifacts/metagenomes) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[genome-taxonomy](artifacts/genome-taxonomy) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[genome-taxonomy-txt](artifacts/genome-taxonomy-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-experimental-organization](programs/anvi-experimental-organization)**</span>. <span markdown="1">Create an experimental clustering dendrogram.</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[clustering-configuration](artifacts/clustering-configuration) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[dendrogram](artifacts/dendrogram) <img src="images/icons/NEWICK.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-export-collection](programs/anvi-export-collection)**</span>. <span markdown="1">Export a collection from an anvi&#x27;o database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[collection](artifacts/collection) <img src="images/icons/COLLECTION.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[collection-txt](artifacts/collection-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-export-contigs](programs/anvi-export-contigs)**</span>. <span markdown="1">Export contigs (or splits) from an anvi&#x27;o contigs database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[contigs-fasta](artifacts/contigs-fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-export-functions](programs/anvi-export-functions)**</span>. <span markdown="1">Export functions of genes from an anvi&#x27;o contigs database for a given annotation source</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[functions](artifacts/functions) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[functions-txt](artifacts/functions-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-export-gene-calls](programs/anvi-export-gene-calls)**</span>. <span markdown="1">Export gene calls from an anvi&#x27;o contigs database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[gene-calls-txt](artifacts/gene-calls-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-export-gene-coverage-and-detection](programs/anvi-export-gene-coverage-and-detection)**</span>. <span markdown="1">Export gene coverage and detection data for all genes associated with contigs described in a profile database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[coverages-txt](artifacts/coverages-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[detection-txt](artifacts/detection-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-export-items-order](programs/anvi-export-items-order)**</span>. <span markdown="1">Export an item order from an anvi&#x27;o database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[pan-db](artifacts/pan-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[misc-data-items-order-txt](artifacts/misc-data-items-order-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[dendrogram](artifacts/dendrogram) <img src="images/icons/NEWICK.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[phylogeny](artifacts/phylogeny) <img src="images/icons/NEWICK.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-export-locus](programs/anvi-export-locus)**</span>. <span markdown="1">This program helps you cut a &#x27;locus&#x27; from a larger genetic context (e.g., contigs, genomes). By default, anvi&#x27;o will locate a user-defined anchor gene, extend its selection upstream and downstream based on the --num-genes argument, then extract the locus to create a new contigs database. The anchor gene must be provided as --search-term, --gene-caller-ids, or --hmm-sources. If --flank-mode is designated, you MUST provide TWO flanking genes that define the locus region (Please see --flank-mode help for more information). If everything goes as plan, anvi&#x27;o will give you individual locus contigs databases for every matching anchor gene found in the original contigs database provided. Enjoy your mini contigs databases!</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[locus-fasta](artifacts/locus-fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-export-misc-data](programs/anvi-export-misc-data)**</span>. <span markdown="1">Export additional data or order tables in pan or profile databases for items or layers</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[pan-db](artifacts/pan-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[misc-data-items](artifacts/misc-data-items) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[misc-data-layers](artifacts/misc-data-layers) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[misc-data-layer-orders](artifacts/misc-data-layer-orders) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[misc-data-nucleotides](artifacts/misc-data-nucleotides) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[misc-data-amino-acids](artifacts/misc-data-amino-acids) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[misc-data-items-txt](artifacts/misc-data-items-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[misc-data-layers-txt](artifacts/misc-data-layers-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[misc-data-layer-orders-txt](artifacts/misc-data-layer-orders-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[misc-data-nucleotides-txt](artifacts/misc-data-nucleotides-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[misc-data-amino-acids-txt](artifacts/misc-data-amino-acids-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-export-splits-and-coverages](programs/anvi-export-splits-and-coverages)**</span>. <span markdown="1">Export split or contig sequences and coverages across samples stored in an anvi&#x27;o profile database. This program is especially useful if you would like to &#x27;bin&#x27; your splits or contigs outside of anvi&#x27;o and import the binning results into anvi&#x27;o using `anvi-import-collection` program</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[contigs-fasta](artifacts/contigs-fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[coverages-txt](artifacts/coverages-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-export-splits-taxonomy](programs/anvi-export-splits-taxonomy)**</span>. <span markdown="1">Export taxonomy for splits found in an anvi&#x27;o contigs database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[splits-taxonomy-txt](artifacts/splits-taxonomy-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-export-state](programs/anvi-export-state)**</span>. <span markdown="1">Export an anvi&#x27;o state into a profile database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[pan-db](artifacts/pan-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[state](artifacts/state) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[state-json](artifacts/state-json) <img src="images/icons/JSON.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-export-structures](programs/anvi-export-structures)**</span>. <span markdown="1">Export .pdb structure files from a structure database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[structure-db](artifacts/structure-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[protein-structure-txt](artifacts/protein-structure-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-gen-contigs-database](programs/anvi-gen-contigs-database)**</span>. <span markdown="1">Generate a new anvi&#x27;o contigs database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-fasta](artifacts/contigs-fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[external-gene-calls](artifacts/external-gene-calls) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-gen-fixation-index-matrix](programs/anvi-gen-fixation-index-matrix)**</span>. <span markdown="1">Generate a pairwise matrix of a fixation indices between samples</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[structure-db](artifacts/structure-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[bin](artifacts/bin) <img src="images/icons/BIN.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[variability-profile-txt](artifacts/variability-profile-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[splits-txt](artifacts/splits-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[fixation-index-matrix](artifacts/fixation-index-matrix) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-gen-gene-consensus-sequences](programs/anvi-gen-gene-consensus-sequences)**</span>. <span markdown="1">Collapse variability for a set of genes across samples</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[genes-fasta](artifacts/genes-fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-gen-gene-level-stats-databases](programs/anvi-gen-gene-level-stats-databases)**</span>. <span markdown="1">A program to compute genes databases for a ginen set of bins stored in an anvi&#x27;o collection. Genes databases store gene-level coverage and detection statistics, and they are usually computed and generated automatically when they are required (such as running anvi-interactive with `--gene-mode` flag). This program allows you to pre-compute them if you don&#x27;t want them to be done all at once</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[collection](artifacts/collection) <img src="images/icons/COLLECTION.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[bin](artifacts/bin) <img src="images/icons/BIN.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[genes-db](artifacts/genes-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-gen-genomes-storage](programs/anvi-gen-genomes-storage)**</span>. <span markdown="1">Create a genome storage from internal and/or external genomes for a pangenome analysis</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[external-genomes](artifacts/external-genomes) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[internal-genomes](artifacts/internal-genomes) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[genomes-storage-db](artifacts/genomes-storage-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-gen-phylogenomic-tree](programs/anvi-gen-phylogenomic-tree)**</span>. <span markdown="1">Generate phylogenomic tree from aligment file</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[concatenated-gene-alignment-fasta](artifacts/concatenated-gene-alignment-fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[phylogeny](artifacts/phylogeny) <img src="images/icons/NEWICK.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-gen-structure-database](programs/anvi-gen-structure-database)**</span>. <span markdown="1">Identifies genes in your contigs database that encode proteins that are homologous to proteins with solved structures. If sufficiently similar homologs are identified, they are used as structural templates to predict the 3D structure of proteins in your contigs database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[pdb-db](artifacts/pdb-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[structure-db](artifacts/structure-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-gen-variability-network](programs/anvi-gen-variability-network)**</span>. <span markdown="1">A program to generate a network description from an anvi&#x27;o variability profile (potentially outdated program)</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[variability-profile](artifacts/variability-profile) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> 
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-gen-variability-profile](programs/anvi-gen-variability-profile)**</span>. <span markdown="1">Generate a table that comprehensively summarizes the variability of nucleotide, codon, or amino acid positions. We call these single nucleotide variants (SNVs), single codon variants (SCVs), and single amino acid variants (SAAVs), respectively</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[structure-db](artifacts/structure-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[bin](artifacts/bin) <img src="images/icons/BIN.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[variability-profile](artifacts/variability-profile) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[splits-txt](artifacts/splits-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[variability-profile-txt](artifacts/variability-profile-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-get-aa-counts](programs/anvi-get-aa-counts)**</span>. <span markdown="1">Fetches the number of times each amino acid occurs from a contigs database in a given bin, set of contigs, or set of genes</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[splits-txt](artifacts/splits-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[collection](artifacts/collection) <img src="images/icons/COLLECTION.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[aa-frequencies-txt](artifacts/aa-frequencies-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-get-codon-frequencies](programs/anvi-get-codon-frequencies)**</span>. <span markdown="1">Get amino acid or codon frequencies of genes in a contigs database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[codon-frequencies-txt](artifacts/codon-frequencies-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[aa-frequencies-txt](artifacts/aa-frequencies-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-get-sequences-for-gene-calls](programs/anvi-get-sequences-for-gene-calls)**</span>. <span markdown="1">A script to get back sequences for gene calls</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[genomes-storage-db](artifacts/genomes-storage-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[genes-fasta](artifacts/genes-fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[external-gene-calls](artifacts/external-gene-calls) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-get-sequences-for-gene-clusters](programs/anvi-get-sequences-for-gene-clusters)**</span>. <span markdown="1">Do cool stuff with gene clusters in anvi&#x27;o pan genomes</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[pan-db](artifacts/pan-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[genomes-storage-db](artifacts/genomes-storage-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[genes-fasta](artifacts/genes-fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[concatenated-gene-alignment-fasta](artifacts/concatenated-gene-alignment-fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[misc-data-items](artifacts/misc-data-items) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-get-sequences-for-hmm-hits](programs/anvi-get-sequences-for-hmm-hits)**</span>. <span markdown="1">Get sequences for HMM hits from many inputs</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[external-genomes](artifacts/external-genomes) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[internal-genomes](artifacts/internal-genomes) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[hmm-source](artifacts/hmm-source) <img src="images/icons/HMM.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[hmm-hits](artifacts/hmm-hits) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[genes-fasta](artifacts/genes-fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[concatenated-gene-alignment-fasta](artifacts/concatenated-gene-alignment-fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-get-short-reads-from-bam](programs/anvi-get-short-reads-from-bam)**</span>. <span markdown="1">Get short reads back from a BAM file with options for compression, splitting of forward and reverse reads, etc</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[bin](artifacts/bin) <img src="images/icons/BIN.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[bam-file](artifacts/bam-file) <img src="images/icons/BAM.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[short-reads-fasta](artifacts/short-reads-fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-get-short-reads-mapping-to-a-gene](programs/anvi-get-short-reads-mapping-to-a-gene)**</span>. <span markdown="1">Recover short reads from BAM files that were mapped to genes you are interested in. It is possible to work with a single gene call, or a bunch of them. Similarly, you can get short reads from a single BAM file, or from many of them</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[bam-file](artifacts/bam-file) <img src="images/icons/BAM.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[short-reads-fasta](artifacts/short-reads-fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-get-split-coverages](programs/anvi-get-split-coverages)**</span>. <span markdown="1">Export splits and the coverage table from database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[collection](artifacts/collection) <img src="images/icons/COLLECTION.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[bin](artifacts/bin) <img src="images/icons/BIN.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[coverages-txt](artifacts/coverages-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-import-collection](programs/anvi-import-collection)**</span>. <span markdown="1">Import an external binning result into anvi&#x27;o</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[pan-db](artifacts/pan-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[collection-txt](artifacts/collection-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[collection](artifacts/collection) <img src="images/icons/COLLECTION.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-import-functions](programs/anvi-import-functions)**</span>. <span markdown="1">Parse and store functional annotation of genes</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[functions-txt](artifacts/functions-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[functions](artifacts/functions) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-import-items-order](programs/anvi-import-items-order)**</span>. <span markdown="1">Import a new items order into an anvi&#x27;o database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[pan-db](artifacts/pan-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[misc-data-items-order-txt](artifacts/misc-data-items-order-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[dendrogram](artifacts/dendrogram) <img src="images/icons/NEWICK.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[phylogeny](artifacts/phylogeny) <img src="images/icons/NEWICK.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[misc-data-items-order](artifacts/misc-data-items-order) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-import-misc-data](programs/anvi-import-misc-data)**</span>. <span markdown="1">Populate additional data or order tables in pan or profile databases for items and layers, OR additional data in contigs databases for nucleotides and amino acids (the Swiss army knife-level serious stuff)</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[pan-db](artifacts/pan-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[misc-data-items-txt](artifacts/misc-data-items-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[dendrogram](artifacts/dendrogram) <img src="images/icons/NEWICK.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[phylogeny](artifacts/phylogeny) <img src="images/icons/NEWICK.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[misc-data-layers-txt](artifacts/misc-data-layers-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[misc-data-layer-orders-txt](artifacts/misc-data-layer-orders-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[misc-data-nucleotides-txt](artifacts/misc-data-nucleotides-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[misc-data-amino-acids-txt](artifacts/misc-data-amino-acids-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[misc-data-items](artifacts/misc-data-items) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[misc-data-layers](artifacts/misc-data-layers) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[misc-data-layer-orders](artifacts/misc-data-layer-orders) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[misc-data-nucleotides](artifacts/misc-data-nucleotides) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[misc-data-amino-acids](artifacts/misc-data-amino-acids) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-import-state](programs/anvi-import-state)**</span>. <span markdown="1">Import an anvi&#x27;o state into a profile database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[pan-db](artifacts/pan-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[state-json](artifacts/state-json) <img src="images/icons/JSON.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[state](artifacts/state) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-import-taxonomy-for-genes](programs/anvi-import-taxonomy-for-genes)**</span>. <span markdown="1">Import gene-level taxonomy into an anvi&#x27;o contigs database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[gene-taxonomy-txt](artifacts/gene-taxonomy-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[gene-taxonomy](artifacts/gene-taxonomy) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-import-taxonomy-for-layers](programs/anvi-import-taxonomy-for-layers)**</span>. <span markdown="1">Import layers-level taxonomy into an anvi&#x27;o additional layer data table in an anvi&#x27;o single-profile database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[single-profile-db](artifacts/single-profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[layer-taxonomy-txt](artifacts/layer-taxonomy-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[layer-taxonomy](artifacts/layer-taxonomy) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-init-bam](programs/anvi-init-bam)**</span>. <span markdown="1">Sort/Index BAM files</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[raw-bam-file](artifacts/raw-bam-file) <img src="images/icons/BAM.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[bam-file](artifacts/bam-file) <img src="images/icons/BAM.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-inspect](programs/anvi-inspect)**</span>. <span markdown="1">Start an anvi&#x27;o inspect interactive interface</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[bin](artifacts/bin) <img src="images/icons/BIN.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[interactive](artifacts/interactive) <img src="images/icons/DISPLAY.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-interactive](programs/anvi-interactive)**</span>. <span markdown="1">Start an anvi&#x27;o server for the interactive interface</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[single-profile-db](artifacts/single-profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[genes-db](artifacts/genes-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[bin](artifacts/bin) <img src="images/icons/BIN.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[view-data](artifacts/view-data) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[dendrogram](artifacts/dendrogram) <img src="images/icons/NEWICK.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[phylogeny](artifacts/phylogeny) <img src="images/icons/NEWICK.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[collection](artifacts/collection) <img src="images/icons/COLLECTION.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[bin](artifacts/bin) <img src="images/icons/BIN.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[interactive](artifacts/interactive) <img src="images/icons/DISPLAY.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[svg](artifacts/svg) <img src="images/icons/SVG.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-matrix-to-newick](programs/anvi-matrix-to-newick)**</span>. <span markdown="1">Takes a distance matrix, returns a newick tree</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[view-data](artifacts/view-data) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[dendrogram](artifacts/dendrogram) <img src="images/icons/NEWICK.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-merge](programs/anvi-merge)**</span>. <span markdown="1">Merge multiple anvio profiles</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[single-profile-db](artifacts/single-profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[misc-data-items-order](artifacts/misc-data-items-order) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-merge-bins](programs/anvi-merge-bins)**</span>. <span markdown="1">Merge a given set of bins in an anvi&#x27;o collection</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[pan-db](artifacts/pan-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[collection](artifacts/collection) <img src="images/icons/COLLECTION.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[bin](artifacts/bin) <img src="images/icons/BIN.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> 
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-meta-pan-genome](programs/anvi-meta-pan-genome)**</span>. <span markdown="1">Convert a pangenome into a metapangenome</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[internal-genomes](artifacts/internal-genomes) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[pan-db](artifacts/pan-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[genomes-storage-db](artifacts/genomes-storage-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[metapangenome](artifacts/metapangenome) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-migrate](programs/anvi-migrate)**</span>. <span markdown="1">Migrate an anvi&#x27;o database or config file to a newer version</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[pan-db](artifacts/pan-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[genes-db](artifacts/genes-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[genomes-storage-db](artifacts/genomes-storage-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[structure-db](artifacts/structure-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[modules-db](artifacts/modules-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[workflow-config](artifacts/workflow-config) <img src="images/icons/JSON.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> 
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-oligotype-linkmers](programs/anvi-oligotype-linkmers)**</span>. <span markdown="1">Takes an anvi&#x27;o linkmers report, generates an oligotyping output</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[linkmers-txt](artifacts/linkmers-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[oligotypes](artifacts/oligotypes) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-pan-genome](programs/anvi-pan-genome)**</span>. <span markdown="1">An anvi&#x27;o program to compute a pangenome from an anvi&#x27;o genome storage</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[genomes-storage-db](artifacts/genomes-storage-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[pan-db](artifacts/pan-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[misc-data-items-order](artifacts/misc-data-items-order) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-profile](programs/anvi-profile)**</span>. <span markdown="1">Creates a single anvi&#x27;o profile database.                    When it is run on a BAM file, depending on the user parameters, the program quantifies                   coverage per nucleotide position (and averages them  per contig), calculates                   single-nucleotide, single-codon, and single-amino acid variants, as well as structural variants                   such as insertion and deletions and stores these data into appropriate tables</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[bam-file](artifacts/bam-file) <img src="images/icons/BAM.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[single-profile-db](artifacts/single-profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[misc-data-items-order](artifacts/misc-data-items-order) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[variability-profile](artifacts/variability-profile) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-refine](programs/anvi-refine)**</span>. <span markdown="1">Start an anvi&#x27;o interactive interactive to manually curate or refine a genome, whether it is a metagenome-assembled, single-cell, or an isolate genome</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[bin](artifacts/bin) <img src="images/icons/BIN.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[bin](artifacts/bin) <img src="images/icons/BIN.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-rename-bins](programs/anvi-rename-bins)**</span>. <span markdown="1">Rename all bins in a given collection (so they have pretty names)</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[collection](artifacts/collection) <img src="images/icons/COLLECTION.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[bin](artifacts/bin) <img src="images/icons/BIN.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[collection](artifacts/collection) <img src="images/icons/COLLECTION.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[bin](artifacts/bin) <img src="images/icons/BIN.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-report-linkmers](programs/anvi-report-linkmers)**</span>. <span markdown="1">Reports sequences stored in one or more BAM files that cover one of more specific nucleotide positions in a reference</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[bam-file](artifacts/bam-file) <img src="images/icons/BAM.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[linkmers-txt](artifacts/linkmers-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-run-hmms](programs/anvi-run-hmms)**</span>. <span markdown="1">This program deals with populating tables that store HMM hits in an anvi&#x27;o contigs database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[hmm-source](artifacts/hmm-source) <img src="images/icons/HMM.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[hmm-hits](artifacts/hmm-hits) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-run-interacdome](programs/anvi-run-interacdome)**</span>. <span markdown="1">Run InteracDome on a contigs database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[interacdome-data](artifacts/interacdome-data) <img src="images/icons/DATA.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[binding-frequencies-txt](artifacts/binding-frequencies-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[misc-data-amino-acids](artifacts/misc-data-amino-acids) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-run-kegg-kofams](programs/anvi-run-kegg-kofams)**</span>. <span markdown="1">Run KOfam HMMs on an anvi&#x27;o contigs database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[kegg-data](artifacts/kegg-data) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[kegg-functions](artifacts/kegg-functions) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[functions](artifacts/functions) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-run-ncbi-cogs](programs/anvi-run-ncbi-cogs)**</span>. <span markdown="1">This program runs NCBI&#x27;s COGs to associate genes in an anvi&#x27;o contigs database with functions. COGs database was been designed as an attempt to classify proteins from completely sequenced genomes on the basis of the orthology concept.</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[cogs-data](artifacts/cogs-data) <img src="images/icons/DATA.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[functions](artifacts/functions) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-run-pfams](programs/anvi-run-pfams)**</span>. <span markdown="1">Run Pfam on Contigs Database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[pfams-data](artifacts/pfams-data) <img src="images/icons/DATA.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[functions](artifacts/functions) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-run-scg-taxonomy](programs/anvi-run-scg-taxonomy)**</span>. <span markdown="1">The purpose of this program is to affiliate single-copy core genes in an anvi&#x27;o contigs database with taxonomic names. A properly setup local SCG taxonomy database is required for this program to perform properly. After its successful run, `anvi-estimate-scg-taxonomy` will be useful to estimate taxonomy at genome-, collection-, or metagenome-level)</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[scgs-taxonomy-db](artifacts/scgs-taxonomy-db) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[scgs-taxonomy](artifacts/scgs-taxonomy) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-run-trna-taxonomy](programs/anvi-run-trna-taxonomy)**</span>. <span markdown="1">The purpose of this program is to affiliate tRNA gene sequences in an anvi&#x27;o contigs database with taxonomic names. A properly setup local tRNA taxonomy database is required for this program to perform properly. After its successful run, `anvi-estimate-trna-taxonomy` will be useful to estimate taxonomy at genome-, collection-, or metagenome-level).</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[trna-taxonomy-db](artifacts/trna-taxonomy-db) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[trna-taxonomy](artifacts/trna-taxonomy) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-run-workflow](programs/anvi-run-workflow)**</span>. <span markdown="1">Execute, manage, parallelize, and troubleshoot entire &#x27;omics workflows and chain together anvi&#x27;o and third party programs</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[samples-txt](artifacts/samples-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[fasta-txt](artifacts/fasta-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[workflow-config](artifacts/workflow-config) <img src="images/icons/JSON.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[contigs-workflow](artifacts/contigs-workflow) <img src="images/icons/WORKFLOW.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[metagenomics-workflow](artifacts/metagenomics-workflow) <img src="images/icons/WORKFLOW.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[pangenomics-workflow](artifacts/pangenomics-workflow) <img src="images/icons/WORKFLOW.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[phylogenomics-workflow](artifacts/phylogenomics-workflow) <img src="images/icons/WORKFLOW.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[trnaseq-workflow](artifacts/trnaseq-workflow) <img src="images/icons/WORKFLOW.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-scan-trnas](programs/anvi-scan-trnas)**</span>. <span markdown="1">Identify and store tRNA genes in a contigs database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[hmm-hits](artifacts/hmm-hits) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-search-functions](programs/anvi-search-functions)**</span>. <span markdown="1">Search functions in an anvi&#x27;o contigs database or genomes storage. Basically, this program searches for one or more search terms you define in functional annotations of genes in an anvi&#x27;o contigs database, and generates multiple reports. The default report simply tells you which contigs contain genes with functions matching to serach terms you used, useful for viewing in the interface. You can also request a much more comprehensive report, which gives you anything you might need to know for each hit and serach term</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[genomes-storage-db](artifacts/genomes-storage-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[functions-txt](artifacts/functions-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-search-sequence-motifs](programs/anvi-search-sequence-motifs)**</span>. <span markdown="1">A program to find one or more sequence motifs in contig or gene sequences, and store their frequencies</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[genes-db](artifacts/genes-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[misc-data-items](artifacts/misc-data-items) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[misc-data-layers](artifacts/misc-data-layers) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-setup-interacdome](programs/anvi-setup-interacdome)**</span>. <span markdown="1">Setup InteracDome data</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> 
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[interacdome-data](artifacts/interacdome-data) <img src="images/icons/DATA.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-setup-kegg-kofams](programs/anvi-setup-kegg-kofams)**</span>. <span markdown="1">Download and setup KEGG KOfam HMM profiles and KEGG MODULE data</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> 
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[kegg-data](artifacts/kegg-data) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[modules-db](artifacts/modules-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-setup-ncbi-cogs](programs/anvi-setup-ncbi-cogs)**</span>. <span markdown="1">Download and setup NCBI&#x27;s Clusters of Orthologous Groups database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> 
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[cogs-data](artifacts/cogs-data) <img src="images/icons/DATA.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-setup-pdb-database](programs/anvi-setup-pdb-database)**</span>. <span markdown="1">Setup or update an offline database of representative PDB structures clustered at 95%</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> 
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[pdb-db](artifacts/pdb-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-setup-pfams](programs/anvi-setup-pfams)**</span>. <span markdown="1">Download and setup Pfam data from the EBI</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> 
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[pfams-data](artifacts/pfams-data) <img src="images/icons/DATA.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-setup-scg-taxonomy](programs/anvi-setup-scg-taxonomy)**</span>. <span markdown="1">The purpose of this program is to download necessary information from GTDB (https://gtdb.ecogenomic.org/), and set it up in such a way that your anvi&#x27;o installation is able to assign taxonomy to single-copy core genes using `anvi-run-scg-taxonomy` and estimate taxonomy for genomes or metagenomes using `anvi-estimate-scg-taxonomy`)</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> 
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[scgs-taxonomy-db](artifacts/scgs-taxonomy-db) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-setup-trna-taxonomy](programs/anvi-setup-trna-taxonomy)**</span>. <span markdown="1">The purpose of this program is to setup necessary databases for tRNA genes collected from GTDB (https://gtdb.ecogenomic.org/), genomes in your local anvi&#x27;o installation so taxonomy information for a given set of tRNA sequences can be identified using `anvi-run-trna-taxonomy` and made sense of via `anvi-estimate-trna-taxonomy`)</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> 
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[trna-taxonomy-db](artifacts/trna-taxonomy-db) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-show-collections-and-bins](programs/anvi-show-collections-and-bins)**</span>. <span markdown="1">A script to display collections stored in an anvi&#x27;o profile or pan database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[pan-db](artifacts/pan-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> 
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-show-misc-data](programs/anvi-show-misc-data)**</span>. <span markdown="1">Show all misc data keys in all misc data tables</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[pan-db](artifacts/pan-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> 
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-split](programs/anvi-split)**</span>. <span markdown="1">Split an anvi&#x27;o pan or profile database into smaller, self-contained pieces. Provide either a genomes-storage and pan database or a profile and contigs database pair, and you&#x27;ll get back directories of individual projects for each bin  that can be treated as smaller anvi&#x27;o projects</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[genomes-storage-db](artifacts/genomes-storage-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[pan-db](artifacts/pan-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[collection](artifacts/collection) <img src="images/icons/COLLECTION.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[split-bins](artifacts/split-bins) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-summarize](programs/anvi-summarize)**</span>. <span markdown="1">Summarizer for anvi&#x27;o pan or profile db&#x27;s. Essentially, this program takes a collection id along with either a profile database and a contigs database or a pan database and a genomes storage and generates a static HTML output for what is described in a given collection. The output directory will contain almost everything any downstream analysis may need, and can be displayed using a browser without the need for an anvi&#x27;o installation. For this reason alone, reporting summary outputs as supplementary data with publications is a great idea for transparency and reproducibility</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[collection](artifacts/collection) <img src="images/icons/COLLECTION.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[pan-db](artifacts/pan-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[genomes-storage-db](artifacts/genomes-storage-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[summary](artifacts/summary) <img src="images/icons/SUMMARY.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-trnaseq](programs/anvi-trnaseq)**</span>. <span markdown="1">A program to process raw tRNA-seq dataset, which is the sequencing of tRNA transcripts in a given sample, to generate an anvi&#x27;o tRNA-seq database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[trnaseq-fasta](artifacts/trnaseq-fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[trnaseq-db](artifacts/trnaseq-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-update-db-description](programs/anvi-update-db-description)**</span>. <span markdown="1">Update the description in an anvi&#x27;o database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[pan-db](artifacts/pan-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[genomes-storage-db](artifacts/genomes-storage-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> 
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-update-structure-database](programs/anvi-update-structure-database)**</span>. <span markdown="1">Add or re-run genes from an already existing structure database. All settings used to generate your database will be used in this program</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[structure-db](artifacts/structure-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> 
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-script-add-default-collection](programs/anvi-script-add-default-collection)**</span>. <span markdown="1">A script to add a &#x27;DEFAULT&#x27; collection in an anvi&#x27;o pan or profile database with a bin named &#x27;EVERYTHING&#x27; that describes all items available in the profile database</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[pan-db](artifacts/pan-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[collection](artifacts/collection) <img src="images/icons/COLLECTION.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[bin](artifacts/bin) <img src="images/icons/BIN.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-script-augustus-output-to-external-gene-calls](programs/anvi-script-augustus-output-to-external-gene-calls)**</span>. <span markdown="1">Takes in gene calls by AUGUSTUS v3.3.3, generates an anvi&#x27;o external gene calls file. It may work well with other versions of AUGUSTUS, too. It is just no one has tested the script with different versions of the program</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[augustus-gene-calls](artifacts/augustus-gene-calls) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[external-gene-calls](artifacts/external-gene-calls) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-script-calculate-pn-ps-ratio](programs/anvi-script-calculate-pn-ps-ratio)**</span>. <span markdown="1">FIXME</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[variability-profile-txt](artifacts/variability-profile-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[pn-ps-data](artifacts/pn-ps-data) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-script-checkm-tree-to-interactive](programs/anvi-script-checkm-tree-to-interactive)**</span>. <span markdown="1">A helper script to convert CheckM trees into anvio interactive with taxonomy information</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[phylogeny](artifacts/phylogeny) <img src="images/icons/NEWICK.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[interactive](artifacts/interactive) <img src="images/icons/DISPLAY.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-script-compute-ani-for-fasta](programs/anvi-script-compute-ani-for-fasta)**</span>. <span markdown="1">Run ANI between contigs in a single FASTA file</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[fasta](artifacts/fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[genome-similarity](artifacts/genome-similarity) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-script-filter-fasta-by-blast](programs/anvi-script-filter-fasta-by-blast)**</span>. <span markdown="1">Filter FASTA file according to BLAST table (remove sequences with bad BLAST alignment)</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-fasta](artifacts/contigs-fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[blast-table](artifacts/blast-table) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[contigs-fasta](artifacts/contigs-fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-script-filter-hmm-hits-table](programs/anvi-script-filter-hmm-hits-table)**</span>. <span markdown="1">Add or remove entries in a contigDB hmm_hits table.</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[hmm-source](artifacts/hmm-source) <img src="images/icons/HMM.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[hmm-hits](artifacts/hmm-hits) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[hmm-hits](artifacts/hmm-hits) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-script-fix-homopolymer-indels](programs/anvi-script-fix-homopolymer-indels)**</span>. <span markdown="1">Corrects homopolymer-region associated INDELs in a given genome based on a reference genome. The most effective use of this script is when the input genome is a genome reconstructed by minION long reads, and the reference genome is one that is of high-quality. Essentially, this script will BLAST the genome you wish to correct against the reference genome you provide, identify INDELs in the BLAST results that are exclusively associated with homopolymer regions, and will take the reference genome as a guide to correct the input sequences, and report a new FASTA file. You can use the output FASTA file that is fixed as the input FASTA file over and over again to see if you can eliminate all homopolymer-associated INDELs</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[fasta](artifacts/fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[fasta](artifacts/fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-script-gen-distribution-of-genes-in-a-bin](programs/anvi-script-gen-distribution-of-genes-in-a-bin)**</span>. <span markdown="1">Quantify the detection of genes in genomes in metagenomes to identify the environmental core. This is a helper script for anvi&#x27;o metapangenomic workflow</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[collection](artifacts/collection) <img src="images/icons/COLLECTION.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[bin](artifacts/bin) <img src="images/icons/BIN.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[view-data](artifacts/view-data) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[misc-data-items-txt](artifacts/misc-data-items-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-script-gen-functions-per-group-stats-output](programs/anvi-script-gen-functions-per-group-stats-output)**</span>. <span markdown="1">Generate a TAB delimited file for the distribution of functions across groups of genomes/metagenomes</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[functions](artifacts/functions) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[genomes-storage-db](artifacts/genomes-storage-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[internal-genomes](artifacts/internal-genomes) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[external-genomes](artifacts/external-genomes) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[interactive](artifacts/interactive) <img src="images/icons/DISPLAY.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-script-gen-genomes-file](programs/anvi-script-gen-genomes-file)**</span>. <span markdown="1">Generate an external genomes or internal genomes file</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[profile-db](artifacts/profile-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[collection](artifacts/collection) <img src="images/icons/COLLECTION.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[external-genomes](artifacts/external-genomes) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[internal-genomes](artifacts/internal-genomes) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-script-gen-hmm-hits-matrix-across-genomes](programs/anvi-script-gen-hmm-hits-matrix-across-genomes)**</span>. <span markdown="1">A simple script to generate a TAB-delimited file that reports the frequency of HMM hits for a given HMM source across contigs databases</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[external-genomes](artifacts/external-genomes) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[internal-genomes](artifacts/internal-genomes) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[hmm-source](artifacts/hmm-source) <img src="images/icons/HMM.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[hmm-hits](artifacts/hmm-hits) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[hmm-hits-matrix-txt](artifacts/hmm-hits-matrix-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-script-gen-pseudo-paired-reads-from-fastq](programs/anvi-script-gen-pseudo-paired-reads-from-fastq)**</span>. <span markdown="1">A script that takes a FASTQ file that is not paired-end (i.e., R1 alone) and converts it into two FASTQ files that are paired-end (i.e., R1 and R2). This is a quick-and-dirty workaround that halves each read from the original FASTQ and puts one half in the FASTQ file for R1 and puts the reverse-complement of the second half in the FASTQ file for R2. If you&#x27;ve ended up here, things have clearly not gone very well for you, and Evan, who battled similar battles and ended up implementing this solution wholeheartedly sympathizes</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[short-reads-fasta](artifacts/short-reads-fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[short-reads-fasta](artifacts/short-reads-fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-script-gen-short-reads](programs/anvi-script-gen-short-reads)**</span>. <span markdown="1">Generate short reads from contigs. Useful to reconstruct mock data sets from already assembled contigs</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[configuration-ini](artifacts/configuration-ini) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[short-reads-fasta](artifacts/short-reads-fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-script-gen_stats_for_single_copy_genes.py](programs/anvi-script-gen_stats_for_single_copy_genes.py)**</span>. <span markdown="1">A simple script to generate info from search tables, given a contigs-db</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> 
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[genes-stats](artifacts/genes-stats) <img src="images/icons/STATS.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-script-get-coverage-from-bam](programs/anvi-script-get-coverage-from-bam)**</span>. <span markdown="1">Get nucleotide-level, contig-level, or bin-level coverage values from a BAM file</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[bam-file](artifacts/bam-file) <img src="images/icons/BAM.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[collection-txt](artifacts/collection-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[coverages-txt](artifacts/coverages-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-script-get-hmm-hits-per-gene-call](programs/anvi-script-get-hmm-hits-per-gene-call)**</span>. <span markdown="1">A simple script to generate a TAB-delimited file gene caller IDs and their HMM hits for a given HMM source</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[hmm-source](artifacts/hmm-source) <img src="images/icons/HMM.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[hmm-hits](artifacts/hmm-hits) <img src="images/icons/CONCEPT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[functions-txt](artifacts/functions-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-script-get-primer-matches](programs/anvi-script-get-primer-matches)**</span>. <span markdown="1">You provide this program with FASTQ files for one or more samples AND one or more short sequences, and it collects reads from FASTQ files that matches to your sequences. This tool can be most powerful if you want to collect all short reads from one or more metagenomes that are downstream to a known sequence. Using the comprehensive output files you can analyze the diversity of seuqences visually, manually, or using established strategies such as oligotyping.</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[samples-txt](artifacts/samples-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[short-reads-fasta](artifacts/short-reads-fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-script-merge-collections](programs/anvi-script-merge-collections)**</span>. <span markdown="1">Generate an additional data file from multiple collections</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[contigs-db](artifacts/contigs-db) <img src="images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[collection-txt](artifacts/collection-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> 
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-script-pfam-accessions-to-hmms-directory](programs/anvi-script-pfam-accessions-to-hmms-directory)**</span>. <span markdown="1">You give this program one or more PFAM accession ids, and it generates an anvi&#x27;o compatible HMM directory to be used with `anvi-run-hmms`</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> 
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[hmm-source](artifacts/hmm-source) <img src="images/icons/HMM.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-script-process-genbank](programs/anvi-script-process-genbank)**</span>. <span markdown="1">This script takes a GenBank file, and outputs a FASTA file, as well as two additional TAB-delimited output files for external gene calls and gene functions that can be used with the programs `anvi-gen-contigs-database` and `anvi-import-functions`</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[genbank-file](artifacts/genbank-file) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[contigs-fasta](artifacts/contigs-fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[external-gene-calls](artifacts/external-gene-calls) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[functions-txt](artifacts/functions-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-script-process-genbank-metadata](programs/anvi-script-process-genbank-metadata)**</span>. <span markdown="1">This script takes the &#x27;metadata&#x27; output of the program `ncbi-genome-download` (see [https://github.com/kblin/ncbi-genome-download](https://github.com/kblin/ncbi-genome-download) for details), and processes each GenBank file found in the metadata file to generate a FASTA file, as well as genes and functions files for each entry. Plus, it autmatically generates a FASTA TXT file descriptor for anvi&#x27;o snakemake workflows. So it is a multi-talented program like that</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> 
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[contigs-fasta](artifacts/contigs-fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[functions-txt](artifacts/functions-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[external-gene-calls](artifacts/external-gene-calls) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-script-reformat-fasta](programs/anvi-script-reformat-fasta)**</span>. <span markdown="1">Reformat FASTA file (remove contigs based on length, or based on a given list of deflines, and/or generate an output with simpler names)</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[fasta](artifacts/fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[contigs-fasta](artifacts/contigs-fasta) <img src="images/icons/FASTA.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-script-snvs-to-interactive](programs/anvi-script-snvs-to-interactive)**</span>. <span markdown="1">Take the output of anvi-gen-variability-profile, prepare an output for interactive interface</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[variability-profile-txt](artifacts/variability-profile-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[interactive](artifacts/interactive) <img src="images/icons/DISPLAY.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-script-transpose-matrix](programs/anvi-script-transpose-matrix)**</span>. <span markdown="1">Transpose a TAB-delimited file</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[view-data](artifacts/view-data) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[functions-txt](artifacts/functions-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[misc-data-items-txt](artifacts/misc-data-items-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[misc-data-layers-txt](artifacts/misc-data-layers-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[gene-calls-txt](artifacts/gene-calls-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r" markdown="1">[linkmers-txt](artifacts/linkmers-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[view-data](artifacts/view-data) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[functions-txt](artifacts/functions-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[misc-data-items-txt](artifacts/misc-data-items-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[misc-data-layers-txt](artifacts/misc-data-layers-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[gene-calls-txt](artifacts/gene-calls-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p" markdown="1">[linkmers-txt](artifacts/linkmers-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>

<div style="width:100%;">
<table class="programs-table">
<tbody>
<tr style="border:none;">
    <td class="program-td">
        <span class="artifact-emoji">🔥</span> <span markdown="1">**[anvi-script-variability-to-vcf](programs/anvi-script-variability-to-vcf)**</span>. <span markdown="1">A script to convert SNV output obtained from anvi-gen-variability-profile to the standard VCF format</span>.
    </td>
</tr>
<tr>
    <td class="artifact-r-td">
    <span class="artifact-emoji">🧀</span> <span class="artifact-r" markdown="1">[variability-profile-txt](artifacts/variability-profile-txt) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
<tr>
    <td class="artifact-p-td">
    <span class="artifact-emoji">🍕</span> <span class="artifact-p" markdown="1">[vcf](artifacts/vcf) <img src="images/icons/TXT.png" class="artifact-icon-mini" /> </span>
    </td>
</tr>
</tbody>
</table>
</div>


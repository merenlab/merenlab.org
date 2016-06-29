---
layout: page
title: Software
modified: 2015-02-05
excerpt: "Software we develop"
redirect_from: /projects/
comments: true
---

{% include _toc.html %}

## Anvi'o

<div class="quotable">
<a href="{{ site.url }}/software/anvio/">Anvi'o</a> is an advanced analysis and visualization platform for ‘omics data. Its interactive interface facilitates the management of metagenomic contigs and associated data for automatic or human-guided identification of genome bins, and their curation. The extensible visualization approach distills multiple dimensions of information for each contig into a single, intuitive display, offering a dynamic and unified work environment for data exploration, manipulation and reporting. Beyond its easy-to-use interface, the advanced modular architecture of anvi’o as a platform allows users with programming skills to implement and test novel ideas with minimal effort. Please see the <a href="{{ site.url }}/software/anvio/">anvi'o project page</a> for details, and the <a href="http://github.com/meren/anvio">codebase</a> to follow the development and/or participate.
</div>

{% include _project-anvio-images.html %}

## Oligotyping

<div class="quotable">
<a href="{{ site.url }}/software/oligotyping/">Oligotyping</a> is a human-guided computational approach that makes it possible to decompose very closely related taxa at one nucleotide resolution. It is generally applied to the high-throughput sequencing of bacterial marker gene amplicons amplified from environmental samples (such as 16S rRNA gene). See the <a href="{{ site.url  }}/software/oligotyping/">project page</a> for details.
</div>

{% include _project-oligotyping-images.html %}

## Minimum Entropy Decomposition

<div class="quotable">
<a href="{{ site.url }}/software/med/">MED</a> is an information theory-based clustering algorithm for sensitive partitioning of high-throughput marker gene sequences. The source code is distributed through the oligotyping pipeline. See the <a href="{{ site.url  }}/software/med/">project page</a> for more.
</div>

{% include _project-med-images.html %}

## Illumina Utilities Library

<div class="quotable">
A lightweight and high-performance library to analyze raw Illumina data. It contains programs for demultiplexing, quality filtering, and mergeing partially or fully overlapping reads. Illumina utils has been a core component of the sequencing operations at the MBL. The source code, installation instructions and examples are available through its GitHub repository:
<br />
<br />
<a href="https://github.com/meren/illumina-utils">https://github.com/meren/illumina-utils</a>
</div>

## BLAST Filtering Pipeline

<div class="quotable">
A metagenomic short read filtering software that uses a flexible configuration format. It allows users to define a chain of genomic filters, each of which perform on the output data provided by the previous filter, to filter reads out from sequencing data. It can exploit Sun Grid Engine and distribute individual processes. The source code is available through its GitHub repository:
<br />
<br />
<a href="https://github.com/meren/BLAST-filtering-pipeline">https://github.com/meren/BLAST-filtering-pipeline</a>
</div>


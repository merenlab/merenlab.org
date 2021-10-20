---
layout: page
title: contigs-db [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/contigs-db
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/DB.png" alt="DB" style="width:100px; border:none" />

A DB-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-gen-contigs-database](../../programs/anvi-gen-contigs-database)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-cluster-contigs](../../programs/anvi-cluster-contigs)</span> <span class="artifact-r">[anvi-compute-completeness](../../programs/anvi-compute-completeness)</span> <span class="artifact-r">[anvi-db-info](../../programs/anvi-db-info)</span> <span class="artifact-r">[anvi-delete-functions](../../programs/anvi-delete-functions)</span> <span class="artifact-r">[anvi-delete-hmms](../../programs/anvi-delete-hmms)</span> <span class="artifact-r">[anvi-display-contigs-stats](../../programs/anvi-display-contigs-stats)</span> <span class="artifact-r">[anvi-display-metabolism](../../programs/anvi-display-metabolism)</span> <span class="artifact-r">[anvi-display-structure](../../programs/anvi-display-structure)</span> <span class="artifact-r">[anvi-estimate-genome-completeness](../../programs/anvi-estimate-genome-completeness)</span> <span class="artifact-r">[anvi-estimate-metabolism](../../programs/anvi-estimate-metabolism)</span> <span class="artifact-r">[anvi-estimate-scg-taxonomy](../../programs/anvi-estimate-scg-taxonomy)</span> <span class="artifact-r">[anvi-estimate-trna-taxonomy](../../programs/anvi-estimate-trna-taxonomy)</span> <span class="artifact-r">[anvi-export-contigs](../../programs/anvi-export-contigs)</span> <span class="artifact-r">[anvi-export-functions](../../programs/anvi-export-functions)</span> <span class="artifact-r">[anvi-export-gene-calls](../../programs/anvi-export-gene-calls)</span> <span class="artifact-r">[anvi-export-gene-coverage-and-detection](../../programs/anvi-export-gene-coverage-and-detection)</span> <span class="artifact-r">[anvi-export-locus](../../programs/anvi-export-locus)</span> <span class="artifact-r">[anvi-export-misc-data](../../programs/anvi-export-misc-data)</span> <span class="artifact-r">[anvi-export-splits-and-coverages](../../programs/anvi-export-splits-and-coverages)</span> <span class="artifact-r">[anvi-export-splits-taxonomy](../../programs/anvi-export-splits-taxonomy)</span> <span class="artifact-r">[anvi-gen-fixation-index-matrix](../../programs/anvi-gen-fixation-index-matrix)</span> <span class="artifact-r">[anvi-gen-gene-consensus-sequences](../../programs/anvi-gen-gene-consensus-sequences)</span> <span class="artifact-r">[anvi-gen-gene-level-stats-databases](../../programs/anvi-gen-gene-level-stats-databases)</span> <span class="artifact-r">[anvi-gen-structure-database](../../programs/anvi-gen-structure-database)</span> <span class="artifact-r">[anvi-gen-variability-profile](../../programs/anvi-gen-variability-profile)</span> <span class="artifact-r">[anvi-get-aa-counts](../../programs/anvi-get-aa-counts)</span> <span class="artifact-r">[anvi-get-codon-frequencies](../../programs/anvi-get-codon-frequencies)</span> <span class="artifact-r">[anvi-get-pn-ps-ratio](../../programs/anvi-get-pn-ps-ratio)</span> <span class="artifact-r">[anvi-get-sequences-for-gene-calls](../../programs/anvi-get-sequences-for-gene-calls)</span> <span class="artifact-r">[anvi-get-sequences-for-hmm-hits](../../programs/anvi-get-sequences-for-hmm-hits)</span> <span class="artifact-r">[anvi-get-short-reads-from-bam](../../programs/anvi-get-short-reads-from-bam)</span> <span class="artifact-r">[anvi-get-short-reads-mapping-to-a-gene](../../programs/anvi-get-short-reads-mapping-to-a-gene)</span> <span class="artifact-r">[anvi-get-split-coverages](../../programs/anvi-get-split-coverages)</span> <span class="artifact-r">[anvi-import-collection](../../programs/anvi-import-collection)</span> <span class="artifact-r">[anvi-import-functions](../../programs/anvi-import-functions)</span> <span class="artifact-r">[anvi-import-misc-data](../../programs/anvi-import-misc-data)</span> <span class="artifact-r">[anvi-import-taxonomy-for-genes](../../programs/anvi-import-taxonomy-for-genes)</span> <span class="artifact-r">[anvi-inspect](../../programs/anvi-inspect)</span> <span class="artifact-r">[anvi-interactive](../../programs/anvi-interactive)</span> <span class="artifact-r">[anvi-merge](../../programs/anvi-merge)</span> <span class="artifact-r">[anvi-migrate](../../programs/anvi-migrate)</span> <span class="artifact-r">[anvi-profile](../../programs/anvi-profile)</span> <span class="artifact-r">[anvi-profile-blitz](../../programs/anvi-profile-blitz)</span> <span class="artifact-r">[anvi-refine](../../programs/anvi-refine)</span> <span class="artifact-r">[anvi-rename-bins](../../programs/anvi-rename-bins)</span> <span class="artifact-r">[anvi-report-inversions](../../programs/anvi-report-inversions)</span> <span class="artifact-r">[anvi-run-hmms](../../programs/anvi-run-hmms)</span> <span class="artifact-r">[anvi-run-interacdome](../../programs/anvi-run-interacdome)</span> <span class="artifact-r">[anvi-run-kegg-kofams](../../programs/anvi-run-kegg-kofams)</span> <span class="artifact-r">[anvi-run-ncbi-cogs](../../programs/anvi-run-ncbi-cogs)</span> <span class="artifact-r">[anvi-run-pfams](../../programs/anvi-run-pfams)</span> <span class="artifact-r">[anvi-run-scg-taxonomy](../../programs/anvi-run-scg-taxonomy)</span> <span class="artifact-r">[anvi-run-trna-taxonomy](../../programs/anvi-run-trna-taxonomy)</span> <span class="artifact-r">[anvi-scan-trnas](../../programs/anvi-scan-trnas)</span> <span class="artifact-r">[anvi-search-functions](../../programs/anvi-search-functions)</span> <span class="artifact-r">[anvi-search-palindromes](../../programs/anvi-search-palindromes)</span> <span class="artifact-r">[anvi-search-sequence-motifs](../../programs/anvi-search-sequence-motifs)</span> <span class="artifact-r">[anvi-show-misc-data](../../programs/anvi-show-misc-data)</span> <span class="artifact-r">[anvi-split](../../programs/anvi-split)</span> <span class="artifact-r">[anvi-summarize](../../programs/anvi-summarize)</span> <span class="artifact-r">[anvi-summarize-blitz](../../programs/anvi-summarize-blitz)</span> <span class="artifact-r">[anvi-update-db-description](../../programs/anvi-update-db-description)</span> <span class="artifact-r">[anvi-update-structure-database](../../programs/anvi-update-structure-database)</span> <span class="artifact-r">[anvi-script-add-default-collection](../../programs/anvi-script-add-default-collection)</span> <span class="artifact-r">[anvi-script-filter-hmm-hits-table](../../programs/anvi-script-filter-hmm-hits-table)</span> <span class="artifact-r">[anvi-script-gen-distribution-of-genes-in-a-bin](../../programs/anvi-script-gen-distribution-of-genes-in-a-bin)</span> <span class="artifact-r">[anvi-script-gen-genomes-file](../../programs/anvi-script-gen-genomes-file)</span> <span class="artifact-r">[anvi-script-gen_stats_for_single_copy_genes.py](../../programs/anvi-script-gen_stats_for_single_copy_genes.py)</span> <span class="artifact-r">[anvi-script-get-hmm-hits-per-gene-call](../../programs/anvi-script-get-hmm-hits-per-gene-call)</span> <span class="artifact-r">[anvi-script-merge-collections](../../programs/anvi-script-merge-collections)</span> <span class="artifact-r">[anvi-script-permute-trnaseq-seeds](../../programs/anvi-script-permute-trnaseq-seeds)</span></p>


## Description

A contigs database is an anvi'o database that **contains key information associated with your sequences**.

In a way, **an anvi'o contigs database is a modern, more talented form of a FASTA file**, where you can store additional information about your sequences in it and others can query and use it. Information storage and access is primarily done by anvi'o programs, however, it can also be done through the command line interface or programmatically.

The information a contigs database contains about its sequences can include the positions of open reading frames, tetra-nucleotide frequencies, functional and taxonomic annotations, information on individual nucleotide or amino acid positions, and more.

### Another (less computation-heavy) way of thinking about it

When working in anvi'o, you'll need to be able to access previous analysis done on a genome or transcriptome. To do this, anvi'o uses tools like contigs databases instead of regular fasta files. So, you'll want to convert the data that you have into a contigs database to use other anvi'o programs (using <span class="artifact-n">[anvi-gen-contigs-database](/software/anvio/help/7.1/programs/anvi-gen-contigs-database)</span>). As seen on the page for <span class="artifact-n">[metagenomes](/software/anvio/help/7.1/artifacts/metagenomes)</span>, you can then use this contigs database instead of your fasta file for all of your anvi'o needs.

In short, to get the most out of your data in anvi'o, you'll want to use your data (which was probably originally in a <span class="artifact-n">[fasta](/software/anvio/help/7.1/artifacts/fasta)</span> file) to create both a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> and a <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span>. That way, anvi'o is able to keep track of many different kinds of analysis and you can easily interact with other anvi'o programs.

## Usage Information

### Creating and populating a contigs database

Contigs databases will be initialized using **<span class="artifact-n">[anvi-gen-contigs-database](/software/anvio/help/7.1/programs/anvi-gen-contigs-database)</span>** using a <span class="artifact-n">[contigs-fasta](/software/anvio/help/7.1/artifacts/contigs-fasta)</span>. This will compute the k-mer frequencies for each contig, soft-split your contigs, and identify open reading frames. To populate a contigs database with more information, you can then run various other programs.

**Key programs that populate an anvi'o contigs database with essential information** include,

* <span class="artifact-n">[anvi-run-hmms](/software/anvio/help/7.1/programs/anvi-run-hmms)</span> (which uses HMMs to annotate your genes against an <span class="artifact-n">[hmm-source](/software/anvio/help/7.1/artifacts/hmm-source)</span>)
* <span class="artifact-n">[anvi-run-scg-taxonomy](/software/anvio/help/7.1/programs/anvi-run-scg-taxonomy)</span> (which associates its single-copy core gene with taxonomic data)
* <span class="artifact-n">[anvi-scan-trnas](/software/anvio/help/7.1/programs/anvi-scan-trnas)</span> (which identifies the tRNA genes)
* <span class="artifact-n">[anvi-run-ncbi-cogs](/software/anvio/help/7.1/programs/anvi-run-ncbi-cogs)</span> (which tries to assign functions to your genes using the COGs database)

Once an anvi'o contigs database is generated and populated with information, it is **always a good idea to run <span class="artifact-n">[anvi-display-contigs-stats](/software/anvio/help/7.1/programs/anvi-display-contigs-stats)</span>** to see a numerical summary of its contents.

Other programs you can run to populate a contigs database with functions include,

* <span class="artifact-n">[anvi-run-kegg-kofams](/software/anvio/help/7.1/programs/anvi-run-kegg-kofams)</span> (which annotates the genes in the database with the KEGG KOfam database)

### Analysis on a populated contigs database

Other essential programs that read from a contigs database and yield key information include <span class="artifact-n">[anvi-estimate-genome-completeness](/software/anvio/help/7.1/programs/anvi-estimate-genome-completeness)</span>, <span class="artifact-n">[anvi-get-sequences-for-hmm-hits](/software/anvio/help/7.1/programs/anvi-get-sequences-for-hmm-hits)</span>, and <span class="artifact-n">[anvi-estimate-scg-taxonomy](/software/anvio/help/7.1/programs/anvi-estimate-scg-taxonomy)</span>.

If you wish to run programs like <span class="artifact-n">[anvi-cluster-contigs](/software/anvio/help/7.1/programs/anvi-cluster-contigs)</span>, <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/7.1/programs/anvi-estimate-metabolism)</span>, and <span class="artifact-n">[anvi-gen-gene-level-stats-databases](/software/anvio/help/7.1/programs/anvi-gen-gene-level-stats-databases)</span>, or view your database with <span class="artifact-n">[anvi-interactive](/software/anvio/help/7.1/programs/anvi-interactive)</span>, you'll need to first use your contigs database to create a <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span>.

## Variants

Contigs databases, like <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span>s, are allowed have different variants, though the only currently implemented variant, the <span class="artifact-n">[trnaseq-contigs-db](/software/anvio/help/7.1/artifacts/trnaseq-contigs-db)</span>, is for tRNA transcripts from tRNA-seq experiments. The default variant stored for "standard" contigs databases is `unknown`. Variants should indicate that substantially different information is stored in the database. For instance, open reading frames are applicable to protein-coding genes but not tRNA transcripts, so ORF data is not recorded for the `trnaseq` variant. The $(trnaseq-workflow)s generates <span class="artifact-n">[trnaseq-contigs-db](/software/anvio/help/7.1/artifacts/trnaseq-contigs-db)</span>s using a very different approach to <span class="artifact-n">[anvi-gen-contigs-database](/software/anvio/help/7.1/programs/anvi-gen-contigs-database)</span>.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/contigs-db.md) to update this information.


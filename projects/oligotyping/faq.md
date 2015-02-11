---
layout: project-oligotyping
title: "Frequently asked Questions"
tags: [faq]
modified: 2015-02-02T20:53:07.573882-04:00
comments: true
image:
  feature: header.jpg
---

{% include _toc.html %}

### What is oligotyping?

Oligotyping is a supervised method to identify closely related but distinct bacterial taxa in high-throughput sequencing data of 16S Ribosomal RNA gene.

The input of oligotyping is the final units of OTU clustering or taxonomy analysis (such as an OTU, or a genus).

If there is an OTU or genus that appears to be present in all datasets, oligotyping can help investigating this critical question: “Is there any unexplained diversity in this unit?“.

The distribution of oligotypes recovered from a given OTU or genus among samples in the datasets may help the better understand the underlying ecology.

### What is an oligotype?

An oligotype is a group of reads that are binned together based on the nucleotides they possess at the location of high variation explained by the Shannon entropy analysis.

When done right, an oligotype is a much more accurate proxy to an organism in an environment than an OTU or a taxon (within the limitations of the marker gene, amplified region, etc). For a better explanation/discussion, please take a look at [this publication](http://www.pnas.org/content/111/28/E2875.abstract).

### What makes oligotyping powerful?

Depending on the environment, an operational taxonomic unit identified through clustering analysis or classification may contain a very large number of distinct organisms. Oligotyping allows its user to detect and utilize very subtle nucleotide variations among the reads that were binned together in one OTU to create oligotypes that are not partitioned based on an arbitrary similarity thresholds.

### What is the weakness of the approach?

Oligotyping is a method that is most efficient when it is used to decompose very closely related taxa. This doesn’t necessarily mean oligotyping can be done only on genus level taxa or 3% OTUs. In some cases an entire phylum can be decomposed successfully (for instance Cyanobacteria phylum found in ocean samples), in other cases one genus could be too diverse for exhaustive explanation of diversity with oligotyping (such as Bacteroides genus in human gut).

This topic deserves more discussion, but I think it is generally safe to say that the current implementation of oligotyping is not the appropriate tool to be used when there is still diversity for taxonomy analysis or OTU clustering can resolve.

### I have other questions!

Please send them to meren at mbl dot edu, or write a comment below.

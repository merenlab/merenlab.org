---
layout: page
title: fasta [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/fasta
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/FASTA.png" alt="FASTA" style="width:100px; border:none" />

A FASTA-type anvi'o artifact. This artifact is typically provided **by the user** for anvi'o to import into its databases, process, and/or use.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-dereplicate-genomes](../../programs/anvi-dereplicate-genomes)</span> <span class="artifact-p">[anvi-script-fix-homopolymer-indels](../../programs/anvi-script-fix-homopolymer-indels)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-dereplicate-genomes](../../programs/anvi-dereplicate-genomes)</span> <span class="artifact-r">[anvi-search-palindromes](../../programs/anvi-search-palindromes)</span> <span class="artifact-r">[anvi-script-compute-ani-for-fasta](../../programs/anvi-script-compute-ani-for-fasta)</span> <span class="artifact-r">[anvi-script-fix-homopolymer-indels](../../programs/anvi-script-fix-homopolymer-indels)</span> <span class="artifact-r">[anvi-script-reformat-fasta](../../programs/anvi-script-reformat-fasta)</span></p>


## Description

A [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file that does not necessarily meet the standards of a <span class="artifact-n">[contigs-fasta](/software/anvio/help/7.1/artifacts/contigs-fasta)</span>. While it is not necessary for all programs, if a given anvi'o program requires a <span class="artifact-n">[contigs-fasta](/software/anvio/help/7.1/artifacts/contigs-fasta)</span>, the program <span class="artifact-n">[anvi-script-reformat-fasta](/software/anvio/help/7.1/programs/anvi-script-reformat-fasta)</span> can turn a regular fasta into a <span class="artifact-n">[contigs-fasta](/software/anvio/help/7.1/artifacts/contigs-fasta)</span> with the flag `--simplify-names`.

### What is a FASTA file?

A FASTA file typically contains one or more DNA, RNA, or amino acid sequences that are formatted as follows:

```
>SEQUENCE_ID VARIOUS_SEQUENCE_DATA
SEQUENCE
(...)
```

The line that starts with the character `>` is also known as the 'defline' for a given sequence. The `VARIOUS_SEQUENCE_DATA` region of the defline can be empty, or contain additional data such as the NCBI taxon ID, GI accession number, a text description of the sequence, or the start and end positions if the sequence is a portion of a larger sample. Because the FASTA file format was designed before there weren't even enough electronic calculators on the planet, there is no actual standard format to organize additional information shared in the defline.

The sequence itself is typically written in standard [IUPAC format](https://en.wikipedia.org/wiki/Nucleic_acid_notation), although you may find FASTA files with sequences that contain lower-case letter, mixed letters, no letters, or pretty much anything really. Over the years we have seen everything, and suggest you to take a careful look at your FASTA files before doing anything with them unless you generated them yourself.

You can learn more about the FASTA format on its [glorious Wikipedia page](https://en.wikipedia.org/wiki/FASTA_format).


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/fasta.md) to update this information.


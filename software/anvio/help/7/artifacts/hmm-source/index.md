---
layout: page
title: hmm-source [artifact]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/HMM.png" alt="HMM" style="width:100px; border:none" />

A HMM-type anvi'o artifact. This artifact is typically provided **by the user** for anvi'o to import into its databases, process, and/or use.

Back to the **[main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-script-pfam-accessions-to-hmms-directory](../../programs/anvi-script-pfam-accessions-to-hmms-directory)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-compute-completeness](../../programs/anvi-compute-completeness)</span> <span class="artifact-r">[anvi-delete-hmms](../../programs/anvi-delete-hmms)</span> <span class="artifact-r">[anvi-get-sequences-for-hmm-hits](../../programs/anvi-get-sequences-for-hmm-hits)</span> <span class="artifact-r">[anvi-run-hmms](../../programs/anvi-run-hmms)</span> <span class="artifact-r">[anvi-script-gen-hmm-hits-matrix-across-genomes](../../programs/anvi-script-gen-hmm-hits-matrix-across-genomes)</span> <span class="artifact-r">[anvi-script-get-hmm-hits-per-gene-call](../../programs/anvi-script-get-hmm-hits-per-gene-call)</span></p>


## Description

An HMM source is a collection of one or more hidden Markov models. HMM sources can be used to identify and recover genes in a <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> that match to those described in the model.

Models in a given HMM source can be searched in a given <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> via the program <span class="artifact-n">[anvi-run-hmms](/software/anvio/help/7/programs/anvi-run-hmms)</span> which would yield an <span class="artifact-n">[hmm-hits](/software/anvio/help/7/artifacts/hmm-hits)</span> artifact. An anvi'o installation will include multiple HMM sources by default. But HMMs for any set of genes can also be put together by the end user and run on any anvi'o <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span>.

HMM hits in a <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> for a given <span class="artifact-n">[hmm-source](/software/anvio/help/7/artifacts/hmm-source)</span> source will be accessible to anvi'o programs globally. Sequences that match to HMM hits can be recovered in an aligned or non-aligned fashion as <span class="artifact-n">[fasta](/software/anvio/help/7/artifacts/fasta)</span> files for downstream analyses including phylogenomics, they can be displayed in anvi'o interfaces, reported in summary outputs, and so on.

Running <span class="artifact-n">[anvi-db-info](/software/anvio/help/7/programs/anvi-db-info)</span> on a <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> will list HMM sources available in it.

### Default HMM sources

An anvi'o installation will include [multiple HMM sources](https://github.com/meren/anvio/tree/master/anvio/data/hmm) by default. These HMM sources can be run on any <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> with <span class="artifact-n">[anvi-run-hmms](/software/anvio/help/7/programs/anvi-run-hmms)</span> to identify and store <span class="artifact-n">[hmm-hits](/software/anvio/help/7/artifacts/hmm-hits)</span>:

<div class="codeblock" markdown="1">
<span class="artifact&#45;n">[anvi&#45;run&#45;hmms](/software/anvio/help/7/programs/anvi&#45;run&#45;hmms)</span> &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span>
</div>

The default HMM sources in anvi'o include:

* **Bacteria_71**: 71 single-copy core genes for domain bacteria that represent a modified version of the HMM profiles published by [Mike Lee](https://doi.org/10.1093/bioinformatics/btz188). The anvi'o collection excludes `Ribosomal_S20p`, `PseudoU_synth_1`, `Exonuc_VII_S`, `5-FTHF_cyc-lig`, `YidD` and `Peptidase_A8` occurred in Lee collection (as they were exceptionally redundant or rare among MAGs from various habitats), and includes `Ribosomal_S3_C`, `Ribosomal_L5`, `Ribosomal_L2` to make it more compatible with [Hug et al](https://www.nature.com/articles/nmicrobiol201648)'s set of ribosomal proteins.
* **Archaea_76**: 76 single-copy core genes for domain archaea by [Mike Lee](https://doi.org/10.1093/bioinformatics/btz188).
* **Protista_83**: 83 single-copy core genes for protists (domain eukarya) by [Tom O. Delmont](http://merenlab.org/delmont-euk-scgs).

Apart from these, anvi'o also includes a number of HMM profiles for individual ribosomal RNA classes derived from [Torsten Seemann's tool](https://github.com/tseemann/barrnap) (we split them into individual classes after [this](https://github.com/merenlab/anvio/issues/1411)):

* **Ribosomal\_RNA\_5S** (eukarya + archaea + bacteria; also includes 5.8S).
* **Ribosomal\_RNA\_12S** (mitochondria)
* **Ribosomal\_RNA\_16S** (bacteria + archaea + mitochondria)
* **Ribosomal\_RNA\_18S** (eukarya)
* **Ribosomal\_RNA\_23S** (bacteria + archaea)
* **Ribosomal\_RNA\_28S** (eukarya)

When <span class="artifact-n">[anvi-run-hmms](/software/anvio/help/7/programs/anvi-run-hmms)</span> is run on an anvi'o <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span> without providing any further arguments, it automatically utilizes all the default HMM sources.

{:.notice}
Similar to Ribosomal RNAs, anvi'o can also identify Transfer RNAs. Even though Transfer RNAs will also appear as an HMM source for all downstream analyses, their initial identification will require running <span class="artifact-n">[anvi-scan-trnas](/software/anvio/help/7/programs/anvi-scan-trnas)</span> program on a <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span>.

### User-defined HMM sources

The user can employ additional HMM sources to identify matching genes in a given <span class="artifact-n">[contigs-db](/software/anvio/help/7/artifacts/contigs-db)</span>.

Any directory with expected files in it will serve as an HMM source:

<div class="codeblock" markdown="1">
anvi&#45;run&#45;hmms &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7/artifacts/contigs&#45;db)</span> \
              &#45;&#45;hmm&#45;source /PATH/TO/USER&#45;HMM&#45;DIRECTORY/
</div>

Anvi'o will expect the HMM source directory to contain six files (see this for [an example directory](https://github.com/merenlab/anvio/tree/master/anvio/data/hmm/Protista_83)). These files are explicitly defined as follows:

* **genes.hmm.gz**: A gzip of concatenated HMM profiles. One can (1) obtain one or more HMMs by computing them from sequence alignments or by downloading previously computed ones from online resources such as [Pfams](https://pfam.xfam.org/family/browse?browse=new), (2) concatenate all profiles into a single file called `genes.hmm`, and finally (3) compress this file using `gzip`.
* **genes.txt**: A TAB-delimited file that must contain three columns: `gene` (gene name), `accession` (gene accession number (can be anything unique)), and `hmmsource` (source of HMM profiles listed in genes.hmm.gz). The list of gene names in this file must perfectly match to the list of gene names in genes.hmm.gz.
* **kind.txt**: A flat text file which contains a single word identifying what type of profile the directory contains. This information will appear in interfaces. Use a single, descriptive word for your collection.
* **reference.txt**: A file containing source information for this profile to cite it properly.
* **target.txt**: A file that specifies the target *alphabet* and  *context* that defines how HMMs should be searched (this is a function of the HMM source that is used). The proper notation is 'alphabet:context'. Alphabet can be `AA`, `DNA`, or `RNA`. Context can be `GENE` or `CONTIG`. The content of this file should be any combination of one alphabet and one context term. For instance, if the content of this file is `AA:GENE`, anvi'o will search genes amino acid sequences, and so on. An exception is `AA:CONTIG`, which is an improper target since anvi'o can't translate contigs to amino acid sequences. See [this](https://github.com/meren/anvio/pull/402) for more details. Please note that HMMs that target `DNA:CONTIG` will result in new gene calls in the contigs database to describe their hits.
* **noise_cutoff_terms.txt**: A file to specify how to deal with noise. [See this comment](https://github.com/merenlab/anvio/issues/498#issuecomment-362115921) for more information on the contents of this file.


### Creating anvi'o HMM sources from ad hoc PFAM accessions

It is also possible to generate an anvi'o compatible HMMs directory for a given set of PFAM accession ids. For instance, the following command will result in a new directory that can be used immediately with the program <span class="artifact-n">[anvi-run-hmms](/software/anvio/help/7/programs/anvi-run-hmms)</span>:

<div class="codeblock" markdown="1">
<span class="artifact&#45;n">[anvi&#45;script&#45;pfam&#45;accessions&#45;to&#45;hmms&#45;directory](/software/anvio/help/7/programs/anvi&#45;script&#45;pfam&#45;accessions&#45;to&#45;hmms&#45;directory)</span> &#45;&#45;pfam&#45;accessions&#45;list PF00705 PF00706 \
                                               &#45;O AD_HOC_HMMs
</div>

These IDs can be given through the command line as a list, or through an input file where every line is a unique accession id.

An example. Let's assume we have a genome or a metagenome that looks like this:

![PFAM example](../../images/p214-wo-upxz.png)

And we wish to identify locations of genes that match to this model: [http://pfam.xfam.org/family/PF06603](http://pfam.xfam.org/family/PF06603)

One can run this command:

<div class="codeblock" markdown="1">
<span class="artifact&#45;n">[anvi&#45;script&#45;pfam&#45;accessions&#45;to&#45;hmms&#45;directory](/software/anvio/help/7/programs/anvi&#45;script&#45;pfam&#45;accessions&#45;to&#45;hmms&#45;directory)</span> &#45;&#45;pfam&#45;accessions&#45;list PF06603 \
                                                &#45;O UpxZ
</div>

which would createa a directory called `UpxZ`. Then, one would run this command to find matches to this model in a given contigs database:

<div class="codeblock" markdown="1">
<span class="artifact&#45;n">[anvi&#45;run&#45;hmms](/software/anvio/help/7/programs/anvi&#45;run&#45;hmms)</span> &#45;c CONTIGS.db \
               &#45;H UpxZ/ \
               &#45;&#45;num&#45;threads 4
</div>

Now it is possible to get the sequences matching to this model:

<div class="codeblock" markdown="1">
<span class="artifact&#45;n">[anvi&#45;get&#45;sequences&#45;for&#45;hmm&#45;hits](/software/anvio/help/7/programs/anvi&#45;get&#45;sequences&#45;for&#45;hmm&#45;hits)</span> &#45;c CONTIGS.db \
                                 &#45;&#45;hmm&#45;source UpxZ \
                                 &#45;o UpxZ.fa
</div>

```
Contigs DB ...................................: Initialized: CONTIGS.db (v. 19)
Hits .........................................: 8 hits for 1 source(s)
Mode .........................................: DNA sequences
Genes are concatenated .......................: False
Output .......................................: UpxZ.fa
```

Or see where they are by visualizing the project using again:

<div class="codeblock" markdown="1">
<span class="artifact&#45;n">[anvi&#45;interactive](/software/anvio/help/7/programs/anvi&#45;interactive)</span> &#45;p PROFILE.db \
                  &#45;c CONTIGS.db
</div>

![PFAM example](../../images/p214-w-upxz.png)


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/hmm-source.md) to update this information.


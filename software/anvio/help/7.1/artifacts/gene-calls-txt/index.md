---
layout: page
title: gene-calls-txt [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/gene-calls-txt
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact can be generated, used, and/or exported **by anvi'o**. It can also be provided **by the user** for anvi'o to import into its databases, process, and/or use.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-export-gene-calls](../../programs/anvi-export-gene-calls)</span> <span class="artifact-p">[anvi-script-transpose-matrix](../../programs/anvi-script-transpose-matrix)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-script-transpose-matrix](../../programs/anvi-script-transpose-matrix)</span></p>


## Description

This file describes all of the gene calls contained in a <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> from a specified list of sources. It is the output of <span class="artifact-n">[anvi-export-gene-calls](/software/anvio/help/7.1/programs/anvi-export-gene-calls)</span>. 

For each gene identified, this file provides various information, including the caller ID, [start and stop position](http://merenlab.org/software/anvio/help/artifacts/external-gene-calls/#gene-startstop-positions), direction, whether or not the gene is partial, the [call type](http://merenlab.org/software/anvio/help/artifacts/external-gene-calls/#call-type), source and version (if available ), and the amino acid sequence. 

{:.notice}
Want more information? This file is in the same format as an <span class="artifact-n">[external-gene-calls](/software/anvio/help/7.1/artifacts/external-gene-calls)</span>, so check out that page. 

Here is an example from the Infant Gut Dataset: 

    gene_callers_id    contig              start    stop    direction    partial    call_type    source     version   aa_sequence
    0                  Day17a_QCcontig1    0        186     f            1          1            prodigal   v2.60     GSSPTAGVEQKQKPTWFLLFLFYSLFFDKLEEGTLKTFIRLKGSYRRMNTSNFSYGIMCLL
    1                  Day17a_QCcontig1    214      1219    f            0          1            prodigal   v2.60     MKILLYFEGEKILAKSGIGRALDHQKRALSEVGIEYTLDADCSDYDILHINTYGVNSHRMVRKARKLGKKVIYHAHSTEEDFRNSFIGSNQLAPLVKKYLISLYSKADHLITPTPYSKTLLEGYGIKVPISAISNGIDLSRFYPSEEKEQKFREYFKIDEEKKVIICVGLFFERKGITDFIEVARQLPEYQFIWFGDTPMYSIPKNIRQLVKEDHPENVIFPGYIKGDVIEGAYAAANLFFFPSREETEGIVVLEALASQQQVLVRDIPVYQGWLVANENCYMGHSIEEFKKYIEGLLEGKIPSTREAGYQVAEQRSIKQIGYELKEVYETVLS
    2                  Day17a_QCcontig1    1265     2489    f            0          1            prodigal   v2.60     MKIGFFTDTYFPQVSGVATSIKTLKDELEKHGHEVYIFTTTDPNATDFEEDVIRMPSVPFVSFKDRRVVVRGMWYAYLIAKELELDLIHTHTEFGAGILGKMVGKKMKIPVIHTYHTMYEDYLHYIAKGKVVRPSHVKFFSRVFTNHTTGVVCPSERVIEKLRDYGVTAPMRIIPTGIEIDKFLRPDITEEMIAGMRQQLGIEEQQIMLLSLSRISYEKNIQAIIQGLPQVIEKLPQTRLVIVGNGPYLEDLKELAEELEVSEYVQFTGEVPNEEVAIYYKAADYFVSASTSETQGLTYTEAMAAGVQCVAEGNAYLNNLFDHESLGKTFKTDSDFAPTLIDYIQANIKMDQTILDEKLFEISSTNFGNKMIEFYQDTLIYFDQLQMEKENADSIKKIKVKFTSLRK
    ...



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/gene-calls-txt.md) to update this information.


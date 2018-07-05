---
layout: page
title: Research
modified: 2015-02-02T20:53:07.573882-04:00
comments: false
excerpt: "An overview of research interests and directions in our lab."
image:
   thumb: holistic-flow.png
---

{% include _toc.html %}

<blockquote>
Study hard what interests you the most, in the most undisciplined, irreverent, and original manner possible.


<div class="blockquote-author">Richard Feynman</div>
</blockquote>

<div style="height: 350px; width: 300px; float: right; padding-left: 20px;">
<img src="images/all.png" style="border:none;" />
</div>

The vast majority of life on our planet is *microbial*. An astonishing number of microbial organisms living in terrestrial and marine habitats represent a biomass that exceeds every living organism that can be seen by naked eye, combined.

The inconceivable diversity of microbes allow them to synthesize or break down a wide array of chemical substrates, and govern biogeochemical cycles that make Earth a habitable planet for much less talented organisms (such as ourselves). Our own body is also home to a diverse assemblage of microbial cells. Bacteria that colonize our gastrointestinal tract help us maintain our health by extracting energy from undigested carbohydrates, synthesizing vitamins, and metabolizing xenobiotics.

Microbes are essential to Earth’s functioning at every scale, and understanding them is imperative for a complete understanding of life.

In our lab, we combine our expertise and interest in microbiology and computation to investigate the diversity and functioning of microbial communities in environments ranging from the human gastrointestinal tract and oral cavity, to sewages, oceans, and soils.

We design experiments, develop software platforms, and work with high-throughput sequencing data for marker genes, metagenomes, and metatranscriptomes to get closer to understanding mechanisms by which microbes interact with their surroundings, evolve, disperse, and initiate and/or adapt to environmental change.

We believe in interdisciplinary, hypothesis-driven, open and collaborative science, and we strive to identify the best computational and experimental practices that follows the inspiration from biology.

The following is an incomplete list of current or future interests of our group.

<div style="clear:both"></div>

## Studying fundamentals of microbial life and evolution

Despite their great importance for the habitability of the planet and maintenance of our health, our understanding of the diversity, functioning, and evolution of microbial life is far from being complete.

This is partially due to our inability to bring microbial life into the lab environment for comprehensive, conventional investigations: in most cases, it is very challenging to successfully isolate individual microbes from their interactive and complex environments and keep them functionally alive in a controlled setting. Even when isolation is possible, understanding how well the isolated members of microbial populations represent their environmental population is not necessarily always clear. Thankfully, in parallel to conventional approaches, we can use state-of-the-art molecular and computational techniques to recover the genomic content of naturally occurring microbes directly from the environment, and investigate some of the most fundamental aspects of their life and evolution. Our lab often uses high-throughput sequencing data to identify known and novel populations of microbes, determine their distribution, functional potential, and genomic heterogeneity to study their relationships with other microbes and environments they reside, which can inform us about their evolution.

Microbial ecology and evolution is happening everywhere, and investigating fundamentals of microbial life is independent of habitats. Hence, we are as happy in our lab working with marine environments as we are working with the human guts, or insect ovaries.

<div style="width: 300px; float: right; padding-left: 20px;">
<a href="http://i.imgur.com/HrlXPOF.jpg"><img src="http://i.imgur.com/HrlXPOF.jpg" style="border:none;" /></a>
</div>

For instance, [one of our recent studies focused on SAR11](https://doi.org/10.1101/170639), which is one of the most abundant microbial lineages in marine habitats. In collaboration with [Stephen Giovannoni](http://microbiology.science.oregonstate.edu/dr-stephen-giovannoni) from the **Oregon State University**, [Mike Rappé](https://rappelab.wordpress.com/) from the **University of Hawai'i**, and [Ismail Uysal](https://www.researchgate.net/profile/Ismail_Uysal3) from the **University of South Florida**, we developed and applied novel approaches to investigate non-synonymous genomic heterogeneity that emerges in a *single SAR11 population* as it goes around the world through large oceanic current movements. Using single-amino acid variants, we observed significantly more protein variants in cold currents and an increased number of protein sweeps in warm currents, exposing a global pattern of alternating genomic diversity for this SAR11 population. From the geographic partitioning of SAAVs we could suggest that natural selection, rather than neutral evolution, is the main driver of the evolution of SAR11 in surface oceans. We will continue to explore approaches to link environmental heterogeneity to the genomic context and protein structure, and generate testable hypothesis from our 'omics studies to better understand determinants of microbial fitness and evolution.


## Metagenomic and experimental insights into Fecal Microbiota Transplantation experiments

Fecal Microbiota Transplantation (FMT) experiments, the quite literal transference of fecal matter from a healthy donor to a recipient, gained recognition as an effective and relatively safe treatment for recurrent or refractory *Clostridium difficile* infection (CDI). Its success in treating otherwise very hard to eradicate CDI sparked great interest in investigating FMT as a treatment also for other medical conditions associated with intestinal dysbiosis, such as ulcerative colitis, Crohn’s disease, irritable bowel syndrome, and even metabolic syndrome or neurodevelopmental and autoimmune disorders.

Despite the excitement due to its therapeutic potential, FMTs also present challenges for researchers and clinicians with potential adverse outcomes, including the transfer of infectious organisms or contaminants from the environment. Deeper insights into FMTs from a basic science perspective will likely help us address these caveats, while teaching us a lot about the microbial dynamics. In fact, in our lab we see FMTs as excellent tools to study the functional basis of human gut colonization and the microbial response to environmental stress. As we learn about microbes by studying FMTs, we are hoping to get closer to understand key microbial populations that are responsible for beneficial outcomes of this procedure.

<div style="width: 300px; float: right; padding-left: 20px;">
<a href="{{ site.url }}/data/2017_Lee_et_al_FMT/figure-01.png"><img src="{{ site.url }}/data/2017_Lee_et_al_FMT/figure-01.png" style="border:none;" /></a>
</div>

[In a recent pilot study](https://doi.org/10.1186/s40168-017-0270-x), we collected fecal samples from two recipients before and after an FMT procedure. By using our metagenomic approaches, we could track individual donor microbial populations in a very highly resolved manner in the recipient guts to study microbial colonization dynamics. Some donor microbes that were initially absent from the recipient guts colonized both individuals successfully, and we could detect them even 8 weeks after the transfer. In contrast, some other donor microbes failed to colonize neither of our recipients. Most interestingly, we found that the colonizers we could identify in our pilot study were common among a much larger cohort of healthy people, and non-colonizers occurred rarely in the same cohort. We use these crucial insights to investigate the functional potential of microbes with high colonization properties, inform cultivation experiments in our lab, and look forward to generate hypotheses that we can translate to model systems.

<div style="clear:both"></div>


## Advanced open-source software platforms for high-resolution microbial 'omics

Computation is at the core of every scientific discipline. The fields of microbiology and microbial ecology, which rely on *big data* more and more, have  also dramatically benefited from the advances in computation during the last decade. The importance of computation in life sciences puts our lab in a lucky situation, but we try to use our skills in computation wisely.

How much of the scientific questions we dare to ask depend on the availability of computational solutions that can facilitate the investigation of those questions in complex datasets we plan to generate? Although the inherent link between the tool and thinking will continue to bind the two together, we believe it must be mostly the intellectual curiosity what drives the direction of science, and not the comfort of what is available.

The agreement we have with ourselves in our lab is to keep the biology as the sole inspiration of our direction, and never let the computational conveniences find their ways into our thinking in the expense of our ability to explore fundamental questions. We strive to create software that would allow users to [get their hands dirty]({% post_url anvio/2015-12-09-musings-over-commamox %}) with their data, without imposing boilerplate analysis practices. In most cases the questions we are interested in require precise computational approaches which can offer enough resolution that would enable us to detect subtle changes. These needs resulted in various software solutions we proposed, including [oligotyping](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12114/full), and [minimum entropy decomposition](http://www.nature.com/ismej/journal/v9/n4/full/ismej2014195a.html) for the analysis of marker gene data, and [anvi’o](https://peerj.com/articles/1319/) for comparative genomics, metagenomics, metatranscriptomics, and visualization of complex data.

We intend to maintain our flexibility, and let the incoming questions shape and re-shape our software.

## Holistic approaches to investigate complex gastrointestinal diseases

Inflammatory bowel diseases (IBD) describe a number of prolonged inflammatory conditions of the human colon and small intestine that affect increasingly more people. Substantial evidence links the occurrence of these chronic, relapsing conditions to aberrant immune responses to microbes that colonize the gut. Although there are remarkable shifts in microbial communities with the presence of IBD, there is no evidence for microbial members or metabolisms that are specific to the guts of individuals who suffer from IBD, and that are absent from healthy guts. The current understanding of mechanisms that lead to the development of these conditions is unfortunately limited. However, there is room for improvement in ways these diseases are studied. In collaboration with [Gene Chang’s group](http://changlab.uchicago.edu/) at the **University of Chicago**, [Mitch Sogin](https://scholar.google.com/citations?user=OQbV_3UAAAAJ)’s and [Hilary Morrison](https://scholar.google.com/citations?user=-cmjxwsAAAAJ)’s groups at the **Marine Biological Laboratory**, and [Dionysios Antonopoulos](http://www.anl.gov/contributors/dionysios-antonopoulos)’ group at the **Argonne National Laboratory**, we study the microbial members of IBD patients and the evolution of these microorganisms through genome-resolved metagenomics using longitudinal sampling strategies that allow individuals to serve as their own controls. Our holistic approach includes combining shotgun metagenomics with cultivation or organisms of interest, associating our findings with host factors, and generating hypotheses to be tested in model systems.


## Genome-resolved understanding of the human oral cavity

The oral cavity represents the receiving end of our digestive tract where the processing of food begins. Like every other mucosal surface on our body, bacteria also colonize the oral cavity, and they play a critical role in health and disease states of the mouth. Some medical conditions in the oral cavity, such as tooth decay, gum diseases, root canal infections, and tonsillitis can result in systemic diseases. Hence, a complete microbial understanding of this environment has always been essential for medical reasons. Besides its immediate relevance for overall health, we believe the oral cavity represents a fascinating environment to study the ecology of microbes. Due to the lack of any physical barriers, and continuous flow of saliva, microbes can disperse everywhere in the mouth. However, their distribution is far from random. While the microbial occupants of niches in the human mouth (such as tongue, cheeks, gums) form distinct communities, a universe of interactions emerges in this relatively small environment. We [previously showed](http://www.pnas.org/content/111/28/E2875) the differential distribution of very closely related microbial organisms in different oral sites at the marker gene-level using our high-resolution computational approaches. Today, in collaboration with [Jessica Mark Welch](http://www.mbl.edu/jbpc/staff/markwelchj/) from the **Marine Biological Laboratory**, we go further in an attempt to develop a genome-resolved understanding of the oral cavity, and characterize the distribution of pangenomic traits across oral sites.

## Studying public health through the guts of the urban ecosystem

Proper removal of waste is one of the most basic requirements of settled human communities. In fact, we owe our ability to live in such small geographical areas with such high population densities primarily to sewer infrastructures and their ability to effectively evacuate and treat human waste. Today, a modern sewer infrastructure is a critical component of every city. With its ubiquitous parts and components (i.e., toilets, pipes, drains, manholes, pumping stations, and treatment centers), this built environment represents a new, and not yet well-characterized ecosystem for microbial life. Although the functioning of microbes in wastewater treatment efforts has been extensively studied due to their industrial applications for bioremediation, understanding the microbial life in the rest of the sewer infrastructure, especially the pipe systems, has not been a major area of interest. In a recent study led by [Ryan Newton](https://scholar.google.com/citations?user=WaVCQwgAAAAJ) and [Sandra McLellan](http://home.freshwater.uwm.edu/mclellanlab/) from the **University of Wisconsin-Milwaukee**, we [demonstrated](http://mbio.asm.org/content/6/2/e02574-14.short) that it is possible to predict the level of obesity in a given US city with more than 80% accuracy by only analyzing the microbial community signatures found in sewage samples. The link between the microbial community structure and the level of obesity as demonstrated by this finding suggests a potentially very important role for sewages to track public health. We pursue a deeper understanding of the ecology of the sewer ecosystems through marker genes and shotgun metagenomes in an attempt to develop baseline metrics for microbial signatures that can identify matters of public health, and environmental change.

---

## Our 2 cents

Some vocabulary we try to use, and promote as much as we can:

- "**Bacteria and archaea**” to describe the two ~~major~~ famous ~~domains~~ whatevers of life, *instead of* "**prokaryotes**”. Because major scientific advances of recent times should not be ignored by scientists. Here is an opinion piece on this from Norman Pace: ["It's time to retire the prokaryote"]({{site.url}}/files/time-to-retire-prokaryote.pdf).

- "**Single-nucleotide variant**” (SNV) to describe nucleotide positions with variation that emerges from the mapping of short reads from environmental shotgun metagenomes to a genomic context, *instead of* "**single-nucleotide polymorphism**” (SNP). [SNP](https://en.wikipedia.org/wiki/Single-nucleotide_polymorphism) has a very specific definition, which makes it inappropriate to use in the context of metagenomic mapping.

- "**Metagenome-assembled genomes**" (MAGs) to describe bins of assembled contigs from shotgun metagenomic data, *instead of* "**draft genomes**". Because the term draft genome have long been used to describe not-yet-finalized genomes acquired from the whole genome sequencing and assembly of cultured microbial isolates.

- "**Marker gene amplicons**" to describe the high-throughput sequencing data of marker genes targeted by (universal) primers, *instead of* "**metagenomes**". Metagenomics is broadly defined as the study of genetic material directly recovered from a sample. With some effort, this broad definition *could* include marker gene amplicon surveys, but it really should not. The authors who coined the term associate the 'metagenome' with having access to the *[collective genomes of microbes](http://www.sciencedirect.com/science/article/pii/S1074552198901089)* in an environment. Today the term metagenomics is mostly used to describe [shotgun sequencing](https://en.wikipedia.org/wiki/Shotgun_sequencing) of the environmental DNA in order to explore the functional potential or community composition of a given community at the level of metagenomic short reads (without assembly), assembled contiguous DNA segments (after assembly), or metagenome-assembled genomes (after assembly and binning). Marker gene surveys, including the ones that amplify hyper-variable regions of the ribosomal RNA genes, does not fit into what metagenomics describes, and the misuse of the term creates a lot of confusion.

- "**kbp**", "**Mbp**", or "**Gbp**" to communicate the number of base pairs in a contig, or a genome, *instead of* "**KB**", "**MB**", or "**GB**". Because the latter are commonly used to describe the amount of digital information in 'bytes', and the alternative use is not appropriate. Although it appears in the literature quite often, the use of "**Kb**", "**Mb**", or "**Gb**" are not very suitable either, since these units are commonly used to quantify digital information in 'bits'.

- *Anything else* to describe microbes that are not yet characterized, *instead of* "**microbial dark matter**". [Just stop it]({% post_url miscellaneous/2017-06-22-microbial-dark-matter %}). [Please](https://twitter.com/merenbey/status/907708177403322368).


If you believe we need corrections, please don't hesitate to [write to us]({{site.url}}/people/).


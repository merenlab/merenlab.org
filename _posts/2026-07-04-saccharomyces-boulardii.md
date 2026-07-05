---
layout: post
authors: [meren]
title: "A probiotic yeast, colorectal cancer, and the gut microbiome"
excerpt: "A rare correlative microbiome study careful enough to think with, and an arrow or two I would redraw."
modified: 2026-07-04
tags: [publication, human-microbiome]
categories: [miscellaneous]
comments: true
noleftpanel: false
notoc: false
thumbnail: /images/thumbnails/2026-07-04-saccharomyces-boulardii.png
image:
  feature: https://merenlab.org/images/thumbnails/2026-07-04-saccharomyces-boulardii.png
  display: false
---

I finally had a chance to read this very recent gut microbiome / colorectal cancer study, "[*Multi-omics analysis of saccharomyces boulardii supplementation reveals coordinated microbiome, metabolic, and immune signaling changes accompanying tumor suppression*](https://doi.org/10.1080/19490976.2026.2690687)" published by [Troels Holger Vaaben](https://scholar.google.com/citations?user=tZbZawMAAAAJ&hl=en), [Morten Otto Alexander Sommer](https://scholar.google.com/citations?user=r2ZoOsMAAAAJ&hl=en), and others at [The Novo Nordisk Foundation Center for Biosustainability](https://www.biosustain.dtu.dk/), [DTU](https://dtu.dk/), and wanted to share my thoughts about it.

This all started as a LinkedIn post, but when I realized that I was going over their limit by over like 4,000 characters, I came here to the blog.

The post is organized the way my LinkedIn post was organized, but I managed to make my arrows much prettier (you will see what I mean by that), and divided the text into sections with the hope that it will be easier to read.

## What the study did

Vaaben and colleagues take an immunocompetent mouse model of colorectal cancer, and feed it *Saccharomyces boulardii*, one of the most widely used yeast probiotics that transiently colonizes the human gastrointestinal tract, and show that the growth of distant tumors slows down even though the yeast itself does not colonize the tumors themselves.

In addition to their extensive omics analyses that help them make this point, they also have data to show correlated changes in the plasma metabolome, in circulating host cytokines, and the host tumor transcriptome.

They are extremely careful to frame all this as *correlative* rather than *mechanistically causal*. A level of care and attention that prevents me from being able to NOT think that if the people at the InnoHK Microbiota I-Center had this kind of data and findings, they would have shoved not only the *S. boulardii*, but the entire mice down the throats of people (here is a [link](https://merenlab.org/blog/unfalsifiable-by-design/) for those of you who don't see where this is coming from).

So I REALLY appreciate the scientific rigor on display here by the authors, and postpone any of my complaints about the fact that [PRJNA1313482](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA1313482) returns an empty page as of this evening when it should have contained the raw sequencing data according to the Data Availability statement in the paper.

## Metabolic independence hypothesis as a tool

The reason this study first caught my eye was due to their very creative use of the 'metabolic independence' hypothesis, which I can summarize as "microbial systems in stress will be dominated by microbes that do not need others for their core metabolic activities".

This emerged a result of our observations in two separate studies we published in 2023 ([Watson, Füssel, Veseli et al](https://doi.org/10.1186/s13059-023-02924-x)) and in 2025 ([Veseli et al](https://doi.org/10.7554/eLife.89862)). In these studies we observed tha one of the hallmarks of stressful gut environments is the following: most microbes that remain in stressful habitats are those that are able to synthesise their own nucleotides, amino acids, and vitamins, while those that are reliant on the ecosystem for these kinds of essential nutrients due to their lacking genomic potential are unable to survive stress. We posit that this neutral phenomenon is the underlying reason for substantial changes in biodiversity in disease states, without any need to make those who survive these environments responsible for the emergence or trajectory of the disease. In a way, this puts microbes back to where they belong: the passenger seat. At least as the primary null hypothesis to be beaten in studies that ascribe them any causal role in non-communicable diseases of complex etiology.

The authors take metabolic independence hypothesis and use it in reverse, which is what I found very creative. Their reasoning is that if *S. boulardii* is handing out amino acids, it should relax the selective pressure favoring high metabolic independence (HMI) populations and allow low metabolic independence (LMI) populations to expand. Which would have implicated *S. boulardii* in the increasing representation of the members of the microbial community that are associated with homeostasis, through which it would implicate *S. boulardii* in other processes that influence tumors (we will get to this later).

To test that, they reconstruct over 600 high-quality genomes from mouse gut metagenomes, and estimate the metabolic independence score of each MAG using {% include PROGRAM name="anvi-script-estimate-metabolic-independence" %} with the default 33-module panel, and classify each genome as LMI or HMI. Their analysis shows that the per-sample percentage of the genomes classified as 'LMI' in the animals treated with *S. boulardii* changes significantly (p < 0.001) but modestly (which will be important later) as this proportion goes from 68% (placebo-treated) to 77% (*S. boulardii* treated).

The authors reason that the increase in LMI populations is likely due to the provisioning of amino acids by *S. boulardii* (which they [demonstrated previously](https://doi.org/10.1093/ismejo/wrae212): *S. boulardii* indeed donates amino acids to its environment, at least in co-cultures): as the members of the microbial community no longer need to be able to make their own amino acids to survive, those who can't make them start to expand in the gut.

## Some arrows of directionality (just to speculate)

As I praised at the beginning of this post, the authors do not go too far from here as they are extremely careful. They literally write "*the multi-omics integration identifies coordinated cross-compartment responses but does not establish directionality or mechanistic hierarchy*". They explicitly mention that all options for hierarchy and questions such as *what affects what* and *in what order that effect takes shape* are all on the table.

But there is nothing wrong with speculating a little after going this far. So when I try very hard to see how the authors would have speculated, I find their most *directional* sentence in the paper, "*transient microbiome reshaping [via S. boulardii treatment] can initiate systemic immunometabolic changes that influence tumor progression*", which points toward the following model:


```
S. boulardii -> microbiome -> host -> tumor
```

In this caricaturish view (for which I apologize), (1) *S. boulardii* reshapes the microbial community towards one that looks more like the microbial communities observed in homeostasis (with the enrichment of LMI populations) through the addition of amino acids, and (2) the *reshaped community* drives the systemic changes on the host side.

The first one is an ecological statement that makes sense, and the follow up is a logical extension of it, *if* we agree that microbiome in homeostasis can affect tumor progression.

But I think there is more to this. And here is where my intuition diverges from theirs: in my opinion, the relative depletion of LMI populations from the cancer environment is extremely unlikely to be only due to the lack of amino acids. As in, the depletion of populations that are unable to synthesise their own stuff in a host system is unlikely to be primarily driven by the depletion of the stuff they are unable to synthesise in the environment. Because if a collapsed public-goods pool were the sole cause of their demise, then simply topping up amino acids would have helped the LMI populations expand in unhealthy hosts. The porportion of LMI populations is merely a marker for a broader phenomenon that underlies the ecosystem beyond the nutrients it can provide. To dig a little further, let's think about the following questions:

* Would replacing *S. boulardii* with amino acid supplements recapitulate the same effects on the microbiome? I don't think so. If amino acid scarcity were the only barrier, LMI populations should not have been almost 70% of the untreated gut. There is a reason for their modest increase, and it is important to think about, in my opinion. 

* Would extending LMI populations in these animals from 68% to 77% on average without *S. boulardii* recapitulate the same effects on the tumors? I really don't think so. If we were talking about 0% to 9%, it would have been a much more interesting change to contend with comapred to a modest shift from 68% to 77%. 9% change makes it unlikely that this particular axis is what drives the changes on the host side.

None of this ignores the likely impact of *S. boulardii* on the microbial community through its amino acid provision, but I think these data assigns a larger role for *S. boulardii* on the system than that when the plasma metabolome and circulating cytokines data as well as the physiology of tumors are included in the picture. That role likely plays through the regulation of immunomodulatory activities, which likely exerts a larger influence on the gut microbiome than the amino acid provision, and likely contributes more to the expansion of LMI populations.

## Arrows I would draw (to speculate even further)

So instead of this linear relationship as a speculation,

```
S. boulardii -> microbiome -> host -> tumor
```

I would speculate that the interactions are likely a multifaceted process and I would model them as the following:

```
                      .-> tumor
             (1)     /
             .-> host
            /        \
S. boulardii          `-> microbiome <-.
            \                           \
             \                          /
              \  (2)                   /
               `----> environment ----`
```

Probably we all would name these arrows similarly regardless of how much we can explain the mechanisms.

The first arrow, (1), is perhaps operates through immunomodulation, which affects both tumors *and* the microbiome ('how' is another question which we don't know). The second named arrow, (2), is likely operates through amino acids, which changes yet another dimension of the environment, and affects the same microbiome. The microbiome in this model is primarily (re)shaped by (1), and the effect is sped up by (2).

Careful readers will likely see that there are no more arrows that connect gut microbes and the tumor in my argument here, which is clearly more of a thought experiment and a hypothesis based on my intuition of the role and relevance of the microbiome than anything else. You are more than welcome to beat it. I would love to change my mind.

## Final words

Of course, as the authors mention it in their work, how *transitionary* any of these effects are is a critical question here. I was much more strict when I was younger, but now I know that transitionary is not necessarily bad, and a lack of mechanistic insights should not be a reason to not explore practical applications as long as all the science is transparent, and we are not walked around by specific interest groups that seem to work against the priorities of the society and science itself.

Troels Holger Vaaben and colleagues put together a comprehensive evaluation of their system (with some limitations of course, such as the lack of tumor-free controls for starters, but I also would argue that 'healthy' controls in these model systems often serves a straw man as using the unhealthy animals as their own controls can be much more rigorous). They used bacteria as a marker of a larger shift in their system, and carefully considered the implications of the yeast that certainly is not a *non-player character* here.

What makes this study much more interesting than a gazillion others I see in the human microbiome field is the fact that it has enough dimensionality and care that inspires deeper thought than just a nod or shrug.

When authors, after generating ton of data and invest hours and hours and hours of work, tell their readers "*we can't say more, or one way or the other*", it is the most genuine of invitations for everyone to think more. It creates space to think more. Makes room for more narratives. Allows additional ideas to join in the fight, rather than going against unnecessarily strong narratives. And I thank them for it.

I am also very proud that our work could support such elaborate science, and I hope the authors could forgive me for saying so much about their work without consulting with them first.
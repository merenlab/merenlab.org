---
layout: post
authors: [meren]
title: "A probiotic yeast, colorectal cancer, and the gut microbiome"
excerpt: "A rare correlative microbiome study careful enough to think with, and an arrow or two I would redraw."
modified: 2026-07-04
tags: [cancer, human microbiome]
categories: [miscellaneous]
comments: true
noleftpanel: false
notoc: false
thumbnail: /images/thumbnails/2026-07-04-saccharomyces-boulardii.png
image:
  feature: https://merenlab.org/images/thumbnails/2026-07-04-saccharomyces-boulardii.png
  display: false
---

This all started as a LinkedIn post, but when I realized that I was going over like 4,000 characters over their limit, I came to our blog.

I finally had a chance to read this very recent gut microbiome / cancer study by Troels Holger Vaaben and colleagues, "Multi-omics analysis of saccharomyces boulardii supplementation reveals coordinated microbiome, metabolic, and immune signaling changes accompanying tumor suppression", and wanted to share my thoughts about it to prove you wrong in case you thought that Meren is someone who complains about everything.

The rest of this post is organized the way my LinkedIn post was organized, but I managed to make my arrows much prettier (you will see what I mean by that), and divided the text into sections with the hope that it will be easier to read.

## What the study did

[Troels Holger Vaaben](https://scholar.google.com/citations?user=tZbZawMAAAAJ&hl=en), [Morten Otto Alexander Sommer](https://scholar.google.com/citations?user=r2ZoOsMAAAAJ&hl=en), and others at [The Novo Nordisk Foundation Center for Biosustainability](https://www.biosustain.dtu.dk/), [DTU](https://dtu.dk/), take an immunocompetent mouse model of colorectal cancer, and feed it *Saccharomyces boulardii*, one of the most widely used yeast probiotics that transiently colonizes the human gastrointestinal tract, and show that the growth of distant tumors slows down even though the yeast itself does not colonize the tumors themselves.

In addition to their extensive omics analyses that help them make this point, they also have data to show correlated changes in the plasma metabolome, in circulating host cytokines, and the host tumor transcriptome.

They are extremely careful to frame all this as *correlative* rather than *mechanistically causal*. A level of care and attention that prevents me from being able to NOT think that if the people at the InnoHK Microbiota I-Center had this kind of data and findings, they would have shoved not only the *S. boulardii*, but the entire mice down the throats of people (here is a [link](https://merenlab.org/blog/unfalsifiable-by-design/) for those of you who don't see where this is coming from).

So I REALLY appreciate the scientific rigor on display here by the authors, and postpone any of my complaints about the fact that [PRJNA1313482](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA1313482) returns an empty page as of this evening when it should have contained the raw sequencing data according to the Data Availability statement in the paper (wink wink).

## Metabolic independence hypothesis as a tool

The reason this study first caught my eye was due to their very creative use of the 'metabolic independence' hypothesis, which I can summarize as "microbial systems in stress will be dominated by microbes that do not need others for their core metabolic activities", as a result of our observations in two separate works we published in 2023 ([Watson, Füssel, Veseli et al](https://doi.org/10.1186/s13059-023-02924-x)) and in 2025 ([Veseli et al](https://doi.org/10.7554/eLife.89862)). Vaaben et al take the metabolic independence hypothesis and use it to track the level of stress in these changing microbiomes by connecting a functional observation with an ecological hypothesis.

One of the hallmarks of stressful environments we observed in the two studies I mentioned above was the following: most of the microbes that remained in stressful habitats were those that were able to synthesise their own nucleotides, amino acids, and vitamins and so on, and all those that were reliant on others for these kinds of essential nutrients due to their lacking genomic potential were unable to survive stress (which is the underlying reason for substantial changes in biodiversity in disease states having to make those who survive these environments responsible for the emergence or trajectory of the disease, putting microbes back to where they belong: the passenger seat as the null hypothesis to be beaten).

The authors use this logic in reverse, which is what I found very creative. Their reasoning is that if *S. boulardii* is handing out amino acids, it should relax the selective pressure favoring high metabolic independence (HMI) populations and allow low metabolic independence (LMI) populations to expand. Which would have implicated *S. boulardii* in the return of the microbial community members that are associated with homeostasis, and would implicate this process in tumor dynamics.

To test that, they reconstruct over 600 high-quality genomes from mouse gut metagenomes, and score each MAG for metabolic independence using {% include PROGRAM name="anvi-script-estimate-metabolic-independence" %} with the default 33-module panel, and classify each genome as LMI or HMI. Which shows that the per-sample percentage of the genomes classified as 'LMI' in the animals treated with *S. boulardii* changes modestly (which will be important later) but significantly (p < 0.001) as it goes from 68% (placebo-treated) to 77% (*S. boulardii* treated).

The authors reason that the increase in LMI populations is likely due to the provisioning of amino acids by *S. boulardii* (which they [demonstrated previously](https://doi.org/10.1093/ismejo/wrae212): *S. boulardii* indeed donates amino acids to its environment, at least in co-cultures). As the members of the microbial community no longer need to be able to make their own amino acids to survive, those who can't make them start to appear.

## Some arrows of directionality (to speculate)

As I praised at the beginning of this post, the authors do not go too far from here as they are extremely careful. They literally write "*the multi-omics integration identifies coordinated cross-compartment responses but does not establish directionality or mechanistic hierarchy*". They explicitly mention that all options for hierarchy and *what affects what* and *in what order* are all on the table. But there is nothing wrong with speculating after going this far, so when I try very hard to see how the authors would have speculated, I find their most *directional* sentence, which points toward the following model as they write "*transient microbiome reshaping can initiate systemic immunometabolic changes that influence tumor progression*":


```
S. boulardii -> microbiome -> host -> tumor
```

In this caricaturish view (for which I apologize), *S. boulardii* reshapes the community towards one that looks more like those that are in homeostasis (with the enrichment of LMI populations) through the addition of amino acids, and the *reshaped community* drives the systemic changes on the host side.

The first one is an ecological statement that makes sense, and the follow up is a logical extension of it if we agree that microbiome in homeostasis can affect tumor progression. But I think there has to be more to this.

Why?

Let me try to explain: in my opinion, the relative depletion of LMI populations from the cancer environment is extremely unlikely to be only due to the lack of amino acids. As in, the depletion of populations that are unable to synthesise their own stuff in a host system is unlikely to be primarily driven by the depletion of the stuff they are unable to synthesise in the environment. Because if a collapsed public-goods pool were the sole cause of their demise, then simply topping up amino acids would have helped the LMI populations expand in unhealthy hosts. Right? So we can perhaps ask the following, and try to answer them with what is already available:

* Would replacing *S. boulardii* with amino acid supplements recapitulate the same effects on the microbiome? I don't think so -- because if amino-acid scarcity were the only barrier, LMI populations shouldn't already be the majority (68%) in the untreated gut. There is a reason for the modest increase, and it is important to think about, in my opinion. 

* Would extending LMI populations in these animals from 68% to 77% on average without *S. boulardii* recapitulate the same effects on the tumors? I really don't think so. If we were talking about 0% to 9%, I may have agreed. But a modest shift from 68% to 77% makes it unlikely that this particular axis -- the fraction of LMI populations -- is what drives the changes on the host side (even if the community is reshaped in other ways, too).

Two caveats are worth stating plainly here. First, this study has no healthy, tumor-free control -- both groups carry tumors -- so the 'depletion' of LMI populations in the cancer state is something I am importing from our earlier work rather than something these data show directly. Second, the tumor here sits under the skin rather than in the gut, so the gut is not the diseased tissue, which if anything makes the host -> microbiome direction *more* plausible than a gut-local tumor loop.

There is no need to ignore the likely impact of *S. boulardii* on the microbial community through its amino acid provision, but I think these data assigns a larger role for *S. boulardii* on the system, likely through the regulation of immunomodulatory activities, which also influence the gut microbiome. And that larger impact of *S. boulardii* through its influence on host likely contributes much more to the expansion of LMI populations. This dual effect is also seen in the plasma metabolome and circulating cytokines data, in my opinion.

## Arrows I would draw (to speculate further)

So instead of this linear relationship as a speculation,

```
S. boulardii -> microbiome -> host -> tumor
```

I would speculate that the interactions are likely a multifaceted process that I would model as the following (thse arrows took much longer to draw than I am willing to admit):

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

Probably we all would name these arrows similarly regardless of how much we can explain the mechanisms. The first arrow, (1), is perhaps immunomodulation, which affects both tumors and microbiome. The second named arrow, (2), is likely amino acids, which changes another dimension of the environment, which in return affects the same microbiome. The microbiome in this model is primarily (re)shaped by (1), and the effect is sped up by (2).

Careful readers will likely see that there are no more arrows that connect gut microbes and the tumor in my argument here, which is clearly more of a thought experiment based on my intuition of the role and relevance of the microbiome to cancer than anything else.

## Final words

Of course, as the authors put it, how *transitionary* any of these effects are is the most critical question here. I was much more strict when I was younger, but now I know transitionary is not necessarily bad, and a lack of mechanistic insights shouldn't stop practical applications as long as all the science is transparent, and we are not walked around by specific interest groups that seem to work against society and science itself.

Troels Holger Vaaben and colleagues put together a comprehensive evaluation of the entire system, used bacteria as a marker of a larger shift in it, with a careful consideration of the implications of the yeast that certainly contributes to the all.

What makes this study much more interesting than a gazillion others I see in the human microbiome field is the fact that it has enough dimensionality and care that inspires deeper thought than just a nod or shrug. When authors, after generating ton of data and invest hours and hours and hours of work and then tell their readers "we can't say more, or one way or the other", it is the most geniune of invitations for everyone to think more. It creates space to think more. Makes room for more narratives. Allows additional ideas to join in the fight, rather than going against.

I am very proud that our work could support such elaborate science, and I hope the authors could forgive me for saying so much about their work without consulting with them.
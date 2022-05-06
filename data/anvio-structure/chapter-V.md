---
layout: page
title: Chapter V - Reproducing Kiefl et al, 2022
modified: 2021-10-21
excerpt: "A complete reproducible workflow of the manuscript 'Structure-informed microbial population genetics elucidate selective pressures that shape protein evolution' by Kiefl et al"
comments: true
authors: [evan]
redirect_from:
  - 
---

{% capture images %}{{site.url}}/data/anvio-structure/images{% endcapture %}
{% capture command_style %}background: #D7484822; border: 4px solid #D74848;{% endcapture %}
{% capture analysis_style %}background: #E6DBE4{% endcapture %}

{:.warning}
This document is **UNDER CONSTRUCTION**. It is not in a state where you can yet reproduce our work. We anticipate this workflow will be finalized by late March, and will remove this message when it is complete.

## Quick Navigation

- [Chapter I: The prologue]({{ site.url }}/data/anvio-structure/chapter-I)
- [Chapter II: Configure your system]({{ site.url }}/data/anvio-structure/chapter-II)
- [Chapter III: Build the data]({{ site.url }}/data/anvio-structure/chapter-III)
- [Chapter IV: Analyze the data]({{ site.url }}/data/anvio-structure/chapter-IV)
- [Chapter V: Reproduce every number]({{ site.url }}/data/anvio-structure/chapter-V) ← _you are here_

## Preface

This Chapter is dedicated exclusively to reproducing all of the numbers found in the text of the manuscript. All numbers are categorized by the various sections they appear in, and are presented in the order they appear in the section.

Each number can be reproduced by running the associated code block.

Each code block begins with a comment that specifies the programming language used, _e.g._ `# python`. The three languages used are `python`, `bash`, and `R`.

If the language is either `python` or `bash`, it is assumed that your working directory is the root directory of the project, _i.e._ `/some/path/that/ends/in/kiefl_2021`.

If the language is `R`, it is assumed that the commands are issued from a fully initialized [GRE](FIXME). Once you've built your GRE, you can ensure your GRE is fully initialized by issuing the following commands:

```R
request_scvs <- TRUE
request_regs <- TRUE
withRestarts(source("load_data.R"), terminate=function() message('load_data.R: data already loaded. Nice.'))
```

Finally, it is assumed you have completed all [Steps](FIXME) from Chapter III (or equivalently downloaded the [last checkpoint datapack](FIXME)). Some code blocks may require you to complete specific [Analyses](FIXME) from Chapter IV. Such steps will specify, with a comment, which Analysis is a prerequisite.

Finally finally, some numbers that are repeatedly multiple times throughout the text. A code block is only provided for the first instance in which it appears.


## Main

### Main

<blockquote>
(...) resulting in 390 million reads that were <span style="color:red">94.5</span>% identical to HIMB83 on average (Figure S1) (...)
</blockquote>

```R
# R
# Prerequisite: Analysis X
args <- list()
args$reads <- "../18_PERCENT_ID.txt"
args$pfam <- "../18_PERCENT_ID_PFAM.txt"
args$pangenome <- "../18_PERCENT_ID_PANGENOME.txt"
reads <- read_tsv(args$reads) %>%
pivot_longer(cols=-gene_callers_id) %>%
group_by(gene_callers_id) %>%
summarise(reads = mean(value))
pfam <- read_tsv(args$pfam) %>%
rename(pfam = percent_id)
pan <- read_tsv(args$pangenome) %>%
rename(pan = percent_id)
df_wide <- reads %>%
left_join(pfam) %>%
left_join(pan)
df <- df_wide %>%
pivot_longer(cols=-gene_callers_id)
options(pillar.sigfig = 5)
print(df %>% group_by(name) %>% summarise(mean=mean(value, na.rm=TRUE), median=median(value, na.rm=TRUE)))
```


----------------------------

<blockquote>
(...) Of the <span style="color:red">1,470</span> genes in HIMB83 (...)
</blockquote>

```bash
# bash
anvi-export-gene-calls -c CONTIGS.db -o temp --gene-caller prodigal
cat temp | tail +2 | wc -l
``` 


----------------------------

<blockquote>
(...) we restricted our analysis to <span style="color:red">799</span> genes that we determined to form the 1a.3.V core genes (...)
</blockquote>

```bash
# bash
wc -l goi
``` 


----------------------------

<blockquote>
(...) and <span style="color:red">74</span> metagenomes in which the average coverage of HIMB83 exceeded 50X (...)
</blockquote>

```bash
# bash
wc -l soi
``` 

### Polymorphism rates reveal intense purification of nonsynonymous mutants


<blockquote>
(...) Within the 1a.3.V core genes, we found a total of <span style="color:red">9,537,022</span> SCVs, or 128,879 per metagenome on average (...)
</blockquote>

```python
# python
import pandas as pd
scv = pd.read_csv("11_SCVs.txt", sep='\t')
scv = scv[scv['departure_from_consensus'] > 0] # remove quince and kiefl entries
scv.shape[0]
``` 


----------------------------

<blockquote>
(...) Within the 1a.3.V core genes, we found a total of 9,537,022 SCVs, or <span style="color:red">128,879</span> per metagenome on average (...)
</blockquote>

```python
# python
import pandas as pd
scv = pd.read_csv("11_SCVs.txt", sep='\t')
scv = scv[scv['departure_from_consensus'] > 0] # remove quince and kiefl entries
scv.groupby("sample_id")['unique_pos_identifier'].count().mean()
``` 


----------------------------

<blockquote>
(...) SCVs distributed throughout the genome such that <span style="color:red">78</span>% of codons (32% of nucleotides) exhibited minor allele frequencies >10% in at least one metagenome (...)
</blockquote>

```R
# R
tmp <- scvs %>% group_by(unique_pos_identifier) %>% summarize(max_var = max(departure_from_consensus, na.rm=T))
tmp %>% filter(max_var > 0.1) %>% dim() %>% .[[1]] / tmp %>% dim() %>% .[[1]]
``` 


----------------------------

<blockquote>
(...) SCVs distributed throughout the genome such that 78% of codons (<span style="color:red">32</span>% of nucleotides) exhibited minor allele frequencies >10% in at least one metagenome (...)
</blockquote>

**first do this**

```bash
# bash
anvi-gen-variability-profile -c CONTIGS.db -p PROFILE.db --samples-of-interest soi --genes-of-interest goi --engine NT -o 11_SNVs.txt
``` 

**then do this**

```R
# R
goi <- read_tsv("../goi", col_names=F) %>% pull(X1)
snvs <- read_tsv("../11_SNVs.txt")
tmp <- snvs %>% filter(corresponding_gene_call %in% goi) %>% group_by(unique_pos_identifier) %>% summarize(max_df = max(departure_from_consensus))
total_nts <- scvs %>% pull(unique_pos_identifier) %>% unique() %>% length() * 3
nts <- tmp %>% filter(max_df > 0.1) %>% dim() %>% .[[1]]
nts/total_nts
```


----------------------------

<blockquote>
(...) our read recruitment strategy is stringent and yields reads that on average differ from HIMB83 in only <span style="color:red">6</span> nucleotides out of 100 (Table S2) (...)
</blockquote>

```R
# R
args <- list()
args$reads <- "../18_PERCENT_ID.txt"
args$pfam <- "../18_PERCENT_ID_PFAM.txt"
args$pangenome <- "../18_PERCENT_ID_PANGENOME.txt"
reads <- read_tsv(args$reads) %>%
pivot_longer(cols=-gene_callers_id) %>%
group_by(gene_callers_id) %>%
summarise(reads = mean(value))
pfam <- read_tsv(args$pfam) %>%
rename(pfam = percent_id)
pan <- read_tsv(args$pangenome) %>%
rename(pan = percent_id)
df_wide <- reads %>%
left_join(pfam) %>%
left_join(pan)
df <- df_wide %>%
pivot_longer(cols=-gene_callers_id)
options(pillar.sigfig = 5)
tmp <- df %>% group_by(name) %>% summarise(mean=mean(value, na.rm=TRUE), median=median(value, na.rm=TRUE))
round(100-tmp[3,2], 0)
``` 


----------------------------

<blockquote>
(...) Overall, we found that the average pS(site) outweighed pN(site) by <span style="color:red">19:1</span> (Table S3) (...)
</blockquote>

```python
# python
import pandas as pd
df = pd.read_csv('11_SCVs.txt', sep='\t')
print(df.pS_popular_consensus.sum()/df.pN_popular_consensus.sum())
``` 

### Nonsynonymous polymorphism avoids buried sites

<blockquote>
(...) We used two independent methods to predict protein structures for the 799 core genes of 1a.3.V: (1) a template-based homology modeling approach with MODELLER (Webb and Sali 2016), which predicted <span style="color:red">346</span> structures, and (2) a transformer-like deep learning approach with AlphaFold (Jumper et al. 2021), which predicted 754 (...)
</blockquote>

```bash
# bash
wc -l 12_GENES_WITH_GOOD_STRUCTURES_MODELLER
``` 

----------------------------

<blockquote>
(...) We used two independent methods to predict protein structures for the 799 core genes of 1a.3.V: (1) a template-based homology modeling approach with MODELLER (Webb and Sali 2016), which predicted 346 structures, and (2) a transformer-like deep learning approach with AlphaFold (Jumper et al. 2021), which predicted <span style="color:red">754</span> (...)
</blockquote>

```bash
# bash
wc -l 12_GENES_WITH_GOOD_STRUCTURES
``` 

----------------------------

<blockquote>
(...) Our evaluation of the <span style="color:red">339</span> genes for which both methods predicted structures (Supplementary Information) revealed a comparable accuracy between AlphaFold and MODELLER (Figure S4, Table S4) (...)
</blockquote>

```bash
# bash
goi_af = set([int(x.strip()) for x in open('12_GENES_WITH_GOOD_STRUCTURES').readlines()])
goi_mod = set([int(x.strip()) for x in open('12_GENES_WITH_GOOD_STRUCTURES_MODELLER').readlines()])
len(goi_af.intersection(goi_mod))
``` 


### Nonsynonymous polymorphism avoids active sites

<blockquote>
(...) [Figure 2:] The pN$^{(site)}$ distribution (red line) and pS$^{(site)}$ distribution (blue line) were created by weighting the RSA values of <span style="color:red">239,528</span> sites (coming from the 754 genes with predicted structures) by the pN$^{(site)}$ and pS$^{(site)}$ values observed in each of the 74 samples, totaling 17,725,072 pN$^{(site)}$ and pS$^{(site)}$ values (...)
</blockquote>

```R
# R
scvs %>% filter(!is.na(rel_solvent_acc)) %>% .$unique_pos_identifier %>% unique() %>% length()
``` 

----------------------------

<blockquote>
(...) [Figure 2:] The pN$^{(site)}$ distribution (red line) and pS$^{(site)}$ distribution (blue line) were created by weighting the RSA values of 239,528 sites (coming from the 754 genes with predicted structures) by the pN$^{(site)}$ and pS$^{(site)}$ values observed in each of the 74 samples, totaling <span style="color:red">17,725,072</span> pN$^{(site)}$ and pS$^{(site)}$ values (...)
</blockquote>

```R
# R
scvs %>% filter(!is.na(rel_solvent_acc)) %>% dim() %>% .[[1]]
``` 

----------------------------

<blockquote>
(...) [Figure 2:] The pN$^{(site)}$ distribution (red line) and pS$^{(site)}$ distribution (blue line) were created by weighting the DTL values of <span style="color:red">155,478</span> sites (coming from 415 genes that had predicted structures and at least one predicted ligand) by the pN$^{(site)}$ and pS$^{(site)}$ values observed in each of the 74 samples, totaling 11,505,372 pN$^{(site)}$ and pS$^{(site)}$ values (...)
</blockquote>

```R
# R
scvs %>% filter(!is.na(ANY_dist)) %>% .$unique_pos_identifier %>% unique() %>% length()
``` 

----------------------------

<blockquote>
(...) [Figure 2:] The pN$^{(site)}$ distribution (red line) and pS$^{(site)}$ distribution (blue line) were created by weighting the DTL values of 155,478 sites (coming from <span style="color:red">415</span> genes that had predicted structures and at least one predicted ligand) by the pN$^{(site)}$ and pS$^{(site)}$ values observed in each of the 74 samples, totaling 11,505,372 pN$^{(site)}$ and pS$^{(site)}$ values (...)
</blockquote>

```R
# R
scvs %>% filter(!is.na(ANY_dist)) %>% .$gene_callers_id %>% unique() %>% length()
``` 

----------------------------

<blockquote>
(...) [Figure 2:] The pN$^{(site)}$ distribution (red line) and pS$^{(site)}$ distribution (blue line) were created by weighting the DTL values of 155,478 sites (coming from 415 genes that had predicted structures and at least one predicted ligand) by the pN$^{(site)}$ and pS$^{(site)}$ values observed in each of the 74 samples, totaling <span style="color:red">11,505,372</span> pN$^{(site)}$ and pS$^{(site)}$ values (...)
</blockquote>

```R
# R
scvs %>% filter(!is.na(ANY_dist)) %>% dim() %>% .[[1]]
``` 

----------------------------

<blockquote>
(...) A model has been fit to each gene-sample pair that passed filtering criteria (see Supplementary Information), resulting in <span style="color:red">16,285</span> nonsynonymous models and 24,553 synonymous models. (...)
</blockquote>

```R
# R
pn_models %>% dim() %>% .[[1]]
``` 

----------------------------

<blockquote>
(...) A model has been fit to each gene-sample pair that passed filtering criteria (see Supplementary Information), resulting in 16,285 nonsynonymous models and <span style="color:red">24,553</span> synonymous models. (...)
</blockquote>

```R
# R
ps_models %>% dim() %>% .[[1]]
``` 

----------------------------

<blockquote>
(...) The average per-site ns-polymorphism rate throughout the 1a.3.V core genome was <span style="color:red">0.0088</span>, however, we observed a nearly 4-fold reduction in this rate to just 0.0024 at predicted ligand binding sites (DTL = 0) (...)
</blockquote>

```R
# R
scvs %>% pull(pN_popular_consensus) %>% mean(na.rm=T)
``` 

----------------------------

<blockquote>
(...) The average per-site ns-polymorphism rate throughout the 1a.3.V core genome was 0.0088, however, we observed a nearly 4-fold reduction in this rate to just <span style="color:red">0.0024</span> at predicted ligand binding sites (DTL = 0) (...)
</blockquote>

```R
# R
scvs %>% filter(ANY_dist == 0) %>% pull(pN_popular_consensus) %>% mean(na.rm=T)
``` 

### Proteomic trends in purifying selection are explained by RSA and DTL

<blockquote>
(...) yada yada yada (...)
<span style="color:red">94.5</span>
</blockquote>

```bash
# bash
<commandhere>
``` 



























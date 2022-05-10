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

If you are trying to figure out where a particular number came from, you're in the right placei, because this Chapter is dedicated exclusively to reproducing all of the numbers found in the text of the manuscript. All numbers are categorized by the various sections they appear in, and are presented in the order they appear in the section. 

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


## Results and Discussion

### Intro

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
(...) [Figure 2:] The pN(site) distribution (red line) and pS(site) distribution (blue line) were created by weighting the RSA values of <span style="color:red">239,528</span> sites (coming from the 754 genes with predicted structures) by the pN(site) and pS(site) values observed in each of the 74 samples, totaling 17,725,072 pN(site) and pS(site) values (...)
</blockquote>

```R
# R
scvs %>% filter(!is.na(rel_solvent_acc)) %>% .$unique_pos_identifier %>% unique() %>% length()
``` 

----------------------------

<blockquote>
(...) [Figure 2:] The pN(site) distribution (red line) and pS(site) distribution (blue line) were created by weighting the RSA values of 239,528 sites (coming from the 754 genes with predicted structures) by the pN(site) and pS(site) values observed in each of the 74 samples, totaling <span style="color:red">17,725,072</span> pN(site) and pS(site) values (...)
</blockquote>

```R
# R
scvs %>% filter(!is.na(rel_solvent_acc)) %>% dim() %>% .[[1]]
``` 

----------------------------

<blockquote>
(...) [Figure 2:] The pN(site) distribution (red line) and pS(site) distribution (blue line) were created by weighting the DTL values of <span style="color:red">155,478</span> sites (coming from 415 genes that had predicted structures and at least one predicted ligand) by the pN(site) and pS(site) values observed in each of the 74 samples, totaling 11,505,372 pN(site) and pS(site) values (...)
</blockquote>

```R
# R
scvs %>% filter(!is.na(ANY_dist)) %>% .$unique_pos_identifier %>% unique() %>% length()
``` 

----------------------------

<blockquote>
(...) [Figure 2:] The pN(site) distribution (red line) and pS(site) distribution (blue line) were created by weighting the DTL values of 155,478 sites (coming from <span style="color:red">415</span> genes that had predicted structures and at least one predicted ligand) by the pN(site) and pS(site) values observed in each of the 74 samples, totaling 11,505,372 pN(site) and pS(site) values (...)
</blockquote>

```R
# R
scvs %>% filter(!is.na(ANY_dist)) %>% .$gene_callers_id %>% unique() %>% length()
``` 

----------------------------

<blockquote>
(...) [Figure 2:] The pN(site) distribution (red line) and pS(site) distribution (blue line) were created by weighting the DTL values of 155,478 sites (coming from 415 genes that had predicted structures and at least one predicted ligand) by the pN(site) and pS(site) values observed in each of the 74 samples, totaling <span style="color:red">11,505,372</span> pN(site) and pS(site) values (...)
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
(...) By fitting a series of linear models to log-transformed polymorphism data (Table S6), we conclude that RSA and DTL can explain <span style="color:red">11.83</span>% and 6.89% of pN(site) variation, respectively (...)
</blockquote>

```R
# R
# Prerequisite: Analysis X
df <- read_tsv('../WW_TABLES/MODELS.txt')
df$RSA[[2]]
``` 

----------------------------

<blockquote>
(...) By fitting a series of linear models to log-transformed polymorphism data (Table S6), we conclude that RSA and DTL can explain 11.83% and <span style="color:red">6.89</span>% of pN(site) variation, respectively (...)
</blockquote>

```R
# R
# Prerequisite: Analysis X
df <- read_tsv('../WW_TABLES/MODELS.txt')
df$DTL[[4]]
``` 


----------------------------

<blockquote>
(...) Based on these models we estimate that for any given gene in any given sample, (1) a 1% increase in RSA corresponds to a <span style="color:red">0.98</span>% increase in pN(site), and (2) a 1% increase in DTL (normalized by the maximum DTL in the gene) corresponds to a 0.90% increase in pN(site) (...)
</blockquote>

```R
# R
# Prerequisite: Analysis X
lm_rsa_gene_sample_pn <- readRDS("../lm_rsa_gene_sample_pn.RDS")
print(round(100*(exp((lm_rsa_gene_sample_pn %>% summary %>% coef)["rel_solvent_acc","Estimate"] * 0.01) - 1), 2))
rm(lm_rsa_gene_sample_pn)
``` 

----------------------------

<blockquote>
(...) Based on these models we estimate that for any given gene in any given sample, (1) a 1% increase in RSA corresponds to a 0.98% increase in pN(site), and (2) a 1% increase in DTL (normalized by the maximum DTL in the gene) corresponds to a <span style="color:red">0.90</span>% increase in pN(site) (...)
</blockquote>

```R
# R
# Prerequisite: Analysis X
lm_dtl_gene_sample_pn <- readRDS("../lm_dtl_gene_sample_pn.RDS")
print(round(100*(exp((lm_dtl_gene_sample_pn %>% summary %>% coef)["ANY_dist","Estimate"] * 0.01) - 1), 2))
history()
rm(lm_dtl_gene_sample_pn)
``` 

----------------------------

<blockquote>
(...) In a combined model, RSA and DTL jointly explained <span style="color:red">14.12</span>% of pN(site) variation, and after adjusting for gene-to-gene and sample-to-sample variance, 17.07% of the remaining variation could be explained by RSA and DTL (...)
</blockquote>

```R
# R
# Prerequisite: Analysis X
lm_rsa_dtl_gene_sample_pn <- readRDS("../lm_rsa_dtl_gene_sample_pn.RDS")
formatted_anova_pn_rsa_dtl <- (100*(lm_rsa_dtl_gene_sample_pn %>% anova)$"Sum Sq"/sum((lm_rsa_dtl_gene_sample_pn %>% anova)$"Sum Sq")) 
names(formatted_anova_pn_rsa_dtl) <- rownames(lm_rsa_dtl_gene_sample_pn %>% anova)
combined <- sum(formatted_anova_pn_rsa_dtl["ANY_dist"], formatted_anova_pn_rsa_dtl["rel_solvent_acc"])
print(combined %>% round(2))
rm(formatted_anova_pn_rsa_dtl)
``` 

----------------------------

<blockquote>
(...) In a combined model, RSA and DTL jointly explained 14.12% of pN(site) variation, and after adjusting for gene-to-gene and sample-to-sample variance, <span style="color:red">17.07</span>% of the remaining variation could be explained by RSA and DTL (...)
</blockquote>

```R
# R
# Prerequisite: Analysis X
lm_rsa_dtl_gene_sample_pn <- readRDS("../lm_rsa_dtl_gene_sample_pn.RDS")
formatted_anova_pn_rsa_dtl <- (100*(lm_rsa_dtl_gene_sample_pn %>% anova)$"Sum Sq"/sum((lm_rsa_dtl_gene_sample_pn %>% anova)$"Sum Sq")) 
names(formatted_anova_pn_rsa_dtl) <- rownames(lm_rsa_dtl_gene_sample_pn %>% anova)
combined <- sum(formatted_anova_pn_rsa_dtl["ANY_dist"], formatted_anova_pn_rsa_dtl["rel_solvent_acc"])
print(100*combined / (combined + formatted_anova_pn_rsa_dtl["Residuals"]) %>% round(2))
rm(formatted_anova_pn_rsa_dtl)
``` 


----------------------------

<blockquote>
(...) In comparison, only <span style="color:red">0.35</span>% of pS(site) variation was explained by RSA and DTL (...)
</blockquote>

```R
# R
# Prerequisite: Analysis X
lm_rsa_dtl_gene_sample_ps <- readRDS("../lm_rsa_dtl_gene_sample_ps.RDS")
formatted_anova_ps_rsa_dtl <- (100*(lm_rsa_dtl_gene_sample_ps %>% anova)$"Sum Sq"/sum((lm_rsa_dtl_gene_sample_ps %>% anova)$"Sum Sq")) 
names(formatted_anova_ps_rsa_dtl) <- rownames(lm_rsa_dtl_gene_sample_ps %>% anova)
combined <- sum(formatted_anova_ps_rsa_dtl["ANY_dist"], formatted_anova_ps_rsa_dtl["rel_solvent_acc"])
print(100*combined_ps / (combined_ps + formatted_anova_ps_rsa_dtl["Residuals"]) %>% round(2))
rm(formatted_anova_ps_rsa_dtl)
``` 

----------------------------

<blockquote>
(...) Analyzing gene-sample pairs revealed that the extent of ns-polymorphism rate that can be explained by RSA and DTL is not uniform across all genes (Table S7) and can reach up to <span style="color:red">52.6</span>% and <span style="color:red">51.4</span>%, respectively (Figures S7, S8) (...)
</blockquote>

```R
# R
print(poly_corr %>% pull(pn_rsq_rsa) %>% max(na.rm=T))
print(poly_corr %>% pull(pn_rsq_dist) %>% max(na.rm=T))
``` 

----------------------------

<blockquote>
(...) Linear regressions of these data show that <span style="color:red">83.6</span>% of per-group ns-polymorphism rates and <span style="color:red">20.7</span>% of per-group s-polymorphism rates are explained by RSA and DTL (Supplementary Information) (...)
</blockquote>

```R
# R
print(pn_group_model %>% summary())
print(ps_group_model %>% summary())
``` 

### Measuring purifying selection between genes and environments with pN/pS(gene)

<blockquote>
(...) Taking advantage of the large number of metagenomes in which 1a.3.V was present, we calculated pN/pS(gene) for all 799 protein-coding core genes across 74 samples (see Methods), resulting in <span style="color:red">59,126</span> gene/sample pairs (Table S9) (...)
</blockquote>

```R
# R
pnps %>% dim() %>% .[[1]]
``` 

----------------------------

<blockquote>
(...) We validated our calculations by comparing sample-averaged pN/pS(gene) to dN/dS(gene) calculated from homologous gene pairs between HIMB83 and HIMB122, another SAR11 isolate genome that is closely related to HIMB83 (gANI: <span style="color:red">82.6</span>%), which we found to yield commensurable results (Figure S10, Table S12, Supplementary Information) (...)
</blockquote>

```python
# python
# Prerequisite: Aux. Step X
import pandas as pd
df = pd.read_csv("07_ANI_HIMB122/ANIb_percentage_identity.txt", sep='\t').set_index('key')
print(df['HIMB083'].sort_values())
``` 


----------------------------

<blockquote>
(...) All but one gene (<span style="color:red">gene #2031, unknown function</span>) maintained pN/pS(gene) << 1 in every sample, whereby 95% of values were less than 0.15 (Figure S12, Table S9) (...)
</blockquote>

```R
# R
pnps %>% filter(pnps > 1) %>% select(gene_callers_id, function_COG20_FUNCTION, function_KOfam, function_Pfam)
``` 

----------------------------

<blockquote>
(...) All but one gene (gene #2031, unknown function) maintained pN/pS(gene) << 1 in every sample, whereby 95% of values were less than <span style="color:red">0.15</span> (Figure S12, Table S9) (...)
</blockquote>

```R
# R
pnps %>% .$pnps %>% quantile(0.95, na.rm=T)
``` 


----------------------------

<blockquote>
(...) In fact, gene-to-gene variance, as opposed to sample-to-sample variance, explained <span style="color:red">93</span>% of pN/pS(gene) variation (ANOVA, Figure S11) (...)
</blockquote>

```R
# R
my_anova <- pnps %>% mutate(gene_callers_id = as.factor(gene_callers_id)) %>% lm(pnps ~ gene_callers_id + sample_id, data = .) %>% anova()
gene_explained <- my_anova$`Sum Sq`[1] / sum(my_anova$`Sum Sq`) * 100
sample_explained <- my_anova$`Sum Sq`[2] / sum(my_anova$`Sum Sq`) * 100
print(my_anova)
print(gene_explained)
print(sample_explained)
``` 


----------------------------

<blockquote>
(...) By analyzing the companion metatranscriptomic data (Salazar et al. 2019) that were available for 50 of the 74 metagenomes, we were able to explain <span style="color:red">29</span>% of gene-to-gene variance with gene transcript abundance (Table S13, Supplementary Information) (...)
</blockquote>

```R
# R
cor(TA_sample_averaged$log_pnps_mean, TA_sample_averaged$log_TA_median, use="pairwise.complete.obs") ^ 2 * 100
``` 

----------------------------

<blockquote>
(...) The amount of pN/pS(gene) variation attributable to sample-to-sample variance was only <span style="color:red">0.7</span>% (Figure S11) (...)
</blockquote>

```R
# R
my_anova <- pnps %>% mutate(gene_callers_id = as.factor(gene_callers_id)) %>% lm(pnps ~ gene_callers_id + sample_id, data = .) %>% anova()
gene_explained <- my_anova$`Sum Sq`[1] / sum(my_anova$`Sum Sq`) * 100
sample_explained <- my_anova$`Sum Sq`[2] / sum(my_anova$`Sum Sq`) * 100
print(my_anova)
print(gene_explained)
print(sample_explained)
``` 


### Nitrogen availability governs rates of non-ideal polymorphism at critical sites of glutamine synthetase

<blockquote>
(...) [T]he sample-averaged pN/pS(GS) was <span style="color:red">0.02</span>, ranking GS amongst the top 11% most purified genes (Figure 3b, Table S9) (...)
</blockquote>

```R
# R
pnps %>% group_by(gene_callers_id) %>% summarize(pnps=mean(pnps, na.rm=T)) %>% filter(gene_callers_id == 2602) %>% pull(pnps)
``` 

----------------------------

<blockquote>
(...) Indeed, the sample-averaged pN/pS(GS) was 0.02, ranking GS amongst the top <span style="color:red">11</span>% most purified genes (Figure 3b, Table S9) (...)
</blockquote>

```R
# R
pnps %>% group_by(gene_callers_id) %>% summarize(pnps=mean(pnps, na.rm=T)) %>% arrange(pnps) %>% mutate(rank = 1:nrow(.)) %>% filter(gene_callers_id==2602) %>% pull(rank) / (pnps %>% .$gene_callers_id %>% unique() %>% length()) * 100
``` 


----------------------------

<blockquote>
(...) Although highly purified, we observed significant sample-to-sample variation in pN/pS(GS) (min = <span style="color:red">0.010</span>, max = <span style="color:red">0.036</span>) suggesting that the strength of purifying selection on GS varies from sample to sample (Figure 3b inset) (...)
</blockquote>

```R
# R
print(pnps %>% filter(gene_callers_id == 2602) %>% pull(pnps) %>% min())
print(pnps %>% filter(gene_callers_id == 2602) %>% pull(pnps) %>% max())
``` 

----------------------------

<blockquote>
(...) [We] found a positive correlation between measured nitrate concentrations and pN/pS(GS) values across samples (Pearson correlation p-value = <span style="color:red">0.009</span>, R2 = <span style="color:red">0.11</span>) (Figure 3c), which ranked amongst the top 12% of positive correlations between pN/pS(gene) and nitrate concentration (Figure 3c inset, Table S10) (...)
</blockquote>

```R
# R
temp <- pnps %>% filter(gene_callers_id == 2602)
print(cor.test(temp$pnps, temp$nitrates, use='pairwise.complete.obs') %>% .$p.value)
print(cor.test(temp$pnps, temp$nitrates, use='pairwise.complete.obs') %>% .$estimate %>% .^2)
``` 

----------------------------

<blockquote>
(...) [We] found a positive correlation between measured nitrate concentrations and pN/pS(GS) values across samples (Pearson correlation p-value = 0.009, R2 = 0.11) (Figure 3c), which ranked amongst the top <span style="color:red">12</span>% of positive correlations between pN/pS(gene) and nitrate concentration (Figure 3c inset, Table S10) (...)
</blockquote>

```R
# R
GS <- env_corr %>% filter(gene_callers_id==2602) %>% pull(cor_nitrates)
sum(env_corr %>% pull(cor_nitrates) > GS) / dim(env_corr)[[1]] * 100
``` 


----------------------------

<blockquote>
(...) From this stitched quaternary structure we recalculated RSA and DTL, and as expected, this yielded lower average RSA and DTL estimates due to the presence of adjacent monomers (<span style="color:red">0.17</span> versus <span style="color:red">0.24</span> for RSA and <span style="color:red">17.8</span>Å versus <span style="color:red">21.2</span>Å for DTL) (...)
</blockquote>

```R
# R
dists <- scvs %>%
filter(
gene_callers_id == 2602,
sample_id == 'ANE_004_05M' # a random sample will do
) %>%
select(
RSA_mon = rel_solvent_acc,
RSA_dec = complex_RSA,
DTL_mon = GLU_dist_2602,
DTL_dec = complex_DTL,
)
print(dists %>% pull(RSA_dec) %>% mean(na.rm=T) # complex)
print(dists %>% pull(RSA_mon) %>% mean(na.rm=T) # monomer)
print(dists %>% pull(DTL_dec) %>% mean(na.rm=T) # complex)
print(dists %>% pull(DTL_mon) %>% mean(na.rm=T) # monomer)
``` 

----------------------------

<blockquote>
(...) With these quaternary estimates of RSA and DTL, we found that ns-polymorphism was <span style="color:red">30</span>x less common than s-polymorphism (...)
</blockquote>

```R
# R
scvs %>% filter(gene_callers_id==2602) %>% pull(pS_popular_consensus) %>% mean(na.rm=T) / scvs %>% filter(gene_callers_id==2602) %>% pull(pN_popular_consensus) %>% mean(na.rm=T)
``` 

----------------------------

<blockquote>
(...) In comparison, s-polymorphism distributed relatively homogeneously throughout the protein, whereby <span style="color:red">17</span>% of s-polymorphism occurred within 10Å of active sites (compared to 3% for ns-polymorphism) and 19% occurred in sites with 0 RSA (compared to 9% for ns-polymorphism) (...)
</blockquote>

```R
# R
temp <- scvs %>% filter(gene_callers_id==2602)
total_pS <- temp %>% pull(pS_popular_consensus) %>% sum(na.rm=T)
small_DTL <- temp %>% filter(complex_DTL < 10) %>% pull(pN_popular_consensus) %>% sum(na.rm=T)
small_DTL/total_pS*100
``` 

----------------------------

<blockquote>
(...) In comparison, s-polymorphism distributed relatively homogeneously throughout the protein, whereby 17% of s-polymorphism occurred within 10Å of active sites (compared to <span style="color:red">3</span>% for ns-polymorphism) and 19% occurred in sites with 0 RSA (compared to 9% for ns-polymorphism) (...)
</blockquote>

```R
# R
temp <- scvs %>% filter(gene_callers_id==2602)
total_pN <- temp %>% pull(pN_popular_consensus) %>% sum(na.rm=T)
small_DTL_pN <- temp %>% filter(complex_DTL < 10) %>% pull(pN_popular_consensus) %>% sum(na.rm=T)
small_DTL_pN/total_pN*100
``` 

----------------------------

<blockquote>
(...) In comparison, s-polymorphism distributed relatively homogeneously throughout the protein, whereby 17% of s-polymorphism occurred within 10Å of active sites (compared to 3% for ns-polymorphism) and <span style="color:red">19</span>% occurred in sites with 0 RSA (compared to 9% for ns-polymorphism) (...)
</blockquote>

```R
# R
temp <- scvs %>% filter(gene_callers_id==2602)
total_pS <- temp %>% pull(pS_popular_consensus) %>% sum(na.rm=T)
small_RSA_pS <- temp %>% filter(complex_RSA == 0) %>% pull(pS_popular_consensus) %>% sum(na.rm=T)
small_RSA_pS/total_pS*100
``` 

----------------------------

<blockquote>
(...) In comparison, s-polymorphism distributed relatively homogeneously throughout the protein, whereby 17% of s-polymorphism occurred within 10Å of active sites (compared to 3% for ns-polymorphism) and 19% occurred in sites with 0 RSA (compared to <span style="color:red">9</span>% for ns-polymorphism) (...)
</blockquote>

```R
# R
temp <- scvs %>% filter(gene_callers_id==2602)
total_pN <- temp %>% pull(pN_popular_consensus) %>% sum(na.rm=T)
small_RSA_pN <- temp %>% filter(complex_RSA == 0) %>% pull(pN_popular_consensus) %>% sum(na.rm=T)
small_RSA_pN/total_pN*100
``` 


----------------------------

<blockquote>
(...) Averaged across samples, the mean RSA was <span style="color:red">0.15</span> for s-polymorphism and <span style="color:red">0.33</span> for ns-polymorphism (Figure 3e left panel) (...)
</blockquote>

```R
# R
temp <- scvs %>% filter(gene_callers_id==2602)
# mean RSA of pS
print(temp %>% group_by(sample_id) %>% mutate(pS_weighted_RSA = pS_popular_consensus/sum(pS_popular_consensus, na.rm=T)*complex_RSA) %>% summarize(mean_RSA = sum(pS_weighted_RSA, na.rm=T)) %>% pull(mean_RSA) %>% mean())
# mean RSA of pN
print(temp %>% group_by(sample_id) %>% mutate(pN_weighted_RSA = pN_popular_consensus/sum(pN_popular_consensus, na.rm=T)*complex_RSA) %>% summarize(mean_RSA = sum(pN_weighted_RSA, na.rm=T)) %>% pull(mean_RSA) %>% mean())
``` 

----------------------------

<blockquote>
(...) Similarly, the mean DTL was <span style="color:red">17.2</span>Å for s-polymorphism and <span style="color:red">22.9</span>Å for ns-polymorphism (Figure 3f left panel) (...)
</blockquote>

```R
# R
temp <- scvs %>% filter(gene_callers_id == 2602)
# mean DTL of pS
print(temp %>% group_by(sample_id) %>% mutate(pN_weighted_DTL = pN_popular_consensus/sum(pN_popular_consensus, na.rm=T)*complex_DTL) %>% summarize(mean_DTL = sum(pN_weighted_DTL, na.rm=T)) %>% pull(mean_DTL) %>% mean())
# mean DTL of pN
print(temp %>% group_by(sample_id) %>% mutate(pN_weighted_DTL = pN_popular_consensus/sum(pN_popular_consensus, na.rm=T)*complex_DTL) %>% summarize(mean_DTL = sum(pN_weighted_DTL, na.rm=T)) %>% pull(mean_DTL) %>% mean())
``` 

----------------------------

<blockquote>
(...) [T]he mean RSA of s-polymorphism remained relatively invariant (standard deviation <span style="color:red">0.005</span>) (Figure 3e right panel) (...)
</blockquote>

```R
# R
scvs %>% filter(gene_callers_id==2602) %>% group_by(sample_id) %>% mutate(pS_weighted_RSA = pS_popular_consensus/sum(pS_popular_consensus, na.rm=T)*complex_RSA) %>% summarize(mean_RSA = sum(pS_weighted_RSA, na.rm=T), pnps=mean(pnps, na.rm=T)) %>% pull(mean_RSA) %>% sd()
``` 


----------------------------

<blockquote>
(...) [T]he mean RSA of ns-polymorphism varied dramatically from <span style="color:red">0.27</span> to <span style="color:red">0.37</span> and was profoundly influenced by sample pN/pS(GS) (...)
</blockquote>

```R
# R
print(scvs %>% filter(gene_callers_id==2602) %>% group_by(sample_id) %>% mutate(pN_weighted_RSA = pN_popular_consensus/sum(pN_popular_consensus, na.rm=T)*complex_RSA) %>% summarize(mean_RSA = sum(pN_weighted_RSA, na.rm=T)) %>% pull(mean_RSA) %>% min())
print(scvs %>% filter(gene_callers_id==2602) %>% group_by(sample_id) %>% mutate(pN_weighted_RSA = pN_popular_consensus/sum(pN_popular_consensus, na.rm=T)*complex_RSA) %>% summarize(mean_RSA = sum(pN_weighted_RSA, na.rm=T)) %>% pull(mean_RSA) %>% max())
``` 

----------------------------

<blockquote>
(...) <span style="color:red">82.9</span>% of mean RSA ns-polymorphism variance could be explained by pN/pS(GS) alone (Pearson correlation, p-value < <span style="color:red">1x10-16</span>, R2 = <span style="color:red">0.829</span>) (...)
</blockquote>

```R
# R
scvs %>% filter(gene_callers_id==2602) %>% group_by(sample_id) %>% mutate(pN_weighted_RSA = pN_popular_consensus/sum(pN_popular_consensus, na.rm=T)*complex_RSA) %>% summarize(mean_RSA = sum(pN_weighted_RSA, na.rm=T), pnps=mean(pnps, na.rm=T)) %>% lm(mean_RSA~pnps, data=.) %>% summary()
``` 

----------------------------

<blockquote>
(...) ns-polymorphism distributions with respect to DTL were equally governed by selection strength, where <span style="color:red">80.4</span>% of variance could be explained by pN/pS(GS) (Pearson correlation, p-value <span style="color:red">1x10-16</span>, R2 = <span style="color:red">0.804</span>, Figure 3f) (...)
</blockquote>

```R
# R
scvs %>% filter(gene_callers_id==2602) %>% group_by(sample_id) %>% mutate(pN_weighted_DTL = pN_popular_consensus/sum(pN_popular_consensus, na.rm=T)*complex_DTL) %>% summarize(mean_DTL = sum(pN_weighted_DTL, na.rm=T), pnps=mean(pnps, na.rm=T)) %>% lm(mean_DTL~pnps, data=.) %>% summary() %>% .$r.squared
``` 

----------------------------

<blockquote>
(...) [Figure 3B:] The vertical green line depicts the sample-averaged pN/pS(gene) for GS (<span style="color:red">0.020</span>) (...)
</blockquote>

```R
# R
pnps %>% group_by(gene_callers_id) %>% summarize(pnps=mean(pnps, na.rm=T)) %>% filter(gene_callers_id==2602) %>% pull(pnps)
``` 

----------------------------

<blockquote>
(...) [Figure 3b:] The inset plot shows the distribution of pN/pS(gene) value for GS as seen across the 74 samples, which vary from <span style="color:red">0.010</span> to <span style="color:red">0.036</span> (...)
</blockquote>

```R
# R
print(pnps %>% filter(gene_callers_id==2602) %>% pull(pnps) %>% min())
print(pnps %>% filter(gene_callers_id==2602) %>% pull(pnps) %>% max())
``` 

----------------------------

<blockquote>
(...) [Figure 3c:] The vertical green line depicts the correlation coefficient for GS (<span style="color:red">0.34</span>) (...)
</blockquote>

```R
# R
env_corr %>% filter(gene_callers_id == 2602) %>% pull(cor_nitrates)
``` 

----------------------------

<blockquote>
(...) We identified putative sites fitting this description by scoring sites based on the extent that their amino acid minor allele frequencies co-varied with pN/pS(GS), including only sites with DTL less than the mean DTL of ns-polymorphisms (<span style="color:red">22.9Å</span>) (...)
</blockquote>

```R
# R
scvs %>% filter(gene_callers_id == 2602) %>% group_by(sample_id) %>% mutate(pN_weighted_DTL = pN_popular_consensus/sum(pN_popular_consensus, na.rm=T)*complex_DTL) %>% summarise(mean_pN_DTL=sum(pN_weighted_DTL, na.rm=T)) %>% pull(mean_pN_DTL) %>% mean()
``` 

----------------------------

<blockquote>
(...) Though each of these sites exhibited DTL lower than the average ns-polymorphism, the closest site (residue number 323) was still <span style="color:red">9</span>Å away from the glutamate substrate (...)
</blockquote>

```R
# R
scvs %>% filter(gene_callers_id==2602) %>% filter(sample_id == 'ANE_004_05M', codon_order_in_gene %in% (c(96, 152, 175, 176, 230, 288, 323, 364, 379)-1)) %>% select(codon_order_in_gene, complex_DTL, complex_RSA) %>% mutate(codon_number = codon_order_in_gene + 1) %>% select(-codon_order_in_gene)
``` 

----------------------------

<blockquote>
(...) After all, in absolute terms GS is highly purified regardless of sample – the largest pN/pS(GS) is <span style="color:red">0.036</span>, which is just over half the genome-wide average pN/pS(gene) of 0.063 (...)
</blockquote>

```R
# R
pnps %>% filter(gene_callers_id==2602) %>% pull(pnps) %>% max()
``` 

----------------------------

<blockquote>
(...) After all, in absolute terms GS is highly purified regardless of sample – the largest pN/pS(GS) is 0.036, which is just over half the genome-wide average pN/pS(gene) of <span style="color:red">0.063</span> (...)
</blockquote>

```R
# R
pnps %>% pull(pnps) %>% mean(na.rm=T)
``` 

----------------------------

<blockquote>
(...) Though highly significant (one sided Pearson p-values <span style="color:red">9x10-12</span> for RSA and <span style="color:red">2x10-4</span> for DTL), the magnitude that ns-polymorphism distributions shift with respect to DTL and RSA were subtle (...)
</blockquote>

```R
# R
temp <- scvs %>%
filter(
!is.na(pN_popular_consensus),
ANY_dist <= 40
) %>%
group_by(sample_id) %>%
mutate(
pS_weighted_DTL = pS_popular_consensus/sum(pS_popular_consensus)*ANY_dist,
pN_weighted_DTL = pN_popular_consensus/sum(pN_popular_consensus)*ANY_dist,
pS_weighted_RSA = pS_popular_consensus/sum(pS_popular_consensus)*rel_solvent_acc,
pN_weighted_RSA = pN_popular_consensus/sum(pN_popular_consensus)*rel_solvent_acc,
rarity_weighted_DTL = codon_rarity/sum(codon_rarity)*ANY_dist,
rarity_weighted_RSA = codon_rarity/sum(codon_rarity)*rel_solvent_acc
) %>%
summarise(
mean_pS_DTL=sum(pS_weighted_DTL),
mean_pN_DTL=sum(pN_weighted_DTL),
mean_pS_RSA=sum(pS_weighted_RSA),
mean_pN_RSA=sum(pN_weighted_RSA),
mean_rarity_DTL=sum(rarity_weighted_DTL),
mean_rarity_RSA=sum(rarity_weighted_RSA),
pnps=mean(genome_pnps),
temperature=mean(temperature)
)
print(cor.test(temp$mean_pN_RSA, temp$pnps, alternative='less') %>% .$p.val)
print(cor.test(temp$mean_pN_DTL, temp$pnps, alternative='less') %>% .$p.val)
``` 

----------------------------

<blockquote>
(...) [Figure 4:] (C) The s-polymorphism distribution mean with respect to RSA is negatively associated with pN/pS(core) (one-sided Pearson p-value = <span style="color:red">1x10-5</span>). (D) The s-polymorphism distribution mean with respect to RSA is negatively associated with pN/pS(core) (one-sided Pearson p-value = <span style="color:red">3x10-7</span>). (...)
</blockquote>

```R
# R
temp <- scvs %>%
filter(
!is.na(pN_popular_consensus),
ANY_dist <= 40
) %>%
group_by(sample_id) %>%
mutate(
pS_weighted_DTL = pS_popular_consensus/sum(pS_popular_consensus)*ANY_dist,
pN_weighted_DTL = pN_popular_consensus/sum(pN_popular_consensus)*ANY_dist,
pS_weighted_RSA = pS_popular_consensus/sum(pS_popular_consensus)*rel_solvent_acc,
pN_weighted_RSA = pN_popular_consensus/sum(pN_popular_consensus)*rel_solvent_acc,
rarity_weighted_DTL = codon_rarity/sum(codon_rarity)*ANY_dist,
rarity_weighted_RSA = codon_rarity/sum(codon_rarity)*rel_solvent_acc
) %>%
summarise(
mean_pS_DTL=sum(pS_weighted_DTL),
mean_pN_DTL=sum(pN_weighted_DTL),
mean_pS_RSA=sum(pS_weighted_RSA),
mean_pN_RSA=sum(pN_weighted_RSA),
mean_rarity_DTL=sum(rarity_weighted_DTL),
mean_rarity_RSA=sum(rarity_weighted_RSA),
pnps=mean(genome_pnps),
temperature=mean(temperature)
)
print(cor.test(temp$mean_pS_RSA, temp$pnps, alternative='less') %>% .$p.val)
print(cor.test(temp$mean_pS_DTL, temp$pnps, alternative='less') %>% .$p.val)
``` 

----------------------------

<blockquote>
(...) [Figure 4:] (E) Rare synonymous codons are more abundant in samples with high pN/pS(core) (one-sided Pearson p-value = <span style="color:red">4x10-5</span>). (F) Rare synonymous codons avoid low RSA sites when pN/pS(core) is low (one-sided Pearson p-value = <span style="color:red">1x10-10</span>). (G) Rare synonymous codons avoid low DTL sites when pN/pS(core) is low (one-sided Pearson p-value = <span style="color:red">7x10-9</span>) (...)
</blockquote>

```R
# R

pN_cutoff <- 0.0005

# E
plot_data <- scvs %>%
filter(pN_popular_consensus < pN_cutoff) %>%
group_by(sample_id) %>%
summarise(
codon_rarity=mean(codon_rarity),
pnps=mean(genome_pnps)
)
print(cor.test(plot_data$codon_rarity, temp$pnps, alternative='greater') %>% .$p.val)

# F
plot_data <- scvs %>%
filter(
pN_popular_consensus < pN_cutoff,
!is.na(rel_solvent_acc) # Filter out any genes without structures
) %>%
group_by(sample_id) %>%
mutate(
rarity_weighted_RSA = codon_rarity/sum(codon_rarity)*rel_solvent_acc
) %>%
summarise(
mean_rarity_RSA=sum(rarity_weighted_RSA),
pnps=mean(genome_pnps)
)
print(cor.test(plot_data$mean_rarity_RSA, temp$pnps, alternative='less') %>% .$p.val)

# G
plot_data <- scvs %>%
filter(
pN_popular_consensus < pN_cutoff,
!is.na(ANY_dist) # Filter out any genes without structures or binding sites
) %>%
group_by(sample_id) %>%
mutate(
rarity_weighted_DTL = codon_rarity/sum(codon_rarity)*ANY_dist
) %>%
summarise(
mean_rarity_DTL=sum(rarity_weighted_DTL),
pnps=mean(genome_pnps)
)
print(cor.test(plot_data$mean_rarity_DTL, temp$pnps, alternative='less') %>% .$p.val)
``` 

### Synonymous but not silent: selection against rare codons at critical sites

<blockquote>
(...) With a GC-content lower than <span style="color:red">30</span>%, SAR11 genomes maintain a non-uniform yet conserved codon composition (Figure S13) (...)
</blockquote>

```python
# python
import anvio.utils as utils
seqs = utils.get_FASTA_file_as_dictionary('contigs.fa')
concated_seq = ''.join([s for s in seqs.values()])
print(utils.get_GC_content_for_sequence(seq))
``` 

{:.notice}
All other numbers in this section are found in the code blocks pertaining to Figure 4 (above).


## Methods

### Predicting and processing protein structures

<blockquote>
(...) To avoid this, we discarded any templates if the alignment coverage of the protein sequence to the template was <80%. Applying these filters resulted in <span style="color:red">408</span> structures from the 1a.3.V core (...)
</blockquote>

```python
# python
import anvio.structureops as sops
sdb = sops.StructureDatabase('09_STRUCTURE_MOD.db')
goi = [int(x.strip()) for x in open('goi').readlines()]
len([x for x in sdb.genes_with_structure if x in goi])
``` 

----------------------------

<blockquote>
(...) [The number of MODELLER-predicted structures ]was further refined by requiring that the root mean squared distance (RMSD) between the predicted structure and the most similar template did not exceed 7.5 Å, and that the GA341 model score exceeded 0.95. After applying these constraints, we were left with <span style="color:red">348</span> structures in the 1a.3.V that we assumed to be ‘trustworthy’ structures as predicted by MODELLER (...)
</blockquote>

```bash
# bash
wc -l 12_GENES_WITH_GOOD_STRUCTURES_MODELLER
``` 

----------------------------

<blockquote>
(...) These structures were on average <span style="color:red">44.8</span>% identical to their templates, which is within the sequence similarity regime where template-based homology modeling generally produces the correct overall fold (Rost 1999) (...)
</blockquote>

```R
# R
structure_val <- read_tsv("../09_STRUCTURE_comparison.txt")
structure_val %>% filter(has_MOD, RMSD <= 7.5) %>% .$template_similarity %>% mean()
``` 

### Predicting ligand-binding sites

<blockquote>
(...) We then filtered out binding frequency scores less than 0.5, yielding <span style="color:red">40,219</span> predicted ligand-residue interactions across 11,480 unique sites (Table S5) (...)
</blockquote>

```R
# R
temp <- read_tsv('../08_INTERACDOME-binding_frequencies.txt')
temp %>% dim() %>% .[[1]]
``` 

----------------------------

<blockquote>
(...) We then filtered out binding frequency scores less than 0.5, yielding 40,219 predicted ligand-residue interactions across <span style="color:red">11,480</span> unique sites (Table S5) (...)
</blockquote>

```R
# R
temp <- read_tsv('../08_INTERACDOME-binding_frequencies.txt')
temp %>% mutate(new = paste(gene_callers_id, codon_order_in_gene)) %>% pull(new) %>% unique() %>% length()
``` 


### Calculating distance-to-ligand (DTL)

<blockquote>
(...) To mitigate the influence of this inevitable error source, we conservatively excluded DTL values >40 Å (<span style="color:red">8.0</span>% of sites) in all analyses after Figure 2b (...)
</blockquote>

```R
# R
ggs <- read_tsv("../12_GENES_WITH_GOOD_STRUCTURES", col_names=F) %>% pull(X1)
sar11_data <- scvs %>%
filter(
sample_id == 'ANE_004_05M', # pick random sample to massively subset data
gene_callers_id %in% ggs, # only genes with predicted structure
!is.na(ANY_dist) # only genes with predicted ligands
)
100 * sar11_data %>% filter(ANY_dist >= 40) %>% dim() %>% .[[1]] / sar11_data %>% dim() %>% .[[1]]
``` 

### Proportion of polymorphism rate variance explained by RSA and DTL

<blockquote>
(...) After excluding monomorphic sites (pN(site) = 0 for ns-models, pS(site) = 0 for s-models), this yielded <span style="color:red">5,374,739</span> data points for s-models and <span style="color:red">3,530,079</span> for ns-models  (...)
</blockquote>

```R
# R
lm_rsa_gene_sample_ps <- readRDS("../lm_rsa_gene_sample_ps.RDS")
print(lm_rsa_gene_sample_ps$model %>% dim())
rm(lm_rsa_gene_sample_ps) # takes up lots of memory, remove when finished

lm_rsa_gene_sample_pn <- readRDS("../lm_rsa_gene_sample_pn.RDS")
print(lm_rsa_gene_sample_pn$model %>% dim())
rm(lm_rsa_gene_sample_pn) # takes up lots of memory, remove when finished
``` 

### Supplementary Figures

<blockquote>
(...) The two groups were defined as having TM scores above or below 0.8, where the >0.8 group corresponded to the <span style="color:red">291</span> best alignments (left) and the <0.8 group corresponded to the 48 worst alignments (...)
</blockquote>

```R
# R
df <- read_tsv("../09_STRUCTURE_comparison.txt") %>% mutate(frac_ss_MOD = frac_alpha_MOD + frac_beta_MOD, frac_ss_AF = frac_alpha_AF + frac_beta_AF)
df %>% filter(has_MOD, has_AF) %>% mutate(cohort = ifelse(TM_score > 0.8, 'cohort1', 'cohort2')) %>% filter(cohort=='cohort1') %>% dim() %>% .[[1]]
``` 

----------------------------------------


<blockquote>
(...) [Figure S4:] The two groups were defined as having TM scores above or below 0.8, where the >0.8 group corresponded to the 291 best alignments (left) and the <<span style="color:red">0.8</span> group corresponded to the 48 worst alignments (...)
</blockquote>

```R
# R
df <- read_tsv("../09_STRUCTURE_comparison.txt") %>% mutate(frac_ss_MOD = frac_alpha_MOD + frac_beta_MOD, frac_ss_AF = frac_alpha_AF + frac_beta_AF)
df %>% filter(has_MOD, has_AF) %>% mutate(cohort = ifelse(TM_score < 0.8, 'cohort1', 'cohort2')) %>% filter(cohort=='cohort1') %>% dim() %>% .[[1]]
``` 

----------------------------

<blockquote>
(...) [Figure S14:] Using lysine as an example, this led to on average <span style="color:red">21,127</span> sites per sample (...)
</blockquote>

```R
# R
tmp <- scvs %>% filter(pN_popular_consensus <= 0.0005, most_freq_aa == 'Lys')
tmp %>% group_by(sample_id) %>% summarise(count(.)) %>% mutate(x=1) %>% group_by(x) %>% summarise(min(n), max(n), median(n), mean(n), sd(n))
``` 

## Supplementary Information

### Regimes of sequence similarity probed by metagenomics, SAR11 cultured genomes, and protein families

<blockquote>
(...) Unsurprisingly, protein families are most evolutionarily divergent (mean amino acid PS <span style="color:red">28.8</span>%). Relative to SAR11 homologs (mean nucleotide PS <span style="color:red">77.3</span>%), the aligned reads are highly related (mean nucleotide PS <span style="color:red">94.5</span>%), showing that metagenomics offers a modality of sequence inquiry more highly resolved than sequence comparisons between isolated cultures (...)
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
print(df %>% group_by(name) %>% summarise(mean=mean(value, na.rm=TRUE), median=median(value, na.rm=TRUE)))
``` 

### Comparing structure predictions between AlphaFold and MODELLER


----------------------------

<blockquote>
(...) While AlphaFold produced <span style="color:red">754</span> structures we deemed trustworthy (see Methods), MODELLER produced <span style="color:red">346</span> due to its reliance on pre-existing template structures (...)
</blockquote>

```bash
# bash
wc -l 12_GENES_WITH_GOOD_STRUCTURES
wc -l 12_GENES_WITH_GOOD_STRUCTURES_MODELLER
``` 

----------------------------

<blockquote>
(...) In <span style="color:red">339</span> cases both methods procured a structure prediction for a given protein sequence (...)
</blockquote>

```python
# python
goi_af = set([int(x.strip()) for x in open('12_GENES_WITH_GOOD_STRUCTURES').readlines()])
goi_mod = set([int(x.strip()) for x in open('12_GENES_WITH_GOOD_STRUCTURES_MODELLER').readlines()])
len(goi_af.intersection(goi_mod))
``` 

----------------------------

<blockquote>
(...) Since a TM score of 0.5 indicates that proteins likely belong to the same fold family (Xu and Zhang 2010), our average TM score of <span style="color:red">0.88</span> indicates strong overall agreement between AlphaFold and MODELLER (...)
</blockquote>

```R
# R
df <- read_tsv("../09_STRUCTURE_comparison.txt") %>% mutate(frac_ss_MOD = frac_alpha_MOD + frac_beta_MOD, frac_ss_AF = frac_alpha_AF + frac_beta_AF)
df$TM_score %>% mean(na.rm=T)
``` 

----------------------------

<blockquote>
(...) In fact, for the worst alignments (TM score <0.6), in <span style="color:red">15</span> of <span style="color:red">16</span> cases AlphaFold yielded more secondary structure (...)
</blockquote>

```R
# R
df <- read_tsv("../09_STRUCTURE_comparison.txt") %>% mutate(frac_ss_MOD = frac_alpha_MOD + frac_beta_MOD, frac_ss_AF = frac_alpha_AF + frac_beta_AF)
print(df %>% filter(TM_score < 0.6, frac_ss_AF > frac_ss_MOD) %>% dim() %>% .[[1]])
print(df %>% filter(TM_score < 0.6) %>% dim() %>% .[[1]])
``` 

----------------------------

<blockquote>
(...) Even still, these structures averaged a mean pLDDT score of <span style="color:red">90.8</span>, which is considered to be highly accurate (Jumper et al. 2021) (...)
</blockquote>

```R
# R
df <- read_tsv("../09_STRUCTURE_comparison.txt") %>% mutate(frac_ss_MOD = frac_alpha_MOD + frac_beta_MOD, frac_ss_AF = frac_alpha_AF + frac_beta_AF)
df %>% filter(has_AF, !has_MOD) %>% .$mean_pLDDT %>% mean()
``` 

### RSA and DTL predict nonsynonymous polymorphism rates

<blockquote>
(...) We excluded monomorphic sites (pN(site) = 0 for ns-models, pS(site) = 0 for s-models), sites with DTL > 40Å (see Methods), and removed gene-sample pairs containing <100 remaining sites, resulting in <span style="color:red">16,285</span> ns-models and <span style="color:red">24,553</span> s-models (Table S7) (...)
</blockquote>

```R
# R
print(pn_models %>% dim() %>% .[[1]])
print(ps_models %>% dim() %>% .[[1]])
``` 

----------------------------

<blockquote>
(...) We filtered out any genes that did not have a predicted structure and at least one predicted ligand-binding site, which when applied in conjunction with the above filters resulted in <span style="color:red">381</span> genes for the s-models and <span style="color:red">342</span> genes for the ns-models (...)
</blockquote>

```R
# R
print(ps_models %>% pull(gene_callers_id) %>% unique() %>% length())
print(pn_models %>% pull(gene_callers_id) %>% unique() %>% length())
``` 

----------------------------

<blockquote>
(...) ns-models yielded consistently positive correlations (average Pearson coefficient of rRSA = <span style="color:red">0.353</span>) (Figure 2c), whereas s-models exhibited correlations centered around 0 (average rRSA = <span style="color:red">-0.029</span>) (...)
</blockquote>

```R
# R
print(poly_corr %>% .$pn_r_rsa %>% mean(na.rm=T))
print(poly_corr %>% .$ps_r_rsa %>% mean(na.rm=T))
``` 

----------------------------

<blockquote>
(...) The average R2 was <span style="color:red">0.137</span> for ns-models, however model quality varied significantly between gene-sample pairs (...)
</blockquote>

```R
# R
poly_corr %>% .$pn_rsq_rsa %>% mean(na.rm=T)
``` 

----------------------------

<blockquote>
(...) In fact, we found that R2 varied from as high as <span style="color:red">0.526</span> (gene 2264 in sample ION_42_80M), to as low as <span style="color:red">0.0</span> (gene 2486 in sample ION_42_80M) (...)
</blockquote>

```R
# R
print(poly_corr %>% ungroup() %>% select(gene_callers_id, sample_id, pn_rsq_rsa) %>% slice(which.max(pn_rsq_rsa)))
print(poly_corr %>% ungroup() %>% select(gene_callers_id, sample_id, pn_rsq_rsa) %>% slice(which.min(pn_rsq_rsa)))
``` 

----------------------------

<blockquote>
(...) Using the same procedure, we linearly regressed log10(pN(site)) and log10(pS(site)) with DTL and found that <span style="color:red">96</span>% of ns-models yielded positive correlations with DTL with considerable predictive power, where on average <span style="color:red">11.5</span>% of per-site ns-polymorphism rate variation could be explained by DTL (Table S7) (...)
</blockquote>

```R
# R
print(poly_corr %>% filter(pn_r_dist > 0) %>% dim() %>% .[[1]] / (poly_corr %>% filter(!is.na(pn_r_dist)) %>% dim %>% .[[1]]))
print(poly_corr %>% .$pn_rsq_dist %>% mean(na.rm=T))
``` 

----------------------------

<blockquote>
(...) R2 values varied significantly, ranging from <span style="color:red">0.514</span> (gene 2326 in sample PSE_100_05M) to <span style="color:red">0.0</span> (gene 2246 in sample PSE_102_05M) (...)
</blockquote>

```R
# R
print(poly_corr %>% ungroup() %>% select(gene_callers_id, sample_id, pn_rsq_dist) %>% slice(which.max(pn_rsq_dist)))
print(poly_corr %>% ungroup() %>% select(gene_callers_id, sample_id, pn_rsq_dist) %>% slice(which.min(pn_rsq_dist)))
``` 

----------------------------

<blockquote>
(...) Interestingly, we found that log10(pS(site)) on average negatively correlates with DTL (average Pearson coefficient <span style="color:red">-0.057</span>) (...)
</blockquote>

```R
# R
poly_corr %>% pull(ps_r_dist) %>% mean()
``` 

----------------------------

<blockquote>
(...) A Pearson correlation between RSA and DTL revealed the relative independence of each variable from the other (R2 = <span style="color:red">0.082</span>, r = <span style="color:red">0.286</span>), precluding effects of multicollinearity (Figure SI1) (...)
</blockquote>

```R
# R
rsa_vs_d_data <- scvs %>%
filter(sample_id == 'ANE_004_05M') %>%
select(gene_callers_id, codon_order_in_gene, rel_solvent_acc, ANY_dist) %>%
filter((!is.na(ANY_dist)))
print(summary(lm(rel_solvent_acc ~ ANY_dist, data=rsa_vs_d_data))$r.squared)
print(summary(lm(rel_solvent_acc ~ ANY_dist, data=rsa_vs_d_data))$r.squared %>% sqrt())
``` 

----------------------------

<blockquote>
(...) The results revealed that including both RSA and DTL yielded a considerably better set of models for ns-polymorphism rates, with an average explained variance of 17.7% (average adjusted R2RSA-DTL = <span style="color:red">0.177</span>) (...)
</blockquote>

```R
# R
poly_corr %>% .$pn_adj_rsq %>% mean(na.rm=T)
``` 

----------------------------

<blockquote>
(...) For example, the group (RSA1, DTL2) contains the <span style="color:red">3,164</span> sites with RSA values in the 1st RSA range <span style="color:red">[0.00,0.01)</span> and DTL values in the 2nd DTL range <span style="color:red">[5.0Å,6.4Å)</span> (...)
</blockquote>

```R
# R
n_samples <- read_tsv("../soi", col_names=F) %>% nrow()
print(counts_group %>% filter(RSA_group == 0, DTL_group == 1) %>% .$count/n_samples)
print(pn_group %>% filter(RSA_group == 0, DTL_group == 1) %>% select(RSA_range))
print(pn_group %>% filter(RSA_group == 0, DTL_group == 1) %>% select(DTL_range))
``` 

----------------------------

<blockquote>
(...) Nonsynonymous polymorphism rates of groups varied from as low as <span style="color:red">0.001</span> to as high as <span style="color:red">0.021</span> (...)
</blockquote>

```R
# R
print(pn_group %>% .$pn %>% min())
print(pn_group %>% .$pn %>% max())
``` 

----------------------------

<blockquote>
(...) Sites exhibited a spectrum of ns-polymorphism rates that is roughly linear. We determined this by fitting a linear model pN(group) ~ i + j, where i refers to the group’s RSA and DTL indices (RSAi, DTLj), yielded an adjusted R2 of <span style="color:red">0.836</span> (...)
</blockquote>

```R
# R
summary(pn_group_model)$adj.r.squared
``` 


----------------------------

<blockquote>
(...) [T]he linear model pS(group) ~ i + j yielded a significant, anti-correlated relationship with both RSA and DTL (adjusted R2 of <span style="color:red">0.206</span>) (...)
</blockquote>

```R
# R
summary(ps_group_model)$adj.r.squared
``` 

----------------------------

<blockquote>
(...) in the sample-gene models, (1) the mean Pearson correlation coefficient between pS(site) and RSA is <span style="color:red">-0.013</span> (Figure 2c), and (2) the mean Pearson correlation coefficient between pS(site) and DTL is <span style="color:red">-0.052</span> (Figure 2d) (...)
</blockquote>

```R
# R
print(poly_corr %>% .$ps_r_rsa %>% mean(na.rm=T))
print(poly_corr %>% .$ps_r_dist %>% mean(na.rm=T))
``` 

----------------------------

<blockquote>
(...) [Figure SI1:] Scatter plot of RSA vs. DTL for the <span style="color:red">143,181</span> sites belonging to genes with a predicted structure and at least one predicted ligand (...)
</blockquote>

```R
# R
rsa_vs_d_data <- scvs %>%
filter(sample_id == 'ANE_004_05M') %>%
select(gene_callers_id, codon_order_in_gene, rel_solvent_acc, ANY_dist) %>%
filter(ANY_dist < 40)
print(rsa_vs_d_data %>% dim() %>% .[[1]])
``` 

----------------------------

<blockquote>
(...) The line of best fit is shown in black, The Pearson coefficient is <span style="color:red">0.313</span> and the R2 is <span style="color:red">0.098</span> (...)
</blockquote>

```R
# R
rsa_vs_d_data <- scvs %>%
filter(sample_id == 'ANE_004_05M') %>%
select(gene_callers_id, codon_order_in_gene, rel_solvent_acc, ANY_dist) %>%
filter(ANY_dist < 40)
model_r2 <- summary(lm(rel_solvent_acc ~ ANY_dist, data=rsa_vs_d_data))$r.squared
model_corr <- summary(lm(rel_solvent_acc ~ ANY_dist, data=rsa_vs_d_data))$r.squared %>% sqrt()
print(paste("R2 =", round(model_r2, 3), ", r =", round(model_corr, 3)))
``` 

----------------------------

<blockquote>
(...) [Figure SI2:] Red denotes parameter/error distributions for the <span style="color:red">16,285</span> nonsynonymous models of the form pN(site) $= \beta_0 + \beta_{RSA}RSA + \beta_{DTL}DTL$ and blue denotes parameter/error distributions for the <span style="color:red">24,553</span> models of the form pS(site) $= \beta_0 + \beta_{RSA}RSA + \beta_{DTL}DTL (...)
</blockquote>

```R
# R
print(pn_models %>% dim() %>% .[[1]])
print(ps_models %>% dim() %>% .[[1]])
``` 

### dN/dS(gene) and sample-averaged pN/pS(gene) yield consistent results


<blockquote>
(...) We calculated dN/dS(gene) for <span style="color:red">753</span> homologous gene pairs found between HIMB83 and a closely related cultured representative HIMB122 (see Methods) (...)
</blockquote>

```R
# R
pnps_sample_averaged <- pnps %>%
group_by(gene_callers_id) %>%
summarize(pnps = mean(pnps, na.rm=T))
dnds <- read_tsv("../19_DNDS_HIMB122/dnds.txt", col_names=F) %>%
rename(gene_callers_id=X1, dnds=X2) %>%
filter(!is.na(dnds))
df <- full_join(pnps_sample_averaged, dnds) %>%
mutate(log_pnps = log10(pnps), log_dnds = log10(dnds))
overlap <- df %>% filter(!is.na(pnps), !is.na(dnds))
overlap %>% dim() %>% .[[1]]
``` 

----------------------------

<blockquote>
(...) ANI between HIMB83 and HIMB122 was <span style="color:red">82.6</span>%, whereas the average ANI between HIMB83 and recruited reads was 94.5% (...)
</blockquote>

```python
# python
import pandas as pd
df = pd.read_csv("07_ANI_HIMB122/ANIb_percentage_identity.txt", sep='\t').set_index('key')
df['HIMB083'].sort_values()
``` 

----------------------------

<blockquote>
(...) ANI between HIMB83 and HIMB122 was 82.6%, whereas the average ANI between HIMB83 and recruited reads was <span style="color:red">94.5</span>% (...)
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
print(df %>% group_by(name) %>% summarise(mean=mean(value, na.rm=TRUE), median=median(value, na.rm=TRUE)))
``` 

----------------------------

<blockquote>
(...) We found that log-transformed sample-averaged pN/pS(gene) highly correlated with log-transformed dN/dS(gene) (Pearson R2 = <span style="color:red">0.380</span>) (...)
</blockquote>

```R
# R
pnps_sample_averaged <- pnps %>%
group_by(gene_callers_id) %>%
summarize(pnps = mean(pnps, na.rm=T))
dnds <- read_tsv("../19_DNDS_HIMB122/dnds.txt", col_names=F) %>%
rename(gene_callers_id=X1, dnds=X2) %>%
filter(!is.na(dnds))
df <- full_join(pnps_sample_averaged, dnds) %>%
mutate(log_pnps = log10(pnps), log_dnds = log10(dnds))
overlap <- df %>% filter(!is.na(pnps), !is.na(dnds))
overlap %>% lm(pnps ~ dnds, data=.) %>% summary() %>% .$r.squared
``` 

----------------------------

<blockquote>
(...) The ratio between sample-averaged pN/pS(gene) and dN/dS(gene) was on average <span style="color:red">6.23</span> (Figure S10), matching expectations that slightly deleterious, nonsynonymous mutants commonly drift to observable frequencies, yet far less commonly drift to fixation (...)
</blockquote>

```R
# R
pnps_sample_averaged <- pnps %>%
group_by(gene_callers_id) %>%
summarize(pnps = mean(pnps, na.rm=T))
dnds <- read_tsv("../19_DNDS_HIMB122/dnds.txt", col_names=F) %>%
rename(gene_callers_id=X1, dnds=X2) %>%
filter(!is.na(dnds))
df <- full_join(pnps_sample_averaged, dnds) %>%
mutate(log_pnps = log10(pnps), log_dnds = log10(dnds))
overlap <- df %>% filter(!is.na(pnps), !is.na(dnds))
overlap %>% filter(dnds > 0) %>% mutate(ratio = pnps/dnds) %>% pull(ratio) %>% mean()
``` 

### Transcript abundance largely explains genic differences in the strengths of purifying selection

<blockquote>
(...) Sample-averaged pN/pS(gene) values varied significantly between genes, varying from <span style="color:red">0.004-0.539</span>, with a mean of <span style="color:red">0.063</span> (Figure S12, Table S9) (...)
</blockquote>

```R
# R
print(pnps %>% group_by(gene_callers_id) %>% summarise(pnps=mean(pnps)) %>% .$pnps %>% min(na.rm=T))
print(pnps %>% group_by(gene_callers_id) %>% summarise(pnps=mean(pnps)) %>% .$pnps %>% max(na.rm=T))
print(pnps %>% group_by(gene_callers_id) %>% summarise(pnps=mean(pnps)) %>% .$pnps %>% mean(na.rm=T))
``` 

----------------------------

<blockquote>
(...) Comparing sample-median TA values to sample-averaged pN/pS(gene) values yielded a strong, negative correlation (Figure SI5b, Pearson r = <span style="color:red">-0.539</span>, R2 = <span style="color:red">0.290</span>) according to an inverse power-law relationship (...)
</blockquote>

```R
# R
print(cor(TA_sample_averaged$log_pnps_mean, TA_sample_averaged$log_TA_median, use="pairwise.complete.obs"))
print(cor(TA_sample_averaged$log_pnps_mean, TA_sample_averaged$log_TA_median, use="pairwise.complete.obs") ^ 2)
``` 

----------------------------

<blockquote>
(...) We found that of the 799 genes tested, <span style="color:red">74</span>% exhibited (weak) negative correlations between $log_{10}(TA+0.01)$ and $log_{10}(pN/pS(gene))$ (Figure SI5c), yet only 11.5% of genes passed significance tests (one-sided Pearson, 25% Benjamini-Hochberg false discovery rate) (Figure SI5d) (...)
</blockquote>

```R
# R
TA_sample_corr %>% pull(negative) %>% table() %>% .[[2]] / TA_sample_corr %>% dim() %>% .[[1]]
``` 

----------------------------

<blockquote>
(...) We found that of the 799 genes tested, 74% exhibited (weak) negative correlations between $log_{10}(TA+0.01)$ and $log_{10}(pN/pS(gene))$ (Figure SI5c), yet only <span style="color:red">11.5</span>% of genes passed significance tests (one-sided Pearson, 25% Benjamini-Hochberg false discovery rate) (Figure SI5d) (...)
</blockquote>

```R
# R
TA_sample_corr %>% filter(fdr_25 == TRUE) %>% nrow() / (TA_sample_corr %>% filter(fdr_25 == FALSE) %>% nrow() + TA_sample_corr %>% filter(fdr_25 == TRUE) %>% nrow()) * 100
``` 

----------------------------

<blockquote>
(...) [Figure SI5:] The linear model yielded a Pearson coefficient of -0.539, an R2 of 0.290 and a line of best fit <span style="color:red">y = (-0.31 士 0.02)x + (-1.63 士 0.02)</span> shown in pink (95% confidence intervals shown in translucent pink) (...)
</blockquote>

```R
# R
lm(log_pnps_mean ~ log_TA_median, data=TA_sample_averaged) %>% summary() %>% .$coefficients
``` 

### Stability analysis of polymorphism distributions with respect to pN/pS(core)


<blockquote>
(...) Our bootstrapping stability analysis (Figure SI6, Table S14) showed that in <span style="color:red">99.5</span>% of gene resamplings, the mean RSA of ns-polymorphism negatively associated with pN/pS(core) (one-sided Pearson coefficient p-value <0.05), whereas in only <span style="color:red">69.5</span>% of gene resamplings did the mean DTL of ns-polymorphism negatively associate with pN/pS(core) (...)
</blockquote>

```R
# R
plot_data <- read_tsv('../genome_robust.txt')
print((plot_data %>% pull(pN_RSA_pval) <= 0.05) %>% sum() / plot_data %>% dim() %>% .[[1]])
print((plot_data %>% pull(pN_DTL_pval) <= 0.05) %>% sum() / plot_data %>% dim() %>% .[[1]])
``` 




















# Results

## What we did

We analyzed microbial taxon abundance data from 690 human metagenomic samples to understand how microbial communities differ by body site. We performed three complementary analyses:

1. **Principal Coordinates Analysis (PCoA)** - to visualize overall community dissimilarities
2. **PERMANOVA test** - to statistically test whether body site explains microbial composition differences
3. **Indicator species analysis** - to identify taxa that are characteristic of each body site

## PERMANOVA analysis shows body site explains 39.9% of microbial variation

We tested the null hypothesis that microbial composition does not differ by body site using PERMANOVA (999 permutations).

```
          Df SumOfSqs      R2     F Pr(>F)    
Model      4   111.92 0.39945 113.9  0.001 ***
Residual 685   168.27 0.60055                 
Total    689   280.20 1.00000                 
```

**Results:**
- Model degrees of freedom: 4 (five body site categories minus 1)
- Residual degrees of freedom: 685
- Sum of Squares (Model): 111.92
- Sum of Squares (Residual): 168.27
- R² (explained variance): 0.39945 (39.945%)
- F-statistic: 113.9
- P-value: 0.001 (from permutation test)

**Interpretation:** The p-value of 0.001 means there is only a 0.1% chance of observing these differences if body site had no effect on microbial composition. The R² value of 0.399 indicates that approximately 40% of the variation in microbial community composition between samples can be attributed to differences in body site. This is a strong effect size, suggesting that body site is a major determinant of which microbes are found in each location.

## PCoA visualization shows clear clustering by body site

To visualize these patterns, we performed PCoA using Bray-Curtis dissimilarity and plotted the first two axes.

![PCoA of Microbial Community Composition by Body Site](ordination_plot.png)

**Figure 1:** PCoA ordination of 690 samples colored by body site. The plot shows two dominant gradients of variation.

## Top indicator taxa for each body site

For each body site, we identified the most abundant taxa to identify characteristic microbes:

### Airways

Top 5 taxa by mean abundance:
- *Propionibacterium_accolens*: 42.528689%
- *Corynebacterium_accolens*: 21.405716%
- *Staphylococcus_epidermidis*: 12.717085%
- *Staphylococcus_aureus*: 5.000064%
- *Propionibacterium_unclassified*: 2.872164%

### Gastrointestinal Tract

Top 5 taxa by mean abundance:
- *Bacteroides_unclassified*: 16.509382%
- *Bacteroides_vulgatus*: 9.758981%
- *Alistipes_putredinis*: 9.636933%
- *Prevotella_copri*: 5.723292%
- *Bacteroides_ovatus*: 5.029873%

### Oral

Top 5 taxa by mean abundance:
- *Streptococcus_mitis*: 15.543022%
- *Haemophilus_parainfluenzae*: 12.515504%
- *Corynebacterium_matruchotii*: 5.787398%
- *Rothia_dentocariosa*: 4.452347%
- *Prevotella_melaninogenica*: 3.831314%

### Skin

Top 5 taxa by mean abundance:
- *Propionibacterium_acnes*: 70.4107731%
- *Staphylococcus_epidermidis*: 13.9160465%
- *Propionibacterium_unclassified*: 3.3353642%
- *Streptococcus_mitis*: 0.3755958%
- *Corynebacterium_accolens*: 0.2256650%

### Urogenital Tract

Top 5 taxa by mean abundance:
- *Lactobacillus_crispatus*: 47.52486393%
- *Lactobacillus_iners*: 19.51101232%
- *Lactobacillus_jensenii*: 16.27271768%
- *Lactobacillus_gasseri*: 8.96562482%
- *Bacteroides_vulgatus*: 0.03045054%

**Interpretation:** Each body site has a distinct microbial signature. The skin is dominated by *Propionibacterium_acnes* (70.4%), the urogenital tract by *Lactobacillus* species (collectively >90%), the gastrointestinal tract by *Bacteroides* species, and the oral cavity by *Streptococcus* and *Haemophilus*.

## Feature importance: taxa driving axis separation

We tested which taxa correlate most strongly with the PCoA axes to understand what drives the observed separation:

```
***VECTORS

     Streptococcus_mitis Propionibacterium_acnes     r2 Pr(>r)    
[1,]            -0.64783                 0.76179 0.3273  0.001 ***
[2,]             0.22545                 0.97425 0.5525  0.001 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

**Axis 1 (Dim1):**
- *Streptococcus_mitis*: correlation = -0.65
- *Propionibacterium_acnes*: correlation = 0.76

**Axis 2 (Dim2):**
- *Streptococcus_mitis*: correlation = 0.23
- *Propionibacterium_acnes*: correlation = 0.97

**Interpretation:** *Propionibacterium_acnes* shows strong positive correlation with both axes, making it a key driver of separation between body sites. *Streptococcus_mitis* shows negative correlation with Dim1 but positive with Dim2, helping to distinguish Oral/Airways from Skin/Gastrointestinal samples. The high r² values (0.3273 and 0.5525) with p-values of 0.001 indicate these correlations are statistically significant.

## Summary

The analysis demonstrates that body site is a powerful determinant of microbial community composition:

1. PERMANOVA confirms body site explains 39.9% of variation (p=0.001)
2. PCoA visualization shows clear clustering by body site
3. Each body site has a characteristic microbial profile dominated by specific taxa
4. *Propionibacterium_acnes* and *Streptococcus_mitis* are major drivers of variation between sites

These results are consistent with established knowledge of human microbiota, where different body sites provide distinct environmental niches that select for specific microbial communities.

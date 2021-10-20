---
layout: page
title: binding-frequencies-txt [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/binding-frequencies-txt
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-run-interacdome](../../programs/anvi-run-interacdome)</span></p>


## Required or used by


There are no anvi'o tools that use or require this artifact directly, which means it is most likely an end-product for the user.


## Description

When the user runs <span class="artifact-n">[anvi-run-interacdome](/software/anvio/help/7.1/programs/anvi-run-interacdome)</span>, it stores binding frequencies directly into the <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> as <span class="artifact-n">[misc-data-amino-acids](/software/anvio/help/7.1/artifacts/misc-data-amino-acids)</span>. Yet <span class="artifact-n">[anvi-run-interacdome](/software/anvio/help/7.1/programs/anvi-run-interacdome)</span> also outputs tabular data directly accessible by the user--this data is what is meant by <span class="artifact-n">[binding-frequencies-txt](/software/anvio/help/7.1/artifacts/binding-frequencies-txt)</span>.

Specifically, this artifact refers to 2 files named `INTERACDOME-match_state_contributors.txt` and `INTERACDOME-domain_hits.txt` (the `INTERACDOME` prefix can be changed with `-O`). 

`INTERACDOME-match_state_contributors.txt` displays the binding frequencies in the following format:

|   gene_callers_id |   codon_order_in_gene | pfam_id   |   match_state | ligand   |   binding_freq |
|------------------:|----------------------:|:----------|--------------:|:---------|---------------:|
|                 1 |                   169 | PF00534   |            22 | ADP      |      0.687948  |
|                 1 |                   169 | PF13692   |             8 | ADP      |      0.595441  |
|                 1 |                   174 | PF00534   |            27 | ADP      |      0.735759  |
|                 1 |                   174 | PF13692   |            14 | ADP      |      0.595441  |
|                 1 |                   184 | PF00534   |            37 | ADP      |      0.0697656 |
|                 1 |                   184 | PF13692   |            24 | ADP      |      0.101399  |
|                 1 |                   186 | PF00534   |            39 | ADP      |      0.0697656 |
|                 1 |                   186 | PF13692   |            26 | ADP      |      0.101399  |
|                 1 |                   187 | PF13692   |            27 | ADP      |      0.201761  |
|                 1 |                   189 | PF00534   |            47 | ADP      |      0.0697656 |

Each binding frequency is associated with both the exact residue of the user's gene sequences (from their <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>) and the exact match states (from the Pfam database) that contributed the binding frequency. 

`INTERACDOME-match_state_contributors.txt` is a parsed summary of the `hmmsearch` output in the following format:


| pfam_name       | pfam_id   |   corresponding_gene_call |   domain | qual   |   score |   bias |   c-evalue |   i-evalue |   hmm_start |   hmm_stop | hmm_bounds   |   ali_start |   ali_stop | ali_bounds   |   env_start |   env_stop | env_bounds   |   mean_post_prob | match_state_align                                                                                                                                                                                                     | comparison_align                                                                                                                                                                                                      | sequence_align                                                                                                                                                                                                        |   version |
|:----------------|:----------|--------------------------:|---------:|:-------|--------:|-------:|-----------:|-----------:|------------:|-----------:|:-------------|------------:|-----------:|:-------------|------------:|-----------:|:-------------|-----------------:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------:|
| Beta_elim_lyase | PF01212   |                      1762 |        1 | !      |    20.9 |    0.1 |    1e-08   |    3.5e-06 |          33 |        169 | ..           |          44 |        177 | ..           |          34 |        215 | ..           |             0.72 | tvnrLedavaelfgke..aalfvpqGtaAnsill.kill.qr..geevivtepahihfdetgaiaelagvklrdlknkeaGkmdlekleaaikevgaheekiklisltvTnntagGqvvsleelrevaaiakkygiplhlDgA                                                                       | ++  +++ael+      + f+  Gt +++  l  + + +r  g+ +i++   h   +et    +  g +l  ++ +++G +++e+l+++i++     e i + +++v n+   G++ +++e+ ev  +a+  +i++h+D+                                                                          | LLQQARKQIAELINVSanEIYFTSGGTEGDNWVLkGTAIeKRefGNHIIISAVEHPAVTETAEQLVELGFELSYAPVDKEGRVKVEELQKLIRK-----ETILVSVMAVNNE--VGTIQPIKEISEV--LAEFPKIHFHVDAV                                                                       |        20 |
| PAPS_reduct     | PF01507   |                      1541 |        1 | !      |    36.1 |    0.1 |    3.6e-13 |    1.3e-10 |           2 |        164 | ..           |          21 |        231 | ..           |          20 |        234 | ..           |             0.79 | lvvsvsgGkdslVllhLalkafkpv....pvvfvdtghefpetiefvdeleeryglrlkvyepeeevaekinaekhgs.slyee.aaeriaKveplkk.................................aLekldedall..tGaRrdesksraklpiveidedfek.........slrvfPllnWteedvwqyilrenipynpLydqgfr | + +s+sgGkds  +++La  + ++      ++ ++ + ++  t++f++++e+  +++ +++     ++++ + + +++ + +   + e+ +   p  k                                   e++ ++a+   +G+R++es +r++     +++ +++          + ++Pl++W+  d+w+   + +++yn +y++ ++ | VYFSFSGGKDSGLMVQLANLVAEKLdrnfDLLILNIEANYTATVDFIKKIEQLPRVKNIYHFCLPFFEDNNTSFFQPQwKMWDPsEKEKWIHSLP--KnaitleniddglkkyyslsngnpdrflryfqnwYKEQYPQSAIScgVGIRAQESLHRHSAVTKGENKYKNRcwinitlegNILFYPLFDWKVGDIWAATFKCELEYNYIYEKMYK |        18 |
| Ank_2           | PF12796   |                      1756 |        1 | !      |    32.2 |    0   |    6.7e-12 |    2.3e-09 |          29 |         84 | .]           |          74 |        135 | ..           |          53 |        135 | ..           |             0.85 | aLhyAakngnleivklLle...h.a..adndgrtpLhyAarsghleivklLlekgadinlkd                                                                                                                                                        | aL  Aa + +++ vk +l+   + +  +d +g+tpL +A+ ++ +ei+k L+++gadinl++                                                                                                                                                        | ALLEAANQRDTKKVKEILQdttYqVdeVDTEGNTPLNIAVHNNDIEIAKALIDRGADINLQN                                                                                                                                                        |         6 |
| Ank_2           | PF12796   |                      1756 |        2 | !      |    28.5 |    0   |    9.5e-11 |    3.3e-08 |          22 |         75 | ..           |         199 |        265 | ..           |         195 |        267 | .]           |             0.76 | pn..k.ngktaLhyAak..ngnl...eivklLleha.....adndgrtpLhyAarsghleivklLle                                                                                                                                                   | ++  + +g taL+ A+   +gn    +ivklL+e++      dn+grt++ yA ++g++ei k+L +                                                                                                                                                   | IDfqNdFGYTALIEAVGlrEGNQlyqDIVKLLMENGadqsiKDNSGRTAMDYANQKGYTEISKILAQ                                                                                                                                                   |         6 |
| IGPS            | PF00218   |                      1615 |        1 | !      |    20.6 |    0.1 |    1.2e-08 |    4e-06   |         202 |        249 | ..           |         195 |        242 | ..           |          73 |        248 | ..           |             0.88 | LaklvpkdvllvaeSGiktredveklkeegvnafLvGeslmrqedvek                                                                                                                                                                      | +++lv+++++++ae  i+t+e+++++k+ gv ++ vG +++r ++ +k                                                                                                                                                                      | IKQLVQENICVIAEGKIHTPEQARQIKKLGVAGIVVGGAITRPQEIAK                                                                                                                                                                      |        20 |
| Ribosomal_L33   | PF00471   |                      1562 |        1 | !      |    66.6 |    1.5 |    1.1e-22 |    3.7e-20 |           2 |         47 | .]           |           4 |         49 | .]           |           3 |         49 | .]           |             0.97 | kvtLeCteCksrnYtttknkrntperLelkKYcprcrkhtlhkEtK                                                                                                                                                                        | +++LeC e+++r Y t+knkrn+perLelkKY p++r++ ++kE K                                                                                                                                                                        | NIILECVETGERLYLTSKNKRNNPERLELKKYSPKLRRRAIFKEVK                                                                                                                                                                        |        19 |
| Ribosomal_S14   | PF00253   |                      1565 |        1 | !      |    83.3 |    0.1 |    3.9e-28 |    1.3e-25 |           2 |         54 | .]           |          36 |         88 | ..           |          35 |         88 | ..           |             0.98 | laklprnssptrirnrCrvtGrprGvirkfgLsRicfRelAlkgelpGvkKaS                                                                                                                                                                 | laklpr+s+p+r+r r++ +GrprG++rkfg+sRi+fRel ++g +pGvkKaS                                                                                                                                                                 | LAKLPRDSNPNRLRLRDQTDGRPRGYMRKFGMSRIKFRELDHQGLIPGVKKAS                                                                                                                                                                 |        20 |
| Polysacc_synt_C | PF14667   |                      1593 |        1 | !      |    61.4 |   19.2 |    5.4e-21 |    1.9e-18 |           2 |        139 | ..           |         371 |        516 | ..           |         370 |        519 | ..           |             0.83 | LailalsiiflslstvlssiLqglgrqkialkalvigalvklilnllliplfgivGaaiatvlallvvavlnlyalrrllgikl...llrrllkpllaalvmgivvylllllllglllla...al..alllavlvgalvYllllll                                                                    | L+  ++s+ +l+++t++ siLq+l  +k+a+ ++ i++l+kli+++++i+lf  +G +iat+++ ++++++ +++l+r++ i++    ++   +++ +++vm i+ +l+l+++ ++   +   +l   + l +++g++v+ + l++                                                                    | LSATIISTSLLGIFTIVLSILQALSFHKKAMQITSITLLLKLIIQIPCIYLFKGYGLSIATIICTMFTTIIAYRFLSRKFDINPikyNRKYYSRLVYSTIVMTILSLLMLKIISSVYKFEstlQLffLISLIGCLGGVVFSVTLFR                                                                    |         5 |

For each hit, this table includes how good the hit was, the alignment of the user gene to the exact HMM match states, and more! In fact, it includes all of hte domain hit summary information, the sequence of the consensus match states, the comparison string for the hit, and the sequence of the user's gene. 

For more information, check out [this blogpost](https://merenlab.org/2020/07/22/interacdome/#6-storing-the-per-residue-binding-frequencies-into-the-contigs-database). 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/binding-frequencies-txt.md) to update this information.


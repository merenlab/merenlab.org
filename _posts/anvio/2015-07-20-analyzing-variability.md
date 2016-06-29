---
layout: post
authors: [meren]
title: "Analyzing variability with anvi'o"
excerpt: "Exploring variability patterns using anvi'o for deeper insights"
modified: 2015-07-20
tags: []
categories: [anvio]
comments: true
---

{% include _toc.html %}

*De novo* characterization of variable nucleotide sites in a given genome bin (or in a given contig, or even in a given gene) is one of the most powerful aspects of anvi'o. This ability opens doors to the use of subtle differences across environments. Through simple programs anvi'o installs, it is possible to focus on a small number of nucleotide positions, and utilize them to infer ecology at a very highly resolved manner.

# What is variation? Why don't we call them SNPs?

I would like to define variation as the extent of disagreement between the aligned nucleotides from short reads that map to the same context.

I do acknowledge that it would be much easier to communicate this entire thing if we called them SNPs. But the term 'SNP' [has a very specific and well-established definition](https://en.wikipedia.org/wiki/Single-nucleotide_polymorphism), and I think it is being a bit overused, especially in microbiology, and even more specifically in metagenomics, where the origin of variation is not always quite clear to us due to the overwhelming complexity of what we sequence, and the bioinformatics approaches we employ to make sense of this data. I do not want to lose you early on, so let's postpone this discussion for another time.

Depending on the complexity in a community, short reads in a metagenomic dataset that map to a contig can generate one or more mismatches if the mapping does not require a 100% sequence identity. The source of a mismatch may be artificial, such as stochastic sequencing or PCR errors, while others may represent ecologically informative variation.

Although the importance of it is probably almost common sense to anyone who understands the steps we follow to recover our contigs for a given environment, and although many amazing studies exist including the ones from Jill Banfield's lab that made use of SNP patterns in metagenomic settings such as [this one](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0060177), working with subtle nucleotide variaiton in metagenomic bins hasn't been so practical to access and/or utilize them to make better sense of the ecology of microbes. There are **many** challenges to address if one wants to focus on "subtle" things in these complex datasets, and in order to set the stage with reasonable expectations I would like to clarify that anvi'o will not exactly offer you a magical wand. However, anvi'o _will_ empower you to work with this aspect of your datasets by making the recovery and query of this information more or less 'straightforward'. I hope at the end of this post it will be clear to you how does it do it.

---

Just so you have a mental picture of how anvi'o visualizes these varaible nucleotide positions, here is a polished screenshot from the interface:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2015-07-20-analyzing-variability/example.png"><img src="{{ site.url }}/images/anvio/2015-07-20-analyzing-variability/example.png" style="border: none; width: 100%;" /></a>
</div>

What you see is the coverage and base frequencies in reported variable positions in a single contig accross two samples. Where these samples are coming from are not really relevant to the rest of this post. But usually, just like it is the case in this example, variable nucleotide positions and their base frequencies do not seem to be random, and they can be extremely reproducible if the population structures are really similar in different environments.

# The anvi'o way

Anvi'o gives access to variable nucleotide information in two steps.

The first step is to identify variable nucleotide positions and reporting them _for each sample separately_ (which takes place during the profiling step of the analysis (see [the tutorial]({% post_url anvio/2016-07-22-anvio-tutorial-v2 %}) if you are not familiar with the metagenomic workflow of anvi'o)).

The second step is to interpret the ecological significance of sample-specific variable positions _across samples_. The first step is agnostic to the experimental design, and/or the variable nucleotide positions in other samples: it simply does its best to find and report variation. The second step is where the user's input about the experimental design, and their questions come into the play.

Next two chapters in this article will detail these two steps, and the third chapter will explain the analysis of variability in the infant gut time-series dataset [from this study](http://genome.cshlp.org/content/23/1/111.long).

Everything you might need to replicate this analysis and to re-generate the third figure in the anvi'o paper is [here]({{ site.url }}/data/). I hope I will manage to give enough details that you will be able to modify my recipe and start exploring patterns in your own datasets if you wish.

Please do not hesitate to ask questions, make suggestions, and/or start discussions regarding this topic here in the comments section down below, or at [the issues section of the anvi'o repository on Github](https://github.com/meren/anvio/issues), or by directly contacting [us]({{ site.url }}/people/).

# *De novo* characterization and reporting of nucleotide variation

During the profiling step that should be done for each sample separately, anvi’o checks the composition of nucleotides that map to the same nucleotide position in a given contig, and keeps track of variation in the form of 'base frequencies'.

Profiler stores this information in the 'profile database' for each sample. and when samples are merged with `anvi-merge`, all reported variable positions from individual samples are also merged into one table. This table is called `variable_positions`, and the structure of it is identical in single and merged profiles. Here is a quick look into a `variable_positions` table (in this example I am showing only the top 30 most highly covered nucleotide positions in the merged infant gut data, which I am going to detail in the last chapter):

 <table id="table_layers">
     <thead>
         <tr>
             <th>entry_id</th>
             <th>sample_id</th>
             <th>split_name</th>
             <th>pos</th>
             <th>coverage</th>
             <th>n2n1ratio</th>
             <th>competing_nts</th>
             <th>consensus</th>
             <th>A</th>
             <th>T</th>
             <th>C</th>
             <th>G</th>
             <th>N</th>
         </tr>
     </thead>
     <tbody>
        <tr><td>24165</td><td>DAY_19</td><td>Day17a_QCcontig25_split_00001</td><td>577</td><td>5664</td><td>0.990838618745596</td><td>GC</td><td>G</td><td>7</td><td>7</td><td>2812</td><td>2838</td><td>0</td></tr>
        <tr><td>24909</td><td>DAY_19</td><td>Day17a_QCcontig14_split_00001</td><td>27562</td><td>5599</td><td>0.929509329647547</td><td>TC</td><td>T</td><td>5</td><td>2894</td><td>2690</td><td>10</td><td>0</td></tr>
        <tr><td>24167</td><td>DAY_19</td><td>Day17a_QCcontig25_split_00001</td><td>808</td><td>5389</td><td>0.799732530926112</td><td>TG</td><td>T</td><td>3</td><td>2991</td><td>3</td><td>2392</td><td>0</td></tr>
        <tr><td>24170</td><td>DAY_19</td><td>Day17a_QCcontig25_split_00001</td><td>1145</td><td>5307</td><td>0.635352286773795</td><td>AT</td><td>A</td><td>3236</td><td>2056</td><td>7</td><td>8</td><td>0</td></tr>
        <tr><td>24166</td><td>DAY_19</td><td>Day17a_QCcontig25_split_00001</td><td>835</td><td>5202</td><td>0.805429864253394</td><td>CT</td><td>C</td><td>5</td><td>2314</td><td>2873</td><td>10</td><td>0</td></tr>
        <tr><td>37483</td><td>DAY_22A</td><td>Day17a_QCcontig25_split_00001</td><td>1145</td><td>5162</td><td>0.625590923416325</td><td>AT</td><td>A</td><td>3173</td><td>1985</td><td>2</td><td>2</td><td>0</td></tr>
        <tr><td>37476</td><td>DAY_22A</td><td>Day17a_QCcontig25_split_00001</td><td>808</td><td>5007</td><td>0.922367409684858</td><td>TG</td><td>T</td><td>1</td><td>2602</td><td>4</td><td>2400</td><td>0</td></tr>
        <tr><td>39469</td><td>DAY_22A</td><td>Day17a_QCcontig14_split_00001</td><td>27562</td><td>4981</td><td>0.984455958549223</td><td>TC</td><td>T</td><td>1</td><td>2509</td><td>2470</td><td>1</td><td>0</td></tr>
        <tr><td>37472</td><td>DAY_22A</td><td>Day17a_QCcontig25_split_00001</td><td>835</td><td>4922</td><td>0.940828402366864</td><td>CT</td><td>C</td><td>2</td><td>2385</td><td>2535</td><td>0</td><td>0</td></tr>
        <tr><td>17042</td><td>DAY_18</td><td>Day17a_QCcontig25_split_00001</td><td>808</td><td>4738</td><td>0.987830465799412</td><td>TG</td><td>T</td><td>1</td><td>2383</td><td>0</td><td>2354</td><td>0</td></tr>
        <tr><td>37471</td><td>DAY_22A</td><td>Day17a_QCcontig25_split_00001</td><td>577</td><td>4688</td><td>0.982233502538071</td><td>CG</td><td>C</td><td>0</td><td>2</td><td>2364</td><td>2322</td><td>0</td></tr>
        <tr><td>17049</td><td>DAY_18</td><td>Day17a_QCcontig25_split_00001</td><td>1145</td><td>4643</td><td>0.562815762883126</td><td>AT</td><td>A</td><td>2969</td><td>1671</td><td>2</td><td>1</td><td>0</td></tr>
        <tr><td>17038</td><td>DAY_18</td><td>Day17a_QCcontig25_split_00001</td><td>835</td><td>4557</td><td>0.973570190641248</td><td>CT</td><td>C</td><td>2</td><td>2247</td><td>2308</td><td>0</td><td>0</td></tr>
        <tr><td>18021</td><td>DAY_18</td><td>Day17a_QCcontig14_split_00001</td><td>27562</td><td>4461</td><td>0.965182899955928</td><td>CT</td><td>C</td><td>1</td><td>2190</td><td>2269</td><td>1</td><td>0</td></tr>
        <tr><td>17037</td><td>DAY_18</td><td>Day17a_QCcontig25_split_00001</td><td>577</td><td>4400</td><td>0.890658631080499</td><td>CG</td><td>C</td><td>7</td><td>1</td><td>2323</td><td>2069</td><td>0</td></tr>
        <tr><td>52211</td><td>DAY_16</td><td>Day17a_QCcontig14_split_00001</td><td>27562</td><td>4057</td><td>0.913597733711048</td><td>TC</td><td>T</td><td>1</td><td>2118</td><td>1935</td><td>3</td><td>0</td></tr>
        <tr><td>51254</td><td>DAY_16</td><td>Day17a_QCcontig25_split_00001</td><td>808</td><td>3929</td><td>0.908560311284047</td><td>TG</td><td>T</td><td>4</td><td>2056</td><td>1</td><td>1868</td><td>0</td></tr>
        <tr><td>24168</td><td>DAY_19</td><td>Day17a_QCcontig25_split_00001</td><td>15856</td><td>3825</td><td>0.0774011299435028</td><td>GT</td><td>G</td><td>6</td><td>274</td><td>5</td><td>3540</td><td>0</td></tr>
        <tr><td>51261</td><td>DAY_16</td><td>Day17a_QCcontig25_split_00001</td><td>1145</td><td>3765</td><td>0.595672464997879</td><td>AT</td><td>A</td><td>2357</td><td>1404</td><td>3</td><td>1</td><td>0</td></tr>
        <tr><td>17046</td><td>DAY_18</td><td>Day17a_QCcontig25_split_00001</td><td>15856</td><td>3734</td><td>0.0739355581127733</td><td>GT</td><td>G</td><td>0</td><td>257</td><td>1</td><td>3476</td><td>0</td></tr>
        <tr><td>51250</td><td>DAY_16</td><td>Day17a_QCcontig25_split_00001</td><td>835</td><td>3699</td><td>0.922557172557173</td><td>CT</td><td>C</td><td>0</td><td>1775</td><td>1924</td><td>0</td><td>0</td></tr>
        <tr><td>51249</td><td>DAY_16</td><td>Day17a_QCcontig25_split_00001</td><td>577</td><td>3648</td><td>0.955984970477724</td><td>GC</td><td>G</td><td>2</td><td>2</td><td>1781</td><td>1863</td><td>0</td></tr>
        <tr><td>51258</td><td>DAY_16</td><td>Day17a_QCcontig25_split_00001</td><td>15856</td><td>3222</td><td>0.0745071834279987</td><td>GT</td><td>G</td><td>4</td><td>223</td><td>2</td><td>2993</td><td>0</td></tr>
        <tr><td>30612</td><td>DAY_15B</td><td>Day17a_QCcontig14_split_00001</td><td>27562</td><td>3179</td><td>0.763333333333333</td><td>TC</td><td>T</td><td>2</td><td>1800</td><td>1374</td><td>3</td><td>0</td></tr>
        <tr><td>37480</td><td>DAY_22A</td><td>Day17a_QCcontig25_split_00001</td><td>15856</td><td>3148</td><td>0.0726466575716235</td><td>GT</td><td>G</td><td>0</td><td>213</td><td>3</td><td>2932</td><td>0</td></tr>
        <tr><td>77430</td><td>DAY_24</td><td>Day17a_QCcontig25_split_00001</td><td>835</td><td>3046</td><td>0.831528279181709</td><td>CT</td><td>C</td><td>1</td><td>1382</td><td>1662</td><td>1</td><td>0</td></tr>
        <tr><td>29463</td><td>DAY_15B</td><td>Day17a_QCcontig25_split_00001</td><td>808</td><td>3022</td><td>0.729621125143513</td><td>TG</td><td>T</td><td>5</td><td>1742</td><td>4</td><td>1271</td><td>0</td></tr>
        <tr><td>77442</td><td>DAY_24</td><td>Day17a_QCcontig25_split_00001</td><td>1145</td><td>3021</td><td>0.661529994496423</td><td>AT</td><td>A</td><td>1817</td><td>1202</td><td>2</td><td>0</td><td>0</td></tr>
        <tr><td>77435</td><td>DAY_24</td><td>Day17a_QCcontig25_split_00001</td><td>808</td><td>3014</td><td>0.838827838827839</td><td>TG</td><td>T</td><td>1</td><td>1638</td><td>1</td><td>1374</td><td>0</td></tr>
        <tr><td>29467</td><td>DAY_15B</td><td>Day17a_QCcontig25_split_00001</td><td>577</td><td>3009</td><td>0.959424083769634</td><td>GC</td><td>G</td><td>10</td><td>5</td><td>1466</td><td>1528</td><td>0</td></tr>
     </tbody>
</table>

And this is the command line that generated this output (by the way, you can download the PROFILE.db file used in this example [from here](http://figshare.com/articles/Sharon_et_al_metagenome/1499236)):

{% highlight bash %}
    sqlite3  PROFILE.db 'select * from variable_positions;' | sed 's/|/      /g' | sort -nrk 5 | head -n 30
{% endhighlight %}

## Variable positions table

Essentially, each line in this table simply reports the base frequencies and coverage information of a single nucleotide position, a single sample maps to. A given line of information has no idea about how many other lines there are coming from the same sample, or how many samples there are in total. So, so far, it is only 'reporting'.

But the reporting itself is a little tricky: If you report base frequencies for each position when there is any variation at all, then the number of reported positions would be enormous. For instance, if you have a nucleotie position with 500X coverage, you can safely assume that there will always be _some_ nucleotides that do not match to the consensus nucleotide for that position --whether it is due to biology, or random sequencing errors, or other types of relevant or irrelevant sources of variation. For the sake of better management of available resources, and to be very quick, anvi'o (with its default settings) ___does not___ report variation in every single nucleotie position during profiling.

Instead, it relies on the following conservative heuristic to determine whether to report the variation at a nucleotide position:
 
<div class="centerimg">
<a href="https://www.desmos.com/calculator/qwocua4zi5"><img src="{{ site.url }}/images/anvio/2015-07-20-analyzing-variability/function.png" style="border: none; width: 200px;" /></a>
</div>

where, `x` represents the coverage, and `b`, `m`, and `c` represent the model parameters equal to 3, 1.45, and 0.05, respectively. Assuming `n1` and `n2` represent the frequency of the most frequent and the second most frequent bases in a given nucleotide position (see the table), base frequencies are reported only if `n2/n1 > y` criterion is satisfied for a given coverage of `x`. It is this simple (and ugly, in a sense). But briefly, this approach sets a dynamic baseline for minimum variation required for reporting _as a function of coverage depth_. According to this heuristic, `y` would be 0.29 for 20X coverage (`x`), 0.13 for 50X coverage, 0.08 for 100X coverage, and ~0.05 for very large values of coverage as `y` approaches to `c`. The goal here is to lessen the impact of sequencing and mapping errors in reported frequencies, and it does it's boring job.

This computation- and storage-efficient strategy reports a rather short list of sample-specific variable nucleotide positions that one can really trust. However, the user has always the option to instruct the profiler to store _all observed frequencies_ by declaring this intention with `--report-variability-full` flag for more statistically appropriate downstream analyses. I am of course talking about the type of analyses Christopher Quince-like people would like to do.

And that's that for "reporting". The next is what we do with what is reported.


# Profiling variability

Say we profiled our samples one by one, and merged a bunch of them. Result is a nicely populated `variable_positions` table in the merged database that contains a list ofunconnected nucleotide positions from every sample. The question is how to make sense of them.

To interpret the ecological significance of sample-specific variable positions across samples, anvi’o installs a helper program called `anvi-gen-variability-profile` (AGVP). This program allows the researcher to specify _filters_ that take the experimental design into consideration, so they can generate a refined variability profile that can answer more specific questions.

A very simple command line for AGVP looks like this:

{% highlight bash %}
    anvi-gen-variability-profile -a ANNOTATION.db -p PROFILE.db -c MY_COLLECTION -b MY_BIN
{% endhighlight %}

If you type `anvi-gen-variability-profile --help` in your terminal you will see that AGVP can process variable positions in a genome bin based on multiple user-defined, optional filters, including the max number of variable positions to use from each split (`-n`; default: 0, which means "use everything"), minimum ratio of the competing nucleotides at a reported variable position (`-r`; default: 0), minimum number of samples in which a nucleotide position is reported as a variable position (`-x`; default: 1), and the minimum _scattering power_ of a variable nucleotide position across samples (`-m`; default: 0).

'Scattering power' needs a bit more explanation. This property belongs to a nucleotide position, and it can give a brief idea about the relevance of a given nucleotide position to infer the ecology. OK. Here is an attempt to explain it further: samples in a merged profile can be organized into one or more groups (`g`) based on the nucleotide identity of the competing bases (`b`) at a given variable position, `p`. Scattering power then represents _the number of samples in the second largest group_. That's it. For instance, at one extreme `b` would be identical in all samples at position `p`, in which case `g` would be 1, and scattering power of `p` would be 0. At the other extreme, `p` would harbor a different `b` in every sample, in which case `g` would equal to the number of samples, and the scattering power of `p` would equal to 1. A value for `g` between these two extremes would yield a scattering power of `>1`. The user can employ scattering power to query only those nucleotide positions that vary _consistently_ across samples, and discard positions that show stochastic behavior that are more likely to result from sequencing or PCR errors, or mapping inconsistencies that are not of interest. Please let me know through private e-mail or public means if this is not clear. I would like to improve it if I can identify which part seems to be confusing.

So we have AGVP to _sample_ a shorter list of variable positions, and start using them.

# An example with the *E. facealis* population in an infant's gut

So, what can we learn by focusing on subtle nucleotide variation?

Here I will briefly go through one of the interesting things we have encountered while testing out our approach: the curious case of *E. faecalis*.

You can read the full story on this in our paper. But here is a *ver* brief background: We re-analyzed the metagenomic data [Sharon *et al.* publihsed with a very meticulous analysis](http://genome.cshlp.org/content/23/1/111.long). Agreeing with their findings, *E. facealis* was the most covered genome bin throughout the sampling period in our re-analysis as well. Here is the Figure 2 panel A from the paper; *E. facealis* bin is the red one around 12 o'clock:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2015-07-20-analyzing-variability/infant-gut.png"><img src="{{ site.url }}/images/anvio/2015-07-20-analyzing-variability/infant-gut.png" style="width: 50%;" /></a>
</div>

So we have a genome bin, that appears to be very abundant in every sample. This is the coverage of each split in this genome bin across samples:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2015-07-20-analyzing-variability/coverage.png"><img src="{{ site.url }}/images/anvio/2015-07-20-analyzing-variability/coverage.png" style="width: 50%;" /></a>
</div>

Clearly, we have no way of knowing the extent of variation through this perspective. But for our paper, we wrote a bunch of scripts to analyze the variation in this genome bin. This quite reproducible analysis is stored [in this Github repository](https://github.com/meren/anvio-methods-paper-analyses/tree/master/SHARON_et_al/VARIABILITY_REPORTS). Please take a look at the shell script that runs the same analysis for three genome bins.

If we were to do it just for *E. faecalis* bin, the command line for AGVP would have looked like this:

{% highlight bash %}
anvi-gen-variability-profile -p PROFILE.db \
                             -a ANNOTATION.db \
                             -c SUPERVISED \
                             -b E_faecalis \
                             -n 5 \
                             -o PROFILE_E_faecalis.txt \
                             -m 3 \
                             --quince
{% endhighlight %}

This command requests AGVP to report base frequencies at variable nucleotide positions that have a minimum scattering power of 3 (`-m 3`), for E_faecalis. To do a more even sampling, we also request only up to 5 nucleotide positions from each split (`-m 3`).

So these filters create an even more scrutinized list of nucleotide positions. The output file is very similar to the structure of `variable_positions` table. To visualize what is reported, [we have an R script](https://github.com/meren/anvio-methods-paper-analyses/blob/master/SHARON_et_al/VARIABILITY_REPORTS/02_GEN_FIGURE_SUMMARY.R), which gives us the non-polished version of the following figure for the *E. faecalis* bin in this dataset:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2015-07-20-analyzing-variability/E_faecalis.png"><img src="{{ site.url }}/images/anvio/2015-07-20-analyzing-variability/E_faecalis.png" style="width: 50%;" /></a>
</div>

Through this figure, we now have a more complete picture of whatis going on within this bin, and it is quite interesting. First, slightly edited version of the figure caption from the paper:

<div class="quotable">
The figure displays for the E. faecalis bin in each sample (from top to bottom), (1) average coverage values for all splits, (2) variation density (number of variable positions reported during the profiling step per kilo base pairs), (3) heatmap of variable nucleotide positions, (4) ratio of variable nucleotide identities, and finally (5) the ratio of transitions (mutations that occur from A to G, or T to C, and vice versa) versus transversions. In the heatmap, each row represents a unique variable nucleotide position, where the color of each tile represents the nucleotide identity, and the shade of each tile represents the square root-normalized ratio of the most frequent two bases at that position (i.e., the more variation in a nucleotide position, the less pale the tile is).
</div>


The figure shows that the variation density changes quite dramatically from one day to another, despite the rather stable coverage. Plus, the emergence of this pattern is not random: the heatmap shows that the nucleotide positions that show variation re-occur, and competing base identities remains the same.

Investigating what causes this, is of course when things start to get exciting. However, I will not go there. Instead, I would like to leave you with this thought: By using patterns of variability, we can start characterizing changes in microbial population structures across environments, and generate better-informed hypotheses to investigate mechanisms that drive these shifts.

---
layout: post
title: "How many bacterial genomes you have in that assembly?"
excerpt: "A quick way to get an insight into the number of genomes your contigs represent."
modified: 2015-12-07
tags: []
categories: [anvio]
comments: true
authors: [meren]
---

{% include _toc.html %}

Your ability to identify genomes appropriately from an assembly will depend on many factors, such as the number of samples you have to exploit the differential coverage patterns of genomes, or the algorithm you use to tease apart that information. But even before doing any of the real binning, you may have a rough answer to this question by just taking a quick look from your contigs:

<blockquote>
How many bacterial genomes should I <em>expect</em> to find in these contigs I assembled from my metagenome?

<div class="blockquote-author">Your Name</div>
</blockquote>

Anvi'o has been very practical for us to answer this, and here you will find the workflow to try on your own dataset.

## The workflow

I will describe the three anvi'o steps using [Sharon et al.'s infant gut](http://www.ncbi.nlm.nih.gov/pubmed/22936250) study as an example. The `contigs.fa` file I will use is this section is the co-assembly of all samples in Sharon et al's study.

As a reminder, Sharon et al. identified 12 bacterial draft genomes in this dataset, **8 of which were complete or near-complete**. Our re-analysis of this dataset in the [anvi'o methods paper](https://peerj.com/articles/1319/) also yielded [similar resulted](http://anvio.org/data/INFANT-CLC-SUMMARY-SUPERVISED/) (you can find more about our re-analysis [here]({{ site.url }}/data/#daily-infant-gut-samples-by-sharon-et-al)).

### Generating a contigs database

First, you will need to introduce anvi'o to your contigs by generating a contigs database --one of the essential files of the anvi'o workflow.

{% highlight bash %}
$ anvi-gen-contigs-database -f contigs.fa -o contigs.db
{% endhighlight %}

And this is how this goes on my screen:

{% highlight bash %}
$ anvi-gen-contigs-database -f contigs.fa -o contigs.db
Contigs database .............................: A new database, new-contigs.db, has been created.
Number of contigs ............................: 4,189
Number of splits .............................: 4,861
Total number of nucleotides ..................: 35,766,167
Split length .................................: 20,000

Finding ORFs in contigs
===============================================
Genes ........................................: /var/folders/nm/jmps9l7j7w1dr8rbpv119zgr0000gn/T/tmpZN2Fwt/contigs.genes
Proteins .....................................: /var/folders/nm/jmps9l7j7w1dr8rbpv119zgr0000gn/T/tmpZN2Fwt/contigs.proteins
Log file .....................................: /var/folders/nm/jmps9l7j7w1dr8rbpv119zgr0000gn/T/tmpZN2Fwt/00_log.txt
Result .......................................: Prodigal (v2.6.2) has identified 32265 genes.

Contigs with at least one gene call ..........: 4174 of 4189 (99.6%)

{% endhighlight %}

### Looking for single-copy genes

Now we have the contigs database, the next thing we will do is to look for bacterial single copy genes. If you have read our article you already know that anvi'o installs with four single-copy gene collections from four different groups.

All you need to do is to run this command for HMM hits for thoese gene collections to be added to the database you just created:

{% highlight bash %}
$ anvi-populate-search-table -c contigs.db
{% endhighlight %}

And this is how this step goes on my screen as anvi'o goes through each single-copy gene collection:

{% highlight bash %}
$ anvi-populate-search-table -c contigs.db
HMM profiles .................................: 4 sources have been loaded: Alneberg_et_al (34 genes), Dupont_et_al (111 genes), Campbell_et_al (139 genes), Creevey_et_al (40 genes)
Sequences ....................................: 32265 sequences reported.
FASTA ........................................: /var/folders/nm/jmps9l7j7w1dr8rbpv119zgr0000gn/T/tmpU2EHjX/sequences.fa

HMM Profiling for Alneberg_et_al
===============================================
Reference ....................................: Alneberg et al, http://www.nature.com/nmeth/journal/v11/n11/full/nmeth.3103.html
Pfam model ...................................: /Users/meren/Desktop/MBL/anvio/anvio/data/hmm/Alneberg_et_al/genes.hmm.gz
Number of genes ..............................: 34
Temporary work dir ...........................: /var/folders/nm/jmps9l7j7w1dr8rbpv119zgr0000gn/T/tmpq6RiH1
HMM scan output ..............................: /var/folders/nm/jmps9l7j7w1dr8rbpv119zgr0000gn/T/tmpq6RiH1/hmm.output
HMM scan hits ................................: /var/folders/nm/jmps9l7j7w1dr8rbpv119zgr0000gn/T/tmpq6RiH1/hmm.hits
Log file .....................................: /var/folders/nm/jmps9l7j7w1dr8rbpv119zgr0000gn/T/tmpq6RiH1/00_log.txt
Number of raw hits ...........................: 368

HMM Profiling for Dupont_et_al
===============================================
Reference ....................................: Dupont et al, http://www.nature.com/ismej/journal/v6/n6/full/ismej2011189a.html
Pfam model ...................................: /Users/meren/Desktop/MBL/anvio/anvio/data/hmm/Dupont_et_al/genes.hmm.gz
Number of genes ..............................: 111
Temporary work dir ...........................: /var/folders/nm/jmps9l7j7w1dr8rbpv119zgr0000gn/T/tmpRPCMEL
HMM scan output ..............................: /var/folders/nm/jmps9l7j7w1dr8rbpv119zgr0000gn/T/tmpRPCMEL/hmm.output
HMM scan hits ................................: /var/folders/nm/jmps9l7j7w1dr8rbpv119zgr0000gn/T/tmpRPCMEL/hmm.hits
Log file .....................................: /var/folders/nm/jmps9l7j7w1dr8rbpv119zgr0000gn/T/tmpRPCMEL/00_log.txt
Number of raw hits ...........................: 886

HMM Profiling for Campbell_et_al
===============================================
Reference ....................................: Campbell et al, http://www.pnas.org/content/110/14/5540.short
Pfam model ...................................: /Users/meren/Desktop/MBL/anvio/anvio/data/hmm/Campbell_et_al/genes.hmm.gz
Number of genes ..............................: 139
Temporary work dir ...........................: /var/folders/nm/jmps9l7j7w1dr8rbpv119zgr0000gn/T/tmp9XgXUp
HMM scan output ..............................: /var/folders/nm/jmps9l7j7w1dr8rbpv119zgr0000gn/T/tmp9XgXUp/hmm.output
HMM scan hits ................................: /var/folders/nm/jmps9l7j7w1dr8rbpv119zgr0000gn/T/tmp9XgXUp/hmm.hits
Log file .....................................: /var/folders/nm/jmps9l7j7w1dr8rbpv119zgr0000gn/T/tmp9XgXUp/00_log.txt
Number of raw hits ...........................: 1,361

HMM Profiling for Creevey_et_al
===============================================
Reference ....................................: Creevey et al, http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0022099
Pfam model ...................................: /Users/meren/Desktop/MBL/anvio/anvio/data/hmm/Creevey_et_al/genes.hmm.gz
Number of genes ..............................: 40
Temporary work dir ...........................: /var/folders/nm/jmps9l7j7w1dr8rbpv119zgr0000gn/T/tmp2GoYlD
HMM scan output ..............................: /var/folders/nm/jmps9l7j7w1dr8rbpv119zgr0000gn/T/tmp2GoYlD/hmm.output
HMM scan hits ................................: /var/folders/nm/jmps9l7j7w1dr8rbpv119zgr0000gn/T/tmp2GoYlD/hmm.hits
Log file .....................................: /var/folders/nm/jmps9l7j7w1dr8rbpv119zgr0000gn/T/tmp2GoYlD/00_log.txt
Number of raw hits ...........................: 432
{% endhighlight %}

So far so good.

### Visualizing the results

Now we have a contigs database with everything we need. It is time to visualize the results. For now this is a two-step process. First we need to generate essential input files for the R program that will do the visualization:

{% highlight bash %}
$ anvi-script-gen_stats_for_single_copy_genes.py contigs.db
{% endhighlight %}

This should generate two new files in the directory:

{% highlight bash %}
$ ls
contigs.db  contigs.db.genes  contigs.db.genes
{% endhighlight %}

And the final step is to visualize the information reported in those files:

{% highlight bash %}
$ anvi-script-gen_stats_for_single_copy_genes.R contigs.db.hits contigs.db.genes
[1] "Alneberg_et_al"
[1] "Creevey_et_al"
[1] "Campbell_et_al"
[1] "Dupont_et_al"
{% endhighlight %}

Which should generate a PDF file in the same directory if you have a proper R installation:

{% highlight bash %}
$ ls
contigs.db  contigs.db.genes  contigs.db.genes  contigs.db.hits_e_1_new.pdf
{% endhighlight %}

Here is the result:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2015-12-07-predicting-number-of-genomes/infant-gut.png"><img src="{{ site.url }}/images/anvio/2015-12-07-predicting-number-of-genomes/infant-gut.png" /></a>
</div>

Which suggests that there are 8 to 10 genomes in this assembly!


## Other examples

Here is the result of the worfklow I described above on an assembly generated from the shotgun seqeuncing of a cultivar:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2015-12-07-predicting-number-of-genomes/b-frag-cultivar.png"><img src="{{ site.url }}/images/anvio/2015-12-07-predicting-number-of-genomes/b-frag-cultivar.png" width="70%" /></a>
</div>

|Source|Number of expected bacterial genomes|
|:-----|:---:|
|Alneberg et al.|**1**|
|Creevey et al.|**1**|
|Campbell et al.|**1**|
|Dupont et al.|**1**|

In contrast an ocean sample we assembled recently:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2015-12-07-predicting-number-of-genomes/ocean.png"><img src="{{ site.url }}/images/anvio/2015-12-07-predicting-number-of-genomes/ocean.png" width="70%" /></a>
</div>

|Source|Number of expected bacterial genomes|
|:-----|:---:|
|Alneberg et al.|**451**|
|Creevey et al.|**451**|
|Campbell et al.|**431**|
|Dupont et al.|**354**|

---

Clearly these estimations are mere approximations at best, and should be taken with a grain of salt. However, they *do* give a rough idea about the complexity of a given metagenome, and what you should expect to get from it.

Please don't hesitate to use the comments section for your suggestions, and please keep an eye on our [repository](http://github.com/meren/anvio).

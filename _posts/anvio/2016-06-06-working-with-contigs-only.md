---
layout: post
title: "Binning without mapping"
excerpt: "So you have an assembly, or a draft genome, or a MAG, but no metagenomic short reads? That's OK."
modified: 2016-06-06
categories: [anvio]
comments: true
authors: [meren]
---

The anvi'o [metagenomic workflow]({% post_url anvio/2015-05-02-anvio-tutorial %}){:target="_blank"} assumes that you have metagenomic short reads. But what if all you have is a bunch of contigs, or a draft genome, or a MAG without any mapping data associated with that?

This need was brought up by one of our early users, and there has been an [open issue](https://github.com/meren/anvio/issues/226) to address this at some point. It is now resolved, and the following functionality is available in the [master branch](https://github.com/meren/anvio/tree/master).

The key is to generate a *blank anvi'o profile database* after generating an anvi'o contigs database for your contigs.

Before I start, let me put this here, so you know what version I am using:

{% highlight bash %}
$ anvi-profile -v
Anvi'o version ...............................: 2.0.0
Profile DB version ...........................: 15
Contigs DB version ...........................: 5
Samples information DB version ...............: 2
Auxiliary HDF5 DB version ....................: 1
Users DB version (for anvi-server) ...........: 1
{% endhighlight %}

## Preparing the FASTA file

You probably already have your data ready, but I need an example file for demonstration purposes. So, as an example, I will download a [Bacillus subtilis](http://bacteria.ensembl.org/Bacillus_subtilis_best7613/Info/Index) genome, which is sequenced by Mitsuhiro Itaya *et al.*, and [made available](https://www.ncbi.nlm.nih.gov/bioproject/?term=ASM32874v1) to the community along with [their publication](http://www.pnas.org/content/102/44/15971.long):

"*Combining two genomes in one cell: Stable cloning of the Synechocystis PCC6803 genome in the Bacillus subtilis 168 genome*"

What an awesome title, by the way. OK. The reason I will use this genome is simple: it is a chimeric genome, and it will be fun to see that visually.

Let's first download the genome:

{% highlight bash %}
$ wget ftp://ftp.ensemblgenomes.org/pub/current/bacteria/fasta/bacteria_25_collection/bacillus_subtilis_best7613/dna/Bacillus_subtilis_best7613.ASM32874v1.31.dna.genome.fa.gz -O Bacillus_subtilis.fa.gz
$ gzip -d Bacillus_subtilis.fa.gz
$ ls -lh Bacillus_subtilis.fa
      -rw-r--r-- 7.4M Jun  6 14:42 Bacillus_subtilis.fa
{% endhighlight %}

So far so good. But as you may know, [anvi'o requires simple deflines in FASTA files]({% post_url anvio/2015-05-02-anvio-tutorial %}/#preparation){:target="_blank"}), and this FASTA file *does not* conform that requiremenet at all:

{% highlight bash %}
$ head -n 1 Bacillus_subtilis.fa
>Chromosome dna:chromosome chromosome:ASM32874v1:Chromosome:1:7585470:1
{% endhighlight %}

So I need to fix that first before I generate the contigs database. I will use `anvi-script-remove-short-contigs-from-fasta` with a minimum length of 0 and `--simplify` flag (which will not remove anything from the FASTA file, while converting deflines to simple ones):

{% highlight bash %}
$ anvi-script-remove-short-contigs-from-fasta Bacillus_subtilis.fa --min-len 0 --simplify-names -o fixed.fa
Input ........................................: Bacillus_subtilis.fa
Output .......................................: fixed.fa
Minimum length ...............................: 0
Total num contigs ............................: 1
Total num nucleotides ........................: 7,585,470
Contigs removed ..............................: 0 (0.00% of all)
Nucleotides removed ..........................: 0 (0.00% of all)
Deflines simplified ..........................: True
$ mv fixed.fa Bacillus_subtilis.fa
{% endhighlight %}

Now it looks much better:

{% highlight bash %}
$ head -n 1 Bacillus_subtilis.fa
>c_000000000001
{% endhighlight %}

So this is my FASTA file, and you have yours. We are golden. The next step is to create an anvi'o contigs database.

## Creating the contigs database

Creating a contigs database is not much different for this workflow than any other. But this time (1) I will set a shorter split size than default (since I want to see a more highly resolved depiction of the genome), and (2) I will skip *mindful splitting*, so anvi'o does not lose time with something that is not useful in this context:

{% highlight bash %}
anvi-gen-contigs-database -f Bacillus_subtilis.fa -o Bacillus_subtilis.db -L 5000 --skip-mindful-splitting
{% endhighlight %}

Once this is done, we can populate the single-copy gene hit tables in it:

{% highlight bash %}
anvi-populate-search-table -c Bacillus_subtilis.db
{% endhighlight %}

If you hapenned to read [this post]({% post_url anvio/2015-12-07-predicting-number-of-genomes}) before, you know that at this point we can take a look at the distribuiton of bacterial single-copy genes in this contigs database and predict the number of genomes in it:

{% highlight bash %}
$ anvi-script-gen_stats_for_single_copy_genes.py Bacillus_subtilis.db
$ anvi-script-gen_stats_for_single_copy_genes.R Bacillus_subtilis.db.hits Bacillus_subtilis.db.genes
{% endhighlight %}

These commands gives me back this PDF image:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2016-06-06-working-with-contigs-only/scgs.png"><img src="{{ site.url }}/images/anvio//2016-06-06-working-with-contigs-only/scgs.png" width="80%" /></a>
</div>

From which we confirm that there are two genomes in this contigs database. This is not a case of contamination, but let's assume it was. You did things and figured out that your draft genome is likely contamianted.

Learning that your genome (or your draft genome, or your MAG, name it) is *contaminated* is always very important, but it gives you nothing but a dark afternoon if you don't know how to make sense of that contamination. When there is mapping data available, anvi'o offers powerful ways to deal with contamination as we recently demonstarted in our publication "[Identifying contamination with advanced visualization and analysis practices](https://peerj.com/articles/1839/)". But none of that works if you don't short reads to map. The next step will enable us to bypass that need.

## Creating a blank profile database

To visualize a contigs database, we need an anvi'o profile database. But what does one do if there is no mapping data? Well, they run anvi-profile with `--blank-profile` parameter! That's what one does:

{% highlight bash %}
$ anvi-profile -c Bacillus_subtilis.db -o Bacillus_subtilis --blank-profile
{% endhighlight %}

Here, a man run the anvi-profile to create a blank profile. Although a man has no desires. But thanks to a blank profile database, a man will be able to visualize what a contigs database holds in it using the `anvi-interactive` program. Now a man can do binning with real time completion / contamination estimates, store and load states, create collections, and even summarize them.

## Visualizing

OK. Everything is ready to simply start the anvi'o interactive interface for a better look:

{% highlight bash %}
anvi-interactive -c Bacillus_subtilis.db -p Bacillus_subtilis/PROFILE.db
{% endhighlight %}

Which gives me this view in the opening window:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2016-06-06-working-with-contigs-only/view.png"><img src="{{ site.url }}/images/anvio//2016-06-06-working-with-contigs-only/view.png" width="80%" /></a>
</div>

The inner tree shows the organization of each 5 kbp split in this genome based on their tetra-nucleotide frequency profiles. As you can see, there are two distinct branches in the tree. One of which probably represents the Synechocystis genome, and the other represents the Bacillus subtilis genome Mitsuhiro Itaya *et al.* masterfully put together.

Although this is the simplest look at this data for the sake of clarity, you can indeed use `--additional-view`, or `--additional-layers` to add more data into the interface, and do pretty much anything you would do with any other *regular* anvi'o project.

As an example, here I will separate those genomes from each other by selecting those distinct branches into different bins, and then store my selections into a 'collection' to summarize it later (click to load the animated gif):

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2016-06-06-working-with-contigs-only/selections.gif"><img src="{{ site.url }}/images/anvio/2016-06-06-working-with-contigs-only/selections.png" width="80%" /></a>
</div>

As you can see, completion / redundancy estimates look much better when they are selected separately.

## Summary

Now I can summarize my collection to create a static output with FASTA files for each selection:

{% highlight bash %}
anvi-summarize -c Bacillus_subtilis.db -p Bacillus_subtilis/PROFILE.db -C merens -o Bacillus_subtilis_summary
{% endhighlight %}

And here is the static HTML output this command generates for your viewing pleasure:

[http://anvio.org/data/Bacillus_subtilis_summary/](http://anvio.org/data/Bacillus_subtilis_summary/)

Mission accomplished!

---

All these aside, as of today we managed to address one of the most sorely missed need of every software platfrom:

<blockquote class="twitter-tweet" data-lang="en"><p lang="en" dir="ltr">We heard that <a href="https://twitter.com/hashtag/anvio?src=hash">#anvio</a> will not be taken seriously unless it has a sticker. Well, we got that covered: <a href="https://t.co/24IVCDQqrC">pic.twitter.com/24IVCDQqrC</a></p>&mdash; A. Murat Eren (@merenbey) <a href="https://twitter.com/merenbey/status/739833021055012866">June 6, 2016</a></blockquote>
<script async src="//platform.twitter.com/widgets.js" charset="utf-8"></script>

Just so you know.
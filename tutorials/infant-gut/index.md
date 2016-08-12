---
layout: page
title: An anvi'o tutorial using the Infant Gut Dataset by Sharon <em>et al.</em>
modified: 2016-07-13
excerpt: "Binning, benchmarking binning approaches, refining poorly identified bins, and other Anvi'o tutorial for metagenomics using the Infant Gut Dataset by Sharon *et al*"
comments: true
image:
  feature: feature-infant-gut-tutorial.png
---

{:.notice}
This tutorial is tailored for anvi'o `v2.0.2` or later. You can learn the version of your installation by typing `anvi-interactive -v`. If you have an older version, some things will not work the way they should.

The purpose of this tutorial is to have a conversation about metagenomic binning (while demonstrating some of the anvi'o capabilities) using the Infant Gut Dataset (IGD), which was generated, analyzed, and published by [Sharon et al. (2013)](http://www.ncbi.nlm.nih.gov/pubmed/22936250), and was re-analyzed in the [anvi'o methods paper](https://peerj.com/articles/1319/).

{:.notice}
While this tutorial will take you through a simple analysis of a real dataset, there also is available a more comprehensive (but more abstract tutorial) on [anvi'o metagenomic workflow]({% post_url anvio/2016-06-22-anvio-tutorial-v2 %}).

{% include _toc.html %}

A typical anvi'o metagenomic workflow [starts with BAM files and a FASTA file]({% post_url anvio/2016-06-22-anvio-tutorial-v2 %}#preparation) of contigs. There are many ways to get your contigs and BAM files for your metagenomes. But we have started implementing a tutorial that describes the workflow _we_ use to generate these files regularly: "[A tutorial on assembly-based metagenomics]({{ site.url }}/tutorials/assembly-based-metagenomics/)". Please feel free to take a look at that one, as well.

## The contigs database & anvi'o merged profile for the Infant Gut Dataset

If you are following this tutorial, you will need the anvi'o merged profile database and the anvi'o contigs database for the IGD. You can download them using this link: [https://figshare.com/articles/Infant_Gut_Data_v2/3502445](https://figshare.com/articles/Infant_Gut_Data_v2/3502445) (click the 'Download all', and unzip the resulting directory, and go into it using your terminal).

{:.notice}
Some crash course on anvi'o databases: An **anvi’o contigs database** keeps all the information related to your contigs: positions of open reading frames, k-mer frequencies for each contig, where splits start and end, functional and taxonomic annotation of genes, etc. The contigs database is an essential component of everything related to anvi’o metagenomic workflow. In contrast to the contigs database, an **anvi'o profile database** stores sample-specific information about contigs. Profiling a BAM file with anvi'o creates a single profile that reports properties for each contig in a single sample based on mapping results. Each profile database automatically links to a contigs database, and anvi’o can merge single profiles that link to the same contigs database into **anvi'o merged profile**s (which is what you have in this directory). If you would like to learn more about these, here are some direct links: [creating an anvi'o contigs database]({% post_url anvio/2016-06-22-anvio-tutorial-v2 %}#creating-an-anvio-contigs-database), [creating single anvi'o profiles]({% post_url anvio/2016-06-22-anvio-tutorial-v2 %}#profiling-bam-files), and [merging anvi'o profiles]({% post_url anvio/2016-06-22-anvio-tutorial-v2 %}#working-with-anvio-profiles)

While in that directory, run this command to unpack the additional files directory,

{% highlight bash %}
 $ tar -zxvf additional-files.tar.gz
{% endhighlight %}

{:.notice}
If you were sent here somewhere from down below, now you can **go back**. If you have no idea what this means, ignore this notice, and continue reading. You're OK :)

### Taking a quick look at the merged profile

If you run the following command using `anvi-interactive`,

{% highlight bash %}
 $ anvi-interactive -p PROFILE.db -c CONTIGS.db 
{% endhighlight %}

The anvi'o interactive interface should pop up with this display:

[![Infant gut merged](images/infant-gut-merged.png)](images/infant-gut-merged.png){:.center-img .width-50}

Once the interactive interface is up and running, you can start binning:

[![Infant gut merged](images/infant-gut-merged.gif)](images/infant-gut-merged.gif){:.center-img .width-50}


When you are tired of the interactive interface, you can go back to the terminal and press `CTRL + C` to kill the server.

### Importing taxonomy

Centrifuge ([code](https://github.com/infphilo/centrifuge), [pre-print](http://biorxiv.org/content/early/2016/05/25/054965.article-info)) is [one of the options]({% post_url anvio/2016-06-18-importing-taxonomy %}#centrifuge-output) to [import taxonomic annotations]({% post_url anvio/2016-06-18-importing-taxonomy %}) into an anvi'o contigs database. Centrifuge files for the IGD are already in the directory `additional-files/centrifuge-files`.

If you import these files into the contigs database the following way,

{% highlight bash %}
 $ anvi-import-taxonomy -c CONTIGS.db -i additional-files/centrifuge-files/*.tsv -p centrifuge 
{% endhighlight %}

And run the interactive interface again,

{% highlight bash %}
 $ anvi-interactive -p PROFILE.db -c CONTIGS.db 
{% endhighlight %}

You will see an additional layer with taxonomy:

{:.notice}
In the Layers tab find the `Taxonomy` layer, set its height to `200`, and click `Draw` again. Then click `Save State` button, and overwrite the `default` state. This will make sure anvi'o remembers to make the height of that layer 200px the next time you run the interactive interface!

[![Infant gut merged](images/infant-gut-with-tax.png)](images/infant-gut-with-tax.png){:.center-img .width-50}


## Binning

### External binning results

The directory `additional-files/external-binning-results` contains a number of files that describe the binning of contigs in the IGD based on various binning approaches. These files are `CONCOCT.txt`, `GROOPM.txt`, `MAXBIN.txt`, `METABAT.txt`, `SHARON_et_al.txt`.

The first three files are courtesy of Elaina Graham, who used [GroopM](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4183954/) (v0.3.5), [MetaBat](https://peerj.com/articles/1165/) (v0.26.3), and [MaxBin](https://microbiomejournal.biomedcentral.com/articles/10.1186/2049-2618-2-26) (v2.1.1) to bin the IGD. For future references, here are the parameters Elaina used for each approach:

{% highlight bash %}
 # GroopM v0.3.5 (followed the general workflow on their manual)
 $ groopm parse groopm.db contigs.fa [list of bam files]
 $ groopm core groopm.db -c 1000 -s 10 -b 1000000
 $ groopm recruit groopm.db -c 500 -s 200
 $ groopm extract groopm.db contigs.fa 

 # MetaBat v0.26.3 (used jgi_summarize_bam_contig_depths to get a depth file from BAM files).
 $ metabat -i contigs.fa -a depth.txt -o bin

 # MaxBin v2.1.1
 $ run_MaxBin.pl -contig contigs.fa -out maxbin_IGM -abund_list [list of all coverage files in associated format]
{% endhighlight %}

[CONCOCT](http://www.nature.com/nmeth/journal/v11/n11/full/nmeth.3103.html) results come from the CONCOCT module embedded within anvi'o. 

Finally, I created a file corresponding to Sharon et al. results by BLAST-searching sequences in bins identified by the authors of the study (see [http://ggkbase.berkeley.edu/carrol](http://ggkbase.berkeley.edu/carrol)) to our contigs to have matching names for our assembly.

Now you have the background information about where these files are coming. Moving on.

### Importing an external binning result

In the anvi'o lingo, a 'collection' is something that describes one or more bins, each of which describe one or more contigs. You can create a collection by using the interactive interface, or you can import external binning results into your profile database as a collection and see how that collection groups contigs. For instance, let's import the CONCOCT collection:

{% highlight bash %}
 $ anvi-import-collection -c CONTIGS.db -p PROFILE.db -C CONCOCT --contigs-mode additional-files/external-binning-results/CONCOCT.txt
{% endhighlight %}

And run the interactive again,

{% highlight bash %}
 $ anvi-interactive -p PROFILE.db -c CONTIGS.db 
{% endhighlight %}

When you click `Bins > Load bin collection > CONCOCT > Load` once the interface is ready, this is what you will see:

[![Infant gut merged](images/infant-gut-concoct.png)](images/infant-gut-concoct.png){:.center-img .width-50}

So far so good.

### Comparing multiple binning approaches

Since we have all these results from different binning approaches, it clearly would have been interesting to compare them to each other (because [benchmarking stuff]({% post_url anvio/2015-06-23-comparing-different-mapping-software %}) is often very insightful). But how to do it? The simplest way to do it is to assume a 'true organization of contigs', and then investigate every other approach with respect to that.

Here we have the organization of contigs based on hierarchical clustering analysis, taxonomy from Centrifuge per contig (which is independent of the organization of contigs so it is a very good validation to see whether the organization makes sense), and results from the original publication from Sharon et al., in which authors did a very careful job to identify every genome in the dataset, including resolving the *Staphylococcus* pangenome (which is extremely hard for automatic binning approaches to resolve with one co-assembly). So these are the things we will assume "true enough" to build upon.

To compare binning results, we could import each collection into the profile database the way we imported CONCOCT. But unfortunately at any given time there could only be one collection that can be displayed in the interface. Luckily there are other things we can do. For instance, as a workaround, we can merge all binning results into a single file, and use that file as an 'additional data file' to visualize them in the interactive interface.

Anvi'o has a script to merge multiple files for external binning results into a single merged file (don't ask why):

{% highlight bash %}
 $ anvi-script-merge-collections -c CONTIGS.db -i additional-files/external-binning-results/*.txt -o collections.tsv
{% endhighlight %}

If you take a look at this file, you will realize that it has a very simple format:

{% highlight bash %}
 $ head collections.tsv | column -t
contig                           CONCOCT  GROOPM     MAXBIN      METABAT               SHARON_et_al
Day17a_QCcontig1000_split_00001  Bin_4    db_bin_11  maxbin.008  metabat_igm.unbinned  Finegoldia_magna
Day17a_QCcontig1001_split_00001  Bin_7    db_bin_46  maxbin.006  metabat_igm.unbinned  Staphylococcus_epidermidis_virus_014
Day17a_QCcontig1002_split_00001  Bin_4    db_bin_11  maxbin.007  metabat_igm.unbinned  Finegoldia_magna
Day17a_QCcontig1003_split_00001  Bin_2    db_bin_1   maxbin.009  metabat_igm.7
Day17a_QCcontig1004_split_00001  Bin_3    db_bin_8   maxbin.008  metabat_igm.10
Day17a_QCcontig1005_split_00001  Bin_5    db_bin_47  maxbin.007  metabat_igm.unbinned  Staphylococcus_epidermidis_viruses
Day17a_QCcontig1006_split_00001  Bin_3    db_bin_8   maxbin.008  metabat_igm.10        Leuconostoc_citreum
Day17a_QCcontig1007_split_00001  Bin_3    db_bin_8   maxbin.008  metabat_igm.10
Day17a_QCcontig1008_split_00001  Bin_2    db_bin_1   maxbin.009  metabat_igm.7         Candida_albcans
{% endhighlight %}

Good. Now you can run the interactive interface to display all collections of bins stored in `collections.tsv` as 'additional layers':

{% highlight bash %}
 $ anvi-interactive -p PROFILE.db -c CONTIGS.db -A collections.tsv
{% endhighlight %}

At this point you should be seeing a display similar to this (after setting the height of each additional layer to 200px):

[![Infant gut merged](images/infant-gut-collections.png)](images/infant-gut-collections.png){:.center-img .width-50}

To emphasize relationships between bins visually. If you import it the following way, 

{% highlight bash %}
 $ anvi-import-state --state additional-files/state-files/state-merged.json --name default -p PROFILE.db
{% endhighlight %}

and run the interactive interface again,

{% highlight bash %}
 $ anvi-interactive -p PROFILE.db -c CONTIGS.db -A collections.tsv
{% endhighlight %}

this time you should get this display:

[![Infant gut merged](images/infant-gut-collections-final.png)](images/infant-gut-collections-final.png){:.center-img .width-50}

So far so good?

Let's also import the collection we published in the [anvi'o methods paper](https://peerj.com/articles/1319/) with the proper colors:

{% highlight bash %}
 $ anvi-import-collection additional-files/collections/merens.txt \
                          -p PROFILE.db \
                          -c CONTIGS.db \
                          -C merens \
                          --bins-info additional-files/collections/merens-info.txt
{% endhighlight %}

Now you can rerun the interactive interface, and click `Bins > Load bin collection > merens > Load` to display our collection in comparison:

[![Infant gut merged](images/infant-gut-collections-final-w-merens.png)](images/infant-gut-collections-final-w-merens.png){:.center-img .width-50}

Much better!

Now we can discuss about the efficacy of different approaches.

As a reminder, you can in fact investigate the taxonomy of contigs by BLASTing them against NCBI's collection using the right-click menu to have a second opinion about what do public databases think they are:

[![Infant gut merged](images/infant-gut-split.png)](images/infant-gut-split.png){:.center-img .width-50}

### Refining automatically identified bins

OK. Let's assume, we didn't see the interactive interface, and we have no idea about the dataset. We didn't do any of the things we did up to this point. We just had profiled and merged the IGD, and we did binning of this dataset using MaxBin. Let's start by importing MaxBin results into the profile database as a collection:

{% highlight bash %}
 $ anvi-import-collection additional-files/external-binning-results/MAXBIN.txt\
                          -c CONTIGS.db -p PROFILE.db -C MAXBIN --contigs-mode
 					
{% endhighlight %}

From here, there are two things we can do very quickly. First, we can create a summary of our new collection (which would generate a comprehensive static output for all the bins that are identified):

{% highlight bash %}
 $ anvi-summarize -c CONTIGS.db -p PROFILE.db -C MAXBIN -o MAXBIN-SuMMARY
{% endhighlight %}

{:.notice}
You can learn more about `anvi-summarize` [here]({% post_url anvio/2016-06-22-anvio-tutorial-v2 %}#anvi-summarize).

Second, we can take a very quick look at the binning results by initiating the interactive interface in the `collection` mode:

{% highlight bash %}
 $ anvi-interactive -p PROFILE.db -c CONTIGS.db -C MAXBIN
{% endhighlight %}

This command should give you a display similar to this:

[![MaxBin results in collection mode](images/maxbin-collection.png)](images/maxbin-collection.png){:.center-img .width-50}

All previous interactive displays were at the contig-level. However, this display is at the *bin-level*. Instead of contigs, this display shows us the distribution of *bins* MaxBin identified. We also have completion and redundancy estimates for each bin, which helps us make some early sense of what is going on.

{:.notice}
Please read this post to learn more about completion and redundancy estimates: [Assessing completion and contamination of metagenome-assembled genomes]({% post_url anvio/2016-06-09-assessing-completion-and-contamination-of-MAGs %})

It is clear that some bins are not as well-resolved as others. For instance, bins `maxbin.007` and `maxbin.008` have redundancy estimates of 90% and 59%, respectively, which suggests each of them describe multiple groups of organisms.

{:.notice}
"Why not `maxbin.009`, Meren?". Yes, I am cheating. But for a good reason. Does anyone have any idea what that reason might be? Hint: what are bacterial single-copy core genes not good for?

Well, clearly we would have preferred those bins to *behave*.

Luckily anvi'o allows the refinement of bins. So in fact you can focus only those two bins, and refine them using the interactive interface. First we need to crete a file with bin names we are interested in to refine:

{% highlight bash %}
 $ python -c 'print("maxbin.007\nmaxbin.008")' > maxbin-bins-to-refine.txt
 $ cat maxbin-bins-to-refine.txt
maxbin.007
maxbin.008
{% endhighlight %}

{:.notice}
Normally we can give a bin name to `anvi-refine` using the `-b` parameter. But since we are interested in more than one bin name, `-B` parameter with a file with names is the right way to go with.

Good. Now we can call `anvi-refine`:

{% highlight bash %}
 $ anvi-refine -p PROFILE.db -c CONTIGS.db -C MAXBIN -B maxbin-bins-to-refine.txt
{% endhighlight %}

Which would give us this:

[![MaxBin refining two bins](images/maxbin-refining.png)](images/maxbin-refining.png){:.center-img .width-50}

Here we see the distribution of each contig originally binned into `maxbin.007` and `maxbin.008` across samples. The hierarchical clustering did pick up some trends, and you can see there are more than two clusters you can identify somewhat easily. Here I made some selections:

[![MaxBin two bins refined](images/maxbin-refined.png)](images/maxbin-refined.png){:.center-img .width-50}

Once you are satisfied, you can store new selections to the database from the `Bins` tab, and kill the server once you see the confirmation on the top-right part of the window.

If you run the interactive interface for the collection `MAXBIN` again,:

{% highlight bash %}
 $ anvi-interactive -p PROFILE.db -c CONTIGS.db -C MAXBIN 
{% endhighlight %}

things will look much better:

[![MaxBin results in collection mode](images/maxbin-collection-refined.png)](images/maxbin-collection-refined.png){:.center-img .width-50}

As you can see, even if the initial results from an automatic binning approach contain poorly identified bins, it is possible to improve the final results through refinement steps.

{:.notice}
You can read more about `anvi-refine` [here]({% post_url anvio/2015-05-11-anvi-refine %}).

Thank you for following the tutorial this far!

## Meren's two cents on binning

Binning is inherently a very challenging task.

In most cases it is absolutely doable, especially when there is a decent assembly, but it *is* very challenging.

The IGD is one of the most friendly metagenomic datasets available to play with (since an astonishing fraction of nucleotides map back to the assembly), and it comes from a well-implemented experimental design (because that's what [Banfield group]((http://geomicrobiology.berkeley.edu/)) does). Yet, you now have seen the extent of disagreement between multiple binning approaches even for this dataset.

You should reming yourself that each of these approaches are implemented by people who are well-trained scientists working with groups of people who are experts in their fields. These tools are benchmarked against others and showed improvements. So each one of them provides the *best result* compared to all others in *at least* one metagenomic dataset. I think understanding what this means is important. There is no silver bullet in the common bioinformatics toolkit that will take care of every dataset when you fire it. In fact, depending on the dataset, even the best tools we have may be as efficient as sticks and stones against the Death Star. Computational people are working very hard to improve things, but they would be the first ones to suggest that their tools should never make users feel free from the fact that it is their own responsibility to make sure the results are meaningful and appropriate.

So which one to choose? How to get out of this situation easily and move on? I know how much desire there is to outsource everything we do to fully automated computational solutions. I also acknowledge that the ability to do that is important to perform large-scale and reproducible analyses without going through too much pain. But we are not at a stage yet with metagenomics where you can rely on any of the available automated binning tools, and expect your MAGs to be safe and sound.

For instance, I think CONCOCT is doing a pretty awesome job identifying MAGs in the IGD, even with the low-abundance organisms. However, it is not perfect, either. In fact if you look carefully, you can see that it creates two bins for one *Candida albicans* genome. Hierarchical clustering will always get you closest to the best organization of contigs with simple distance metrics and linkage algorithms. But there are major challenges associated with that approach, including the fact that it is simply an exploratory method and can't give you "bins" out-of-the-box. Even more importantly, it has tremendous limitations come from its computational complexity (~*O*(*m*<sup>2</sup> log *m*), where *m* is the number of data points). So in most cases it is not even a remote possibility to organize contigs using a hierarchical clustering approach in an assembly in reasonable amount of time (and there is no way to visualize that even if you were to get a dendrogram for 200,000 contigs (you can create simple 2D ordinations with that number of items, but you really shouldn't, but that's another discussion)). Except assemblies with rather smaller number of contigs like the IGD, we are always going to use automated ways to identify bins, at least *initially*, knowing that resulting bins may be, and in most cases will be, crappy. That's why in anvi'o we implemented ways to quickly look into automatically identified bins (i.e., the `collection` mode of `anvi-interactive`), and even refine those with poor redundancy scores to improve final results (i.e., `anvi-refine`).

So we can fix crappy bins to an extent since [we know more or less how things should look like]({% post_url anvio/2016-06-09-assessing-completion-and-contamination-of-MAGs %}), and [we have tools to do that]({% post_url anvio/2015-05-11-anvi-refine %}). That being said, there is one more consideration that is very easy to miss. Although it is somewhat possible to recover from **conflation error** (i.e., more than one genome ends up in one bin), it is much harder to recover from the **fragmentation error** (i.e., one genome is split into multiple bins). You can see an example for fragmentation error if you take a careful look from this figure (i.e., CONCOCT bins between 9:30 and 12:00 o'clock, or MaxBin bins between 5:00 and 7:00 o'clock):

[![Infant gut merged](images/infant-gut-collections-final.png)](images/infant-gut-collections-final.png){:.center-img .width-50}

This is a problem that likely happens quite often, and very hard to deal with once the bins are identified. But we *can* recover from that.

#### From fragmentation to conflation error: A Meren Lab Heuristic to fight back

One of the heuristics we recently started using in our lab to avoid fragmentation error is to confine CONCOCT's clustering space to a much smaller number of clusters than the expected number of bacterial genomes in a given dataset, and then curate resulting contaminated bins manually. Let's say we expect to find `n` bacterial genomes, so we run CONCOCT with a maximum number of clusters of about `n/2` (no judging! I told you it was a heuristic!).

<blockquote>
Well, how do you even know how many bacterial genomes you should expect to find in a metagenome?

<div class="blockquote-author">You</div>
</blockquote>

Thanks for the great question. Although this may sound like a challenging problem to some, we have a very simple way to resolve it (which I described in this [blog post]({% post_url anvio/2015-12-07-predicting-number-of-genomes %})). If you still have access to the IGD, you can run these commands,

{% highlight bash %}
 $ anvi-script-gen_stats_for_single_copy_genes.py CONTIGS.db
 $ anvi-script-gen_stats_for_single_copy_genes.R CONTIGS.db.hits CONTIGS.db.genes
{% endhighlight %}

to get an output PDF file, `contigs.db.hits_e_1_new.pdf`, which would indicate that one should expect 8 to 10 near-complete genomes in this dataset.

Fine. Using `anvi-cluster-with-concoct` program, we ask CONCOCT to naively identify 5 clusters in this dataset, and store the results in the profile database as a collection:

{% highlight bash %}
 $ anvi-cluster-with-concoct -p PROFILE.db -c CONTIGS.db --num-clusters 5 -C CONCOCT_C5
{% endhighlight %}

Now you can run the interface again,

{% highlight bash %}
 $ anvi-interactive -p PROFILE.db -c CONTIGS.db --title 'Infant Gut Time Series by Sharon et al. [w/ 5 CONCOCT clusters]'
{% endhighlight %}

and you would see this after loading the new `CONCOCT_C5` collection from the `Bins` tab:

[![CONCOCT w/ 5 clusters](images/concoct-5-clusters.png)](images/concoct-5-clusters.png){:.center-img .width-50}

As you can see, there aren't fragmentation error anymore, and in fact CONCOCT did an amazing job to identify general patterns in the dataset. Now refining these bins to fix all the conflation errors would be much more easier. If you would like to try, here is an example:

{% highlight bash %}
 $ anvi-refine -p PROFILE.db -c CONTIGS.db -C CONCOCT_C5 -b Bin_1
{% endhighlight %}

---

There are more ways to improve bins and binning results. But although we have seen major improvements in our research by exploring these directions, there are also many other cases nothing is quite enough.

Then it is time to increase the depth of sequencing, implement a different assembly strategy, rethink the sampling strategy, or change the experimental approach to do what seems to be undoable. Here is an example from Tom Delmont et al. to that last point with soil metagenomics: [doi:10.3389/fmicb.2015.00358](http://dx.doi.org/10.3389/fmicb.2015.00358).

We all just have to continue working, and enjoy this revolution.


## Profiling SNVs in a Bin

Here we will profile the single-nucleotide variations (SNVs) in the *E. faecalis* bin found in Sharon et al.'s Infant Gut Dataset (IGD).

{:.notice}
This is more of a practical tutorial for hands on experience to recover and make sense of SNVs. For a more theoretical one on the same topic, please consider first reading the tutorial [Analyzing single nucleotide variations (SNVs) with anvi'o]({% post_url anvio/2015-07-20-analyzing-variability %}).

{:.notice}
**If you haven't followed the previous sections of the tutorial**, you will need the anvi'o merged profile database and the anvi'o contigs database for the IGD available to you. Before you continue, please [click here](#the-contigs-database--anvio-merged-profile-for-the-infant-gut-dataset), do everything mentioned there, and come back here on this page when you are asked to **go back**.

Please run following commands in the IGD dir. They will set the stage for us to take a look at the *E. faecalis* bin:

{% highlight bash %}
 $ anvi-import-taxonomy -c CONTIGS.db -i additional-files/centrifuge-files/*.tsv -p centrifuge
 $ anvi-import-state --state additional-files/state-files/state-merged.json --name default -p PROFILE.db
 $ anvi-import-collection additional-files/collections/e-faecalis.txt \
                          -p PROFILE.db \
                          -c CONTIGS.db \
                          -C E_faecalis \
                          --bins-info additional-files/collections/e-faecalis-info.txt
{% endhighlight %}

OK. Let's first see where this *E. faecalis* bin again. If you now run the interactive interface the following way,

{% highlight bash %}
 $ anvi-interactive -p PROFILE.db -c CONTIGS.db
{% endhighlight %}

And then click `Bins > Load bin collection > E_faecalis > Load`, you will see this:

[![E. faecalis bin](images/e-faecalis-bin.png)](images/e-faecalis-bin.png){:.center-img .width-50}

The red selection in the most outer layer represents the *E. faecalis* bin, which is very abundant in every sample.

Here is a different representation of the coverage of this bin across samples (which is a screenshot from our paper):

<div class='centerimg'>
<a href='{{ site.url }}/images/anvio/2015-07-20-analyzing-variability/coverage.png'><img src='{{ site.url }}/images/anvio/2015-07-20-analyzing-variability/coverage.png' style='width: 50%;' /></a>
</div>

See, it *really* is abundant (every dot here is the coverage of a nucleotide position that was reported as a variable position).

OK. Clearly, we have no way of knowing the extent of variation within this bin through this perspective. But `anvi-gen-variability-profile` is exactly for that purpose, and that's what we will do here using two different methods to visualize its report (using R, and using anvi'o).

Let's first generate the SNV profile output file, [details of which were described extensively here]({% post_url anvio/2015-07-20-analyzing-variability %}#the-output-matrix). Here it is:

{% highlight bash %}
$ anvi-gen-variability-profile -c CONTIGs.db \
                              -p PROFILE.db \
                              -C E_faecalis \
                              -b E_faecalis \
                              --samples-of-interest additional-files/samples.txt \
                              --min-coverage-in-each-sample 20 \
                              --include-split-names \
                              --min-scatter 3 \
                              --quince-mode \
                              -o E-faecalis-SNVs.txt
{% endhighlight %}

This command simply requests `anvi-gen-variability-profile` to select all *E. faecalis* **nucleotide positions** that were identified as 'variable' in at least one sample, and (1) covered by more than 20X in **all** 8 samples of interest (2) and display a minimum scattering power of 3 (minimum scattering power is a very simple attribute of a nucleotide position, and described [here]({% post_url anvio/2015-07-20-analyzing-variability %}#parameters-to-filter-the-output)). Instead of `--min-scatter 3` you could use `--min-occurrence 3`, however, in that case the program would have returned every SNV position that was reported in more than 3 samples, which would have included nucleotide positions that were variable in every sample. Finally we also used the `--samples-of-interest` parameter with the following file:

{% highlight bash %}
 $ cat additional-files/samples.txt
DAY_15B
DAY_16
DAY_17B
DAY_18
DAY_19
DAY_22A
DAY_23
DAY_24
{% endhighlight %}

which contains only the sample names with larger number of reads for every duplicate day to simplify things.

OK. Our beloved `anvi-gen-variability-profile` will kindly store results from 466 nucleotide positions in the output file, `E-faecalis-SNVs.txt`. Because we used the flag `--quince-mode`, there will be a total of 3,728 entries in the file (=466 x 8).

At this point we now know that the *E. faecalis* is not a monoclonal bin, and does maintain some level of heterogeneity. What we don't know is whether it is just noise, or if there is any signal that may tell us something about the nature of this variability.

Now we will visualize the information two different ways.

### Visualizing SNV profiles using R

First, we will use R to recapitulate what we did in our [paper](https://peerj.com/articles/1319/). Of course, please feel absolutely free to look at the Figure 3 behind that link if you are fine with RUINING THE SURPRISE :( 

For this step we will need [this R script](https://github.com/meren/anvio-methods-paper-analyses/blob/master/SHARON_et_al/VARIABILITY_REPORTS/02_GEN_FIGURE_SUMMARY.R) Meren had written before. You can download it the following way:

{% highlight bash %}
 $ wget https://raw.githubusercontent.com/meren/anvio-methods-paper-analyses/master/SHARON_et_al/VARIABILITY_REPORTS/02_GEN_FIGURE_SUMMARY.R \
        -O visualize-SNVs.R
 $ chmod +x visualize-SNVs.R
{% endhighlight %}

Here is an extra step, and why we need it: `E-faecalis-SNVs.txt` contains filtered SNVs, but we also want to know about each sample's *variation density*, the number of nucleotide positions reported as a variable position for each 1,000 nts. The R script can already do it, but we need another file that reports ALL SNVs for every sample, without any filters to get the raw counts. Let's first create that the following way:

{% highlight bash %}
$ anvi-gen-variability-profile -c CONTIGs.db \
                               -p PROFILE.db \
                               -C E_faecalis \
                               -b E_faecalis \
                               --samples-of-interest additional-files/samples.txt \
                               -o E-faecalis-SNV-density.txt
{% endhighlight %}

Now we can use the R script to visualize this information the following way,

{:.notice}
This R script will require some libraries to be installed. You can install all of them by typing `R` in your terminal, and then entering this command: "install.packages(c('ggplot2', 'reshape2', 'reshape', 'gridExtra', 'grid', 'plyr', 'gtools'))". Finally, my R version that worked with this script was `v3.2.2`.

{% highlight bash %}
$ ./visualize-SNVs.R E-faecalis-SNVs.txt E-faecalis-SNV-density.txt 150 2870000
{% endhighlight %}

where `150` is the number of random positions we want to show in the heatmap, and `2780000` is the genome size to calculate the SNV density per sample. If you have all necessary R libraries installed, you should see a new PDF file in the directory that looks like this:

[![E. faecalis bin](images/e-faecalis-SNVs-R.png)](images/e-faecalis-SNVs-R.png){:.center-img .width-40}

<div class='quotable'>
The figure displays for the E. faecalis bin in each sample (from top to bottom), (1) average coverage values for all splits, (2) variation density (number of variable positions reported during the profiling step per kilo base pairs), (3) heatmap of variable nucleotide positions, (4) ratio of variable nucleotide identities, and finally (5) the ratio of transitions (mutations that occur from A to G, or T to C, and vice versa) versus transversions. In the heatmap, each row represents a unique variable nucleotide position, where the color of each tile represents the nucleotide identity, and the shade of each tile represents the square root-normalized ratio of the most frequent two bases at that position (i.e., the more variation in a nucleotide position, the less pale the tile is).
</div>


With some rearrangement and polishing using Inkscape, you can see how this output fit into our [Figure 3](https://peerj.com/articles/1319/#fig-3):

[![Eren et al., Figure 3](images/eren_et_al_fig_3.png)](images/eren_et_al_fig_3.png){:.center-img .width-70}

The figure shows that the variation density changes quite dramatically from one day to another, despite the rather stable coverage. Plus, the emergence of this pattern is not random: the heatmap shows that the nucleotide positions that show variation re-occur, and competing base identities remains the same.

Investigating what causes this, is of course when things start to get exciting. However, we will not go there. Instead, we would like to leave you with this thought: by using patterns of variability, we can start characterizing changes in microbial population structures across environments, and generate better-informed hypotheses to investigate mechanisms that drive these shifts.


### Visualizing SNV profiles using anvi'o

R visualization is useful, but the heatmap in that figure can't effectively visualize more than a couple hundred positions. That's why there is a random subsampling step. But we can use the anvi'o interactive interface to display up to 25,000 nucleotide positions easily.

For this, I wrote a little program to create an *anvi'o interactive interface-compatible* output. Although this is very preliminary now, probably it will do much more in the future versions of anvi'o depending on your or our needs. Please download the script first:

{% highlight bash %}
 $ wget https://raw.githubusercontent.com/meren/anvio/master/sandbox/anvi-script-snvs-to-interactive
{% endhighlight %}

When you run this script the following way,

{% highlight bash %}
 $ python anvi-script-snvs-to-interactive E-faecalis-SNVs.txt -o e_faecalis_snvs
{% endhighlight %}

it will do its magic, and create an output directory with material that can directly be used with `anvi-interactive` with the `--manual` flag.

{:.notice}
A little note for people who are interested in programming: Feel free to take a look at the [relevant line of the source code](https://github.com/meren/anvio/blob/master/sandbox/anvi-script-snvs-to-interactive#L68) of this script to see how easy it is to generate an anvi'o-compatible visualizable output from any TAB-delimited matrix file.

If you run the interactive interface on these results the following way,

{% highlight bash %}
 $ anvi-interactive -d e_faecalis_snvs/view_data.txt \
                    -s e_faecalis_snvs/samples.db \
                    -t e_faecalis_snvs/tree.txt \
                    -p e_faecalis_snvs/profile.db \
                    -A e_faecalis_snvs/additional_view_data.txt \
                    --title "SNV Profile for the E. faecalis bin" \
                    --manual
{% endhighlight %}

You will get this view:

[![E. faecalis SNVs](images/e-faecalis-SNVs-anvio.png)](images/e-faecalis-SNVs-anvio.png){:.center-img .width-50}

This view can definitely be improved. I prepared a state file to match colors of competing nucleotides to the R results. If you import that state file and run the interactive interface the following way,

{% highlight bash %}
 $ anvi-import-state -p e_faecalis_snvs/profile.db --state additional-files/state-files/state-snvs.json --name default
 $ anvi-interactive -d e_faecalis_snvs/view_data.txt \
                    -s e_faecalis_snvs/samples.db \
                    -t e_faecalis_snvs/tree.txt \
                    -p e_faecalis_snvs/profile.db \
                    -A e_faecalis_snvs/additional_view_data.txt \
                    --title "SNV Profile for the E. faecalis bin" \
                    --manual
{% endhighlight %}

This time you will get this display:

[![E. faecalis SNVs](images/e-faecalis-SNVs-anvio-state.png)](images/e-faecalis-SNVs-anvio-state.png){:.center-img .width-50}

As we've seen before, occurrence of SNVs follow a bi-daily fashion. Not that it needs any further convincing, but just to show off here, if you were to click `Samples > Sample order > view_data > Draw`, you can see that even days and odd days nicely separate from each other:

[![E. faecalis SNVs](images/e-faecalis-SNVs-anvio-state-clustered.png)](images/e-faecalis-SNVs-anvio-state-clustered.png){:.center-img .width-50}

For instance while I was on this page, I selected two SNV positions that showed different patterns:

[![E. faecalis SNVs](images/e-faecalis-SNVs-anvio-state-clustered-selected.png)](images/e-faecalis-SNVs-anvio-state-clustered-selected.png){:.center-img .width-50}

Position 80 (the selection towards the right) shows variability only in odd days, position 686 shows variability only in even days. This will probably at some point will become a right-click menu function, but even now we can do some tricks to explore the context of these SNVs. Because that is what anvi'o is all about. Exploring, sometimes to disturbing depths. Well, here are two tiny AWK one-liner to get some specific information about these positions from our files. This is for the one on the left:

{% highlight bash %}
 $ awk '{if(NR == 1 || $2 == 686) print $2 " " $3 " " $4 " " $6 " " $14 " " $15 " " $25}' E-faecalis-SNVs.txt | column -t
unique_pos  sample_id  pos   gene_call  departure_from_ref  competing_nts  split_name
686         DAY_18     6333  1826       0.0617647058824     CT             Day17a_QCcontig4_split_00022
686         DAY_22A    6333  1826       0.0697674418605     CT             Day17a_QCcontig4_split_00022
686         DAY_16     6333  1826       0.0805369127517     CT             Day17a_QCcontig4_split_00022
686         DAY_24     6333  1826       0.0863636363636     CT             Day17a_QCcontig4_split_00022
686         DAY_17B    6333  1826       0                   CC             Day17a_QCcontig4_split_00022
686         DAY_19     6333  1826       0                   CC             Day17a_QCcontig4_split_00022
686         DAY_15B    6333  1826       0                   CC             Day17a_QCcontig4_split_00022
686         DAY_23     6333  1826       0                   CC             Day17a_QCcontig4_split_00022
{% endhighlight %}

And this is for the one on the right:

{% highlight bash %}
 $ awk '{if(NR == 1 || $2 == 80) print $2 " " $3 " " $4 " " $6 " " $14 " " $15 " " $25}' E-faecalis-SNVs.txt | column -t
unique_pos  sample_id  pos   gene_call  departure_from_ref  competing_nts  split_name
80          DAY_17B    7955  233        0.122591943958      GT             Day17a_QCcontig1_split_00012
80          DAY_19     7955  233        0.0752688172043     GT             Day17a_QCcontig1_split_00012
80          DAY_15B    7955  233        0.109271523179      GT             Day17a_QCcontig1_split_00012
80          DAY_23     7955  233        0.11377245509       GT             Day17a_QCcontig1_split_00012
80          DAY_18     7955  233        0                   TT             Day17a_QCcontig1_split_00012
80          DAY_22A    7955  233        0                   TT             Day17a_QCcontig1_split_00012
80          DAY_16     7955  233        0                   TT             Day17a_QCcontig1_split_00012
80          DAY_24     7955  233        0                   TT             Day17a_QCcontig1_split_00012
{% endhighlight %}

Good, everything checks out. Now since we know the split names and positions in splits, we can in fact see where they actually are using the interactive interface to visualize the merged profile database again, and look at the wider context using the 'inspect' option. Which I have already done for you:

[![E. faecalis SNVs](images/e-faecalis-SNV-context.gif)](images/e-faecalis-SNV-context.gif){:.center-img .width-50}

There are many directions you can go once you have the gene caller IDs associated with a question you have. Just take a look at this post and see some of the hidden talents of anvi'o: [Musings over a *Nitrospira* genome that can do complete nitrification]({% post_url anvio/2015-12-09-musings-over-commamox %}).

Here I will stop, but still we have a lot to talk about!


### Visualizing SNV profiles as a network

Finally, you can generate an XML description of the SNV profiles you have generated using `anvi-gen-variability-profile` program, using the program `anvi-gen-variability-network`:

anvi-gen-variability-network -i E-faecalis-SNVs.txt -o E-faecalis-SNVs.gexf

You can use [Gephi](https://gephi.org/) to play with the resulting file to visualize or to analyze the network properties of nucleotide positions further.

Here is a screenshot from Gephi for SNV profiles in the E. faecalis genome bin:

[![E. faecalis SNVs network](images/network.png)](images/network.png){:.center-img .width-100}


## A pangenomic analysis

You can also use anvi'o to perform pangenomic analyses, and here I will give a small demo using the now infamous *E. faecalis* bin.

{:.notice}
**If you haven't followed the previous sections of the tutorial**, you will need the anvi'o merged profile database and the anvi'o contigs database for the IGD available to you. Before you continue, please [click here](#the-contigs-database--anvio-merged-profile-for-the-infant-gut-dataset), do everything mentioned there, and come back here on this page when you are asked to **go back**.

Please run following commands in the IGD dir. They will set the stage for us to take a look at the *E. faecalis* bin:

{% highlight bash %}
 $ anvi-import-collection additional-files/collections/e-faecalis.txt \
                          -p PROFILE.db \
                          -c CONTIGS.db \
                          -C E_faecalis \
                          --bins-info additional-files/collections/e-faecalis-info.txt
{% endhighlight %}

---

For this example I downloaded 6 *E. facealis*, and 5 *E. faecium* genomes to analyze them together with our *E. facealis* bin. For each of these 11 *external genomes*, I generated anvi'o contigs databases. You can find all of them in the additional files directory:

{% highlight bash %}
 $ ls -lh additional-files/pangenomics/external-genomes/
total 185984
-rw-r--r--  1 meren  staff   2.6K Jul 31 23:56 00_INFO_ABOUT_EXTERNAL_GENOMES.txt
-rw-r--r--  1 meren  staff   4.7M Jul 31 23:56 Enterococcus_faecalis_6240.db
-rw-r--r--  1 meren  staff   859K Jul 31 23:56 Enterococcus_faecalis_6240.fa.gz
-rw-r--r--  1 meren  staff   2.9M Jul 31 23:56 Enterococcus_faecalis_6240.h5
-rw-r--r--  1 meren  staff   4.7M Jul 31 23:56 Enterococcus_faecalis_6250.db
-rw-r--r--  1 meren  staff   863K Jul 31 23:56 Enterococcus_faecalis_6250.fa.gz
-rw-r--r--  1 meren  staff   2.9M Jul 31 23:56 Enterococcus_faecalis_6250.h5
-rw-r--r--  1 meren  staff   4.5M Jul 31 23:56 Enterococcus_faecalis_6255.db
-rw-r--r--  1 meren  staff   835K Jul 31 23:56 Enterococcus_faecalis_6255.fa.gz
-rw-r--r--  1 meren  staff   2.8M Jul 31 23:56 Enterococcus_faecalis_6255.h5
-rw-r--r--  1 meren  staff   4.2M Jul 31 23:56 Enterococcus_faecalis_6512.db
-rw-r--r--  1 meren  staff   772K Jul 31 23:56 Enterococcus_faecalis_6512.fa.gz
-rw-r--r--  1 meren  staff   2.6M Jul 31 23:56 Enterococcus_faecalis_6512.h5
-rw-r--r--  1 meren  staff   5.1M Jul 31 23:56 Enterococcus_faecalis_6557.db
-rw-r--r--  1 meren  staff   946K Jul 31 23:56 Enterococcus_faecalis_6557.fa.gz
-rw-r--r--  1 meren  staff   3.2M Jul 31 23:56 Enterococcus_faecalis_6557.h5
-rw-r--r--  1 meren  staff   4.3M Jul 31 23:56 Enterococcus_faecalis_6563.db
-rw-r--r--  1 meren  staff   793K Jul 31 23:56 Enterococcus_faecalis_6563.fa.gz
-rw-r--r--  1 meren  staff   2.7M Jul 31 23:56 Enterococcus_faecalis_6563.h5
-rw-r--r--  1 meren  staff   4.6M Jul 31 23:56 Enterococcus_faecium_6589.db
-rw-r--r--  1 meren  staff   852K Jul 31 23:56 Enterococcus_faecium_6589.fa.gz
-rw-r--r--  1 meren  staff   2.9M Jul 31 23:56 Enterococcus_faecium_6589.h5
-rw-r--r--  1 meren  staff   4.9M Jul 31 23:56 Enterococcus_faecium_6590.db
-rw-r--r--  1 meren  staff   915K Jul 31 23:56 Enterococcus_faecium_6590.fa.gz
-rw-r--r--  1 meren  staff   3.1M Jul 31 23:56 Enterococcus_faecium_6590.h5
-rw-r--r--  1 meren  staff   4.6M Jul 31 23:56 Enterococcus_faecium_6601.db
-rw-r--r--  1 meren  staff   860K Jul 31 23:56 Enterococcus_faecium_6601.fa.gz
-rw-r--r--  1 meren  staff   2.9M Jul 31 23:56 Enterococcus_faecium_6601.h5
-rw-r--r--  1 meren  staff   4.3M Jul 31 23:56 Enterococcus_faecium_6778.db
-rw-r--r--  1 meren  staff   806K Jul 31 23:56 Enterococcus_faecium_6778.fa.gz
-rw-r--r--  1 meren  staff   2.7M Jul 31 23:56 Enterococcus_faecium_6778.h5
-rw-r--r--  1 meren  staff   4.2M Jul 31 23:56 Enterococcus_faecium_6798.db
-rw-r--r--  1 meren  staff   773K Jul 31 23:56 Enterococcus_faecium_6798.fa.gz
-rw-r--r--  1 meren  staff   2.6M Jul 31 23:56 Enterococcus_faecium_6798.h5
{% endhighlight %}

There also are two files in the `additional-files/pangenomics` directory to describe how to access to the external genomes:

|name|contigs_db_path|
|:--|:--|
|E_faecalis_6240|external-genomes/Enterococcus_faecalis_6240.db|
|E_faecalis_6250|external-genomes/Enterococcus_faecalis_6250.db|
|E_faecalis_6255|external-genomes/Enterococcus_faecalis_6255.db|
|E_faecalis_6512|external-genomes/Enterococcus_faecalis_6512.db|
|E_faecalis_6557|external-genomes/Enterococcus_faecalis_6557.db|
|E_faecalis_6563|external-genomes/Enterococcus_faecalis_6563.db|
|E_faecium_6589|external-genomes/Enterococcus_faecium_6589.db|
|E_faecium_6590|external-genomes/Enterococcus_faecium_6590.db|
|E_faecium_6601|external-genomes/Enterococcus_faecium_6601.db|
|E_faecium_6778|external-genomes/Enterococcus_faecium_6778.db|
|E_faecium_6798|external-genomes/Enterococcus_faecium_6798.db|

and the internal one:

|name|bin_id|collection_id|profile_db_path|contigs_db_path|
|:--|:--:|:--:|:--|:--|
|E_faecalis_SHARON|E_faecalis|E_faecalis|../../PROFILE.db|../../CONTIGS.db|

It is this simple to combine MAGs and isolates.

{:.notice}
You can find a comprehensive tutorial on the anvi'o pangenomic workflow [here]({% post_url anvio/2015-11-14-pangenomics %}).

So everything is ready for an analysis. To start the pangenomic workflow, you can run `anvi-pan-genome` program the following way (which takes about 5 minutes on my laptop):

{% highlight bash %}
 $ anvi-pan-genome -i additional-files/pangenomics/internal-genomes.txt \
                   -e additional-files/pangenomics/external-genomes.txt \
                   --num-threads 4 \
                   -o pan
{% endhighlight %}

This will generate an output directory with many files, including a sub-directory that contains all the files to run the anvi'o interactive interface to visualize results.

If you initiate the the interactive interface the following way,

{% highlight bash %}
 $ anvi-interactive -t pan/pan/tree.txt \
                    -d pan/pan/view_data.txt \
                    -s pan/pan/samples.db \
                    -A pan/pan/additional_view_data.txt \
                    -p pan/pan/profile.db \
                    --title "Enterococcus Pan" \
                    --manual
{% endhighlight %}

To get this ugly looking display:

[![E. facealis pan](images/e-faecalis-pan.png)](images/e-faecalis-pan.png){:.center-img .width-50}

I did it prettier for you. If you kill the server, and import the state file the following way, and re-run the server,

{% highlight bash %}
 $ anvi-import-state -p pan/pan/profile.db --state additional-files/state-files/state-pan.json --name default
 $ anvi-interactive -t pan/pan/tree.txt \
                    -d pan/pan/view_data.txt \
                    -s pan/pan/samples.db \
                    -A pan/pan/additional_view_data.txt \
                    -p pan/pan/profile.db \
                    --title "Enterococcus Pan" \
                    --manual
{% endhighlight %}

It will look much more reasonable:

[![E. facealis pan](images/e-faecalis-pan-state.png)](images/e-faecalis-pan-state.png){:.center-img .width-70}

Now not only we can see how our E. facealis genome to what is available, we can also see that it is not missing or carrying a great number of proteins compared to other genomes. The clustering of genomes based on protein clusters indicate that it is most similar to the genome `Enterococcus faecalis 6250`, which, according to the `00_INFO_ABOUT_EXTERNAL_GENOMES.txt` under `05_PANGENOMICS` directory, corresponds to the assembly ID [ASM28119v1](https://www.ncbi.nlm.nih.gov/gquery/?term=ASM28119v1) if you were to be interested in exploring further.


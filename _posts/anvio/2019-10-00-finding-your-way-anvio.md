---
layout: post
title: "Finding Your Way Around anvi'o"
excerpt: "`anvi-search` and other ways to navigate this beautiful maze of a software."
modified: 2019-10-05
authors: [iva]
categories: [anvio]
comments: true
---

{% include _toc.html %}

{% include _join-anvio-slack.html %}


Anvi'o does a lot. And sometimes, it can be a little overwhelming. We get it. We often get overwhelmed by it all, too. That's why we are writing this post, to remind everyone (including ourselves) how to navigate around anvi'o.

If you need help finding a command, exploring anvi'o's capabilities, or solving an issue, this post is for you. <3

# The `anvi-search` command you never knew you needed
You know what you want to do with your data, and you think you can do it with anvi'o. But how do you find an anvi'o program to do it?

Enter `anvi-search`. Recently implemented in anvi'o v-6 by our very own Evan Kiefl, this program allows you to query available anvi'o programs using keywords. It is really quite beautiful.

Let's see an illustrative example. Suppose you want to annotate some contigs with hits from the <a href=https://pfam.xfam.org/>Pfam database</a>. To search for a program that does this in anvi'o, you could use the command `anvi-search pfam`. And after thinking about itself for a bit, anvi'o will return to you a table that looks like the following:
```
╒══════════════════╤═════════════════════════════════╤═════════════╤════════════╤════════════╕
│ Program          │ Description                     │ Resources   │ Provides   │ Requires   │
╞══════════════════╪═════════════════════════════════╪═════════════╪════════════╪════════════╡
│ anvi-run-pfams   │ Run Pfam on Contigs Database.   │             │ functions  │ contigs-db │
│                  │                                 │             │            │            │
│                  │                                 │             │            │ pfams-data │
├──────────────────┼─────────────────────────────────┼─────────────┼────────────┼────────────┤
│ anvi-setup-pfams │ Download Pfam data from the EBI │             │ pfams-data │            │
╘══════════════════╧═════════════════════════════════╧═════════════╧════════════╧════════════╛
```
As you can see, there are two anvi'o programs related to the keyword "pfam". From the descriptions in the table, you can see that `anvi-run-pfams` will do what you want, but that it requires you to have a contigs-db and pfams-data available. You then might suppose that you would need to run `anvi-setup-pfams` first to get the pfams-data (and you would be right). If you didn't already know how to create a contigs-db, you could use `anvi-search --provides contigs-db` to find out, and you would receive the following:
```
╒═══════════════════════════╤════════════════════════════════════════╤═════════════╤════════════╤═══════════════╕
│ Program                   │ Description                            │ Resources   │ Provides   │ Requires      │
╞═══════════════════════════╪════════════════════════════════════════╪═════════════╪════════════╪═══════════════╡
│ anvi-gen-contigs-database │ Generate a new anvio contigs database. │             │ contigs-db │ contigs-fasta │
╘═══════════════════════════╧════════════════════════════════════════╧═════════════╧════════════╧═══════════════╛
```
And there you have it! From two quick searches, you know how to go from the FASTA file of your contigs to a Pfam annotation, using three anvi'o programs.

By default, the output of `anvi-search` tells you, for each matching anvi'o program:
- the command (which you can query directly using the `-n` or `--name` flag)
- a brief description of what the program does
- resources like tutorials that may help you understand or use the command (if available)
- a list of what the program produces (which you can query directly using the `-p` or `--provides` flag)
- a list of what the program needs (which you can query directly using the `-r` or `--requires` flag)

You can modify this output as you like by providing a comma-separated list of headers to the `-R` or `--report` flag. For example, if you want to see all the anvi'o programs that are tagged for "pangenomics" (along with what they provide and require), you would use the following command: `anvi-search -R Requires,Provides,Tags pangenomics`
This would result in a table like the following:
```
╒═══════════════════════════════════════════╤════════════════════╤════════════════╤═════════════════╕
│ Program                                   │ Requires           │ Provides       │ Tags            │
╞═══════════════════════════════════════════╪════════════════════╪════════════════╪═════════════════╡
│ anvi-get-enriched-functions-per-pan-group │ misc-data-layers-  │ functional-    │ pangenomics     │
│                                           │ category           │ enrichment-txt │                 │
│                                           │                    │                │ functions       │
│                                           │ pan-db             │                │                 │
│                                           │                    │                │                 │
│                                           │ genomes-storage-db │                │                 │
├───────────────────────────────────────────┼────────────────────┼────────────────┼─────────────────┤
│ anvi-run-workflow                         │                    │                │ metagenomics    │
│                                           │                    │                │                 │
│                                           │                    │                │ phylogenomics   │
│                                           │                    │                │                 │
│                                           │                    │                │ contigs         │
│                                           │                    │                │                 │
│                                           │                    │                │ pangenomics     │
│                                           │                    │                │                 │
│                                           │                    │                │ metapangenomics │
╘═══════════════════════════════════════════╧════════════════════╧════════════════╧═════════════════╛
```


# Exploring the anvi'o network
Anvi'o is a web of interconnected programs and the databases, data sources, and anvi'o-specific concepts that these programs depepnd upon. We have a nice <a href="http://merenlab.org/software/anvio/vignette/">visualization of this network</a> that is useful for exploring what anvi'o can do and what it needs to do it.

At first glance, the network looks complicated:

<a href="/software/anvio/network/" target="_blank"><img src="/images/anvio-network.png" width="100%" /></a>

But you can click on any one icon to learn more about how it fits into the rest of the anvi'o universe. For example, let's say you wanted to learn more about what you can do with contigs databases. You'd click on the contigs-db icon, and peruse the vast web of anvi'o programs that use or generate these databases:

<img src="/images/anvio-network-contigs-db.png" width="100%" />

# Learning resources
We try our best to keep our online documentation up-to-date. Sometimes, our best is not very good (we're sorry! we have lots of things to do :disappointed:), but most of the time, you should be able to find plenty of information right here on our website. This section gives a brief overview of different kinds of learning resources we provide.

## Anvi'o vocabulary
In most of the documents referenced below, you will find a whole host of anvi'o-specific terms. These can make learning anvi'o initially confusing and seemingly impenetrable. Which is why we developed a <a href="http://merenlab.org/vocabulary/#all-things-anvio">vocabulary page</a> to help you get started. Trust us - just as with any scientific discipline, picking up the terminology will ultimately deepen your understanding of and facilitate your scientific conversations about anvi'o-related analyses.

## The anvi'o vignette
The <a href="http://merenlab.org/software/anvio/vignette/">vignette</a> is the place to go for a comprehensive list of information on all programs in the latest stable release of anvi'o. It shows much of the same information that `anvi-search` and individual program help pages (`-h`) would give you, sometimes with extra tidbits of explanation. Googling an anvi'o program usually directs you to the vignette page. But if you are one of those adventurous folks who is using the development version of anvi'o (good for you!), some of the newer programs may not be on this page yet.

## Tutorials
Anvi'o includes various sets of programs that are meant to be used together to accomplish a larger analysis task. We often write tutorials to explain these common use cases. These are step-by-step instructions on how to get from your initial data (FASTA files, BAM files, etc) to visualizations of biological insights, in anvi'o. They can be accessed from the "Tutorials" pull-down menu on each page of the anvi'o section of our website.

Tutorials under the 'More Abstract' category are generic recipes. They focus on the commands that should be run to get from A to B, without discussion of the finer details of data interpretation and visualization at each step. These are better suited for people who know what they want to do and generally how they could do it, but need to figure out how to perform their analysis in anvi'o.

In contrast, tutorials under the 'More Practical' category are geared towards more comprehensive instruction, using a specific dataset as an example. The datasets are downloadable, allowing you to follow along by running the commands on your own computer and comparing the output you get to that shown in the tutorial. These are suitable for people who know what they want to do but not how to do it (for people who don't know what they want at all, please see the Microbial 'Omics section below. Or a therapist.) Caveats, expert tips, and visualization strategies that go along with each step are often discussed, making these tutorials excellent for those who want to understand each step of the process and need some context to get it.

{:.notice}
Please note that tutorials are often written around the time when their core programs are first published. This (unfortunately) often makes them a snapshot of how to do an analysis with a particular version of anvi'o, rather than an evolving document. We update them occasionally, but often don't have time to go through all the finer details. The consequences of this are that specific commands might not work for you, depending on how much anvi'o has been updated since the tutorial was written. If that happens, don't worry! The overall analysis workflow usually remains the same, and there are plenty of resources to help you figure out how to modify the command. Please feel free to send us a comment, email, or slack message about it - we can help you out, and then update our tutorial so that other people don't struggle with the same issue. :)

## Workshops
We periodically give workshops for in-person training on anvi'o. Yes! You can come meet us! We will be thrilled to make your acquaintance. See a list of upcoming workshops, and a description of past workshops, <a href="http://merenlab.org/2016/08/18/events/">here</a.

## Microbial 'Omics
We are always searching for ways to make our field of science more accessible to everyone. For those who want to learn more about microbial 'omics in general, you can check out Meren's <a href="http://merenlab.org/momics/">classroom resources</a> on the topic. The 'omics-specific <a href="http://merenlab.org/vocabulary/#all-things-omics">vocabulary page</a> may also be helpful.


# HALP I NEED ADVICE T.T (How to get help on using anvi'o)
Anvi'o has a fairly extensive user community. There are a few different ways to reach out to us (the developers) or others in this community if you need help with using anvi'o.

These options are for situations in which you:
- are not sure how to do a particular analysis in anvi'o, and cannot find any hints or documentation about it in our online resources. You can check with us to see if what you want is possible in anvi'o. Maybe it is, and we just haven't documented it yet (or well). In this case, we can usually suggest some ways to do it. Or maybe it isn't possible (yet), and in this case, we may be open to implementing something for you (if you ask nicely :) ).
- are having trouble running or installing anvi'o. These are situations in which you are getting nicely formatted error messages from anvi'o, which usually means something is not right with your command usage or computational environment. (If you are not getting nice error messages, please refer to the technical help section below). To make things easier for us, please tell us at least the full command that you were using and the entire error message. It also often useful for us to know your anvi'o version, for instance via the output of `anvi-interactive -v`.
- are not sure how to set parameters of an anvi'o program to accomplish your particular goal, or have any other such question about more nuanced usage of anvi'o. These questions often benefit from the collective experience of the anvi'o community.
- are confused by any part of our online documentation. We are sincerely sorry if this is the case, and we will be happy to clarify for you.

If you don't receive an answer in a timely fashion using any of the methods below, please do not be discouraged. Send another request for help. You will get an answer eventually.

## Google group
We have a discussion group, served by Google Groups, wherein you can ask your anvi'o related questions and receive answers via email. The nice thing about this group is that all messages are archived, so everyone can benefit from past discussions. The answers you seek may already be out there. :)

Please refer to <a href="http://merenlab.org/2015/10/14/anvio-discussion-group/">this post</a> for details on how to join and participate in the Google group.

## Slack
For those of us who are into IM-ing (that's instant messaging for those who are not up on their old school internet slang), join our Slack channel.

{% include _join-anvio-slack.html %}

## Web page comments
Most of our blog posts, tutorials, and other web pages have a comments section at the bottom where you can ask your questions (or give your feedback) on issues related to the topic discussed in the post/tutorial/page. This is often very useful to the rest of the community, as it collects relevant discussions right there where other people will go to look for that information.

# AAAH I FOUND A BUG 0.0 (How to get help with technical issues)
If anvi'o has crashed on you, then you have found a bug :bug: :sob:. Now you can help us to fix it. Please <a href="https://github.com/merenlab/anvio/issues">submit an issue</a> on Github to inform us of the problem. It is most helpful to us if you include the following information in the issue:
- your anvi'o environment (ie, version information)
- your computational platform (ie, operating system)
- step-by-step instructions on how to reproduce the problem
- complete error messages/crash reports
We will fix the bug as soon as we can, and close the Github issue once it has been resolved.

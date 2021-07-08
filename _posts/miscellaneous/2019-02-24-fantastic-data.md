---
layout: post
title: "Fantastic Data and How to Share Them: A Plea to Journal Editors and Reviewers"
excerpt: "A small checklist for reviewers to follow to make public data sources public and forever"
modified: 2019-02-24
tags: []
authors: [meren]
categories: [miscellaneous]
comments: true
redirect_from: /sharing-data/
---


{:.notice}
A version of this blog post [first appeared](https://researchdata.springernature.com/users/113363-a-murat-eren-meren/posts/44023-fantastic-data-and-how-to-share-them-a-plea-to-journal-editors-and-reviewers) at [the Research Data Blog at Springer Nature](https://researchdata.springernature.com/) as a part of the "[Love Data Week 2019](https://researchdata.springernature.com/users/184393-roza-sakellaropoulou/posts/44241-a-week-of-lovin-on-data-content-for-love-data-week-2019)". When I was asked to contribute why data sharing is important from our perspective, I thought it would be easy for me to write about some positive success stories, yet I couldn't steer away from writing one of the serious issues we have been recently struggling with instead: poorly shared data (unintentionally, or intentionally).

Love for data comes naturally but sharing data can be tricky. Besides technical challenges, there are a lot of emotions involved. Imagine the journey of a scientist with their novel data. Creating with so much care something so unique, enjoying the privilege to be the first to interpret their sophisticated nature, constantly updating a trusted few about every surprising observation, but then finally having to let them go out to the world where the data can have their own adventures with others. With wisdom comes great data, with love comes great anxiety. I can’t be the only one who sees parallels between generating data and becoming a parent. Tough business. But while my mother gets no citations for who I am or any of the stories in which I partake, fortunately those who let their data roam the world freely are given a lot of credit: the data-mediated love is reciprocal between scientists. I for one remember every scientist whose public data helped me understand more and am thankful for their efforts.

My <a href="http://merenlab.org/">group</a> at the University of Chicago Department of Medicine studies the ecology and evolution of microbes. We strive to understand how microbes do all the things they do in their natural habitats, whether those habitats are mammalian guts, oral cavities, surface oceans, or ovaries of tiny mosquitoes. This is an extremely data-heavy discipline. We frequently see new publications with billions of sequences coming from hundreds of samples. In comparison to the heroes of our field, we consume considerably more data than we generate. Yes, we are among those who live on others’ data. But we are proud. When we hear some scientists describe our kind as <a href="http://doi.org/10.1056/NEJMe1516564">research parasites</a>, we <a href="http://researchparasite.com/">laugh</a> and continue to do our best to bring new perspectives and reveal previously unrecognized insights in data that should be publicly available for <a href="https://www.statnews.com/2015/12/23/sharing-data-science/">honest science</a>.

While my group may not be generating as much data for others, we do create computational tools for everyone. Our open-source software platform, anvi’o, <a href="https://naturemicrobiologycommunity.nature.com/users/113363-a-murat-eren-meren/posts/34040-microbiologists-vs-shotgun-metagenomes-surface-ocean">empowers us and other scientists</a> to study complex ‘omics data. We write <a href="http://merenlab.org/software/anvio/">tutorials</a> and implement detailed <a href="http://merenlab.org/data/">reproducible workflows</a> to show the community how we mix our data with data from others. Besides the immense satisfaction we get from what we do, our place in all this helps us to better understand how to generate accessible software and accessible data that are meant to be public. In that respect data and software have some remarkable parallels as well as puzzling differences.

For instance, this sentence will work equally well either with data or software: 

<blockquote>
We got to create this particular ___________ because we have convinced a committee somewhere that spending our scarce public resources doing this was going to address questions, the significance of which were beyond our personal needs or curiosity

<div class="blockquote-author">data or software</div>
</blockquote>

In contrast, the following applies most often to data but almost never to software:

<blockquote>
We created this particular ___________ promising that it would be publicly available, but the way we shared it lacks a crucial piece of information and requires anyone intending to use it to contact us first.

<div class="blockquote-author">data or software</div>
</blockquote>

While this outcome is not always intentional, it dramatically limits the accessibility of a resource despite the intentions of its funders and benefactors. Luckily, in most cases the unintentional omission of critical details of data can be addressed if the reviewers of data-releasing studies are more generous with their guidance.

Those who are not familiar with science may rightfully assume that in this competitive world scientists who let their data take part in stories of others without expecting anything in return must be so rare. Fortunately, this is not true. There are many scientists who share their data completely and enthusiastically, and we know who they are. There are perhaps even more scientists who would be willing to share their data, but they lack proper guidance. And finally, there is always a smaller group of those who shall not be named.

Years of mining large numbers of publications for data to meet the requirements of our ‘parasitic’ lifestyle taught us a couple of things. If we set aside the first group towards whom we are forever thankful and the last group towards whom we are generally politically correct, what is it that we can do as editors and reviewers to help those in the middle group to honor the efforts of the first, and the investment of public resources into science? What follows is our two cents answer to this question in the context of molecular sequencing data.

## Dear Editors

If you are handling manuscripts for a journal with appropriate open data policies, please do not solicit reviewers for a manuscript unless all sequencing data is properly cited with previous publications or new accession numbers. The editorial desk is the first place that can avoid many issues downstream in less than a minute. If the data seem to be in place, please consider including a data-savvy reviewer in your solicitation and asking them specifically to evaluate whether the data meet the requirements of what constitutes publicly available data that is truly reusable without requiring any input from its creators in the future.

## Dear Reviewers

Having accepted to evaluate a study with sequencing data, you are the ambassador of the science community and its future generations (so, no pressure at all). The effort you will invest during your evaluation will guarantee the availability of data even if its creators quit science, decide to not respond to their e-mails ever again, or are simply not interested in seeing other interpretations of their data despite our best interest. Please put on your Batman outfit, and consider following this checklist:

* **Step 1**: Do not review for journals that lack open data policies. Let them gracefully become extinct.
 
* **Step 2**: Decline to review a manuscript that does not make available all the data it uses. No data no cake. You are not there to review fiction.

* **Step 3**: If all accession numbers are present, make sure they are in fact functional by downloading at least one sample, and please take a quick look at the file in your terminal.

* **Step 4**: If the accession numbers resolve to meaningful files, please make sure the manuscript includes a Supplementary Data Table where each accession number is associated with a sample name, and each sample name is associated with environmental variables. And while at it, please make sure this information is not embedded into a PDF file as text, but it can be copy-pasted, sorted, or searched with ease using standard tools.

* **Step 5**: Please take an extra minute to make sure the sample names in the Supplementary Data Table match to sample names shown in figures and mentioned in the text.

These steps could save hundreds of hours of extra work graduate students and post-doctoral researchers invest into mining public data and prevent tremendous amount of waste of our public resources as most of the data that is not shared properly never get to be reanalyzed for additional insights.

Love for data comes naturally, but sharing data appropriately is not so easy. Scientists who generate datasets from valuable samples often let their precious resource go on to their own journey without them, even when they are certain that they have not learned everything that can be learned from those data. Some do it better than others, but those of us who benefit from public data are eternally thankful to all who try.

---

<blockquote class="twitter-tweet" data-lang="en"><p lang="en" dir="ltr">New Blog Post: &quot;Fantastic Data and How to Share Them: A Plea to Journal Editors and Reviewers&quot;.<br><br>Or &quot;How to be an Excellent Reviewer and Guarantee the Accessibility of Public Data&quot; 😇<a href="https://t.co/HcgUWkSSU1">https://t.co/HcgUWkSSU1</a> <a href="https://twitter.com/hashtag/OpenData?src=hash&amp;ref_src=twsrc%5Etfw">#OpenData</a> <a href="https://twitter.com/hashtag/OpenScience?src=hash&amp;ref_src=twsrc%5Etfw">#OpenScience</a></p>&mdash; A. Murat Eren (Meren) (@merenbey) <a href="https://twitter.com/merenbey/status/1099859811682959360?ref_src=twsrc%5Etfw">February 25, 2019</a></blockquote>
<script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script>

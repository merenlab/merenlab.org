---
layout: post
title: "High-Resolution Microbiome Informatics Meeting, 2015"
excerpt: "A small-scale meeting to make a statement about high-resolution."
modified: 2015-03-14 
tags: [high-resolution]
categories: [oligotyping, med]
comments: true
---

I believe there is now a broad understanding among microbial ecologists that the analysis of high-throughput sequencing of marker gene amplicons from microbiomes with 3% OTUs is somewhat limiting, and we can all benefit from more advanced tools that go beyond the apparent limitations of OTU clustering.

Probably the vast majority of scientist who work with microbiomes would concur that the 3% cut-off is an arbitrary threshold, and due to their phylogenetically-mixed nature 3% OTUs can trap multiple ecologically distinct organisms into one unit, which can be quite detrimental in some cases. Although in general it seems to be working OK, in fact one nucleotide difference at the marker gene-level is not any less important than say, ten. We were told a long time ago that one nucleotide difference is not something to ignore [[1]](http://mmbr.asm.org/content/62/4/1353.full), and further evidence supporting this suggestion from the high-throughput sequencing-based studies (along with alternative solutions) is accumulating rapidly [[2](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12114/full), [3](http://www.pnas.org/content/111/28/E2875.full), [4](http://www.nature.com/ismej/journal/v9/n1/full/ismej2014117a.html), [5](http://www.nature.com/ismej/journal/vaop/ncurrent/full/ismej2014195a.html)]. Indeed 3% cut-off serves more to bioinformaticians than microbiologists and microbial ecologists, however, there is a growing number of bioinformaticians who are working hard to correct this while mitigating the artificial inflation of observed units in high-throughput sequencing datasets.

A while ago [Gary Borisy](http://forsyth.org/person/scientist/gary-borisy) had adviced that this may be a good time for microbial ecologists and methods developers to come together and discuss advantages and pitfalls of emerging high-resolution approaches that attempt to circumvent the limitations of 3% OTUs. Understanding that this would be a long-term goal, Jessica Mark Welch, myself, and Jack Gilbert took his advice, and applied for a collaboration grant to at least start the conversation. Our proposal to host a small group of scientist in Woods Hole for a short meeting was [funded](http://www.mbl.edu/blog/first-recipients-of-mbl-uchicagoargonne-collaboration-award-announced/) by a [Marine Biological Laboratory - University of Chicago Collaboration Award](http://www.mbl.edu/research/mbl-uchicagoargonne-collaboration-awards/), and now it is happening. Very soon.

### The Meeting

The meeting will take place on the 19<sup>th</sup> and 20<sup>th</sup> of March in Woods Hole. I know this is quite a short notice, but if you are around, you are welcome to attend and be a part of the discussion.

Here is a list of speakers and their presentations (you can download the meeting program [here]({{ site.url }}/files/hrmi-2015-program.pdf) as well):

* **[Sarah P. Preheim](http://engineering.jhu.edu/dogee/faculty/sarah-preheim/)** (Johns Hopkins University)<br />“*Distribution-based clustering: using ecology to refine the operational taxonomic unit*”
* **Mikhail Tikhonov** (Harvard University)<br />“*Interpreting 16S metagenomic data without clustering to achieve sub-OTU resolution*”* **Amnon Amir** (University of Colorado) <br />“*Deblurring 16S sequencing data*”* **A. Murat Eren** (Marine Biological Laboratory)<br />“*Oligotyping and MED: Using information theory to find ecologically meaningful units in marker gene amplicons*”

* **[Jed A. Fuhrman](http://dornsife.usc.edu/labs/fuhrmanlab)** (University of Southern California)<br /> “*Benefits of higher phylogenetic resolution in marine microbial ecology*”* **[Jack A. Gilbert](http://pondside.uchicago.edu/ecol-evol/people/gilbert.html)** (University of Chicago)<br /> “*Invisible influence: examining microbial activity*”* **[Pawel Gajer](http://medschool.umaryland.edu/facultyresearchprofile/viewprofile.aspx?id=20207)** (University of Maryland) <br />“*Ribotyping with Restricted Boltzmann Machines*”* **Inés Martínez** (University of Alberta)<br />“*What determines microbiome responsiveness towards a non-digestible oligosaccharide?*”* **[Floyd E. Dewhirst](http://forsyth.org/person/scientist/floyd-dewhirst)** (Forsyth Institute)<br />“*The human oral microbiome database: Importance of speciation in NGS era*”

* **[Jeremiah J. Faith](http://www.mountsinai.org/profiles/jeremiah-faith)** (The Mount Sinai Hospital)<br />“*Experimental advances for improving microbiome strain-tracking resolution and the computational opportunities they present*”* **[Jessica L. Mark Welch](http://www.mbl.edu/jbpc/staff/markwelchj/)** (Marine Biological Laboratory)<br />“*Connecting high-resolution informatics to spatial organization through imaging*”* **[Christopher Quince](http://search.warwick.ac.uk/profile?id=MTM3NTAzOGYxODJlYmU%3D)** (University of Warwick)<br />“*Reconstruction of genomes from metagenomes using coverage and composition*”



Yes, marker gene-related talks take a good portion of the schedule. And yes, the need for methods that can deliver highly resolved depictions of microbial ecosystems through sequencing data is not limited to marker gene analyses. However, there are multiple published alternative methods available today for the analysis of marker gene data that can compete or complement the canonical ones, and this is a good time to discuss how can we improve them and make them more accessible.

I hope we will all be convinced at the end of this meeting that we should do this again, and agree that the next schedule should cover more base in the context of metagenomics, and imaging.

---

I am planning to recap everything and post it online shortly after the meeting.

Meanwhile, in case you have suggestions for us to discuss during this meeting, or you know methods that are not yet widely-known, please send me an e-mail and let me know, so I can do my best to include them and their developers into our future conversations.

As a separate note, this is Woods Hole just a week ago:

<a href="{{ site.url }}/images/eel-pond-5.jpg"><img src="{{ site.url }}/images/eel-pond-5.jpg"></a>

The yellow lines that secure dinghies during the summer halfheartedly hold the ice that covers Eel Pond. There is still snow on the ground here, but I hope there will be none of it left before March 19<sup>th</sup>.

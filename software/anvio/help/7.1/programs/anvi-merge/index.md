---
layout: page
title: anvi-merge [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-merge
image:
  featurerelative: ../../../images/header.png
  display: true
---

Merge multiple anvio profiles.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[single-profile-db](../../artifacts/single-profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[misc-data-items-order](../../artifacts/misc-data-items-order) <img src="../../images/icons/CONCEPT.png" class="artifact-icon-mini" /></span></p>


## Usage


The main function of `anvi-merge` is to convert multiple <span class="artifact-n">[single-profile-db](/software/anvio/help/7.1/artifacts/single-profile-db)</span>s into a single <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span> (also called a merged profile database). Basically, this takes the alignment data from each sample (each contained in its own <span class="artifact-n">[single-profile-db](/software/anvio/help/7.1/artifacts/single-profile-db)</span>) and combines them into a single database that anvi'o can look through more easily. 

### Overview: How to run anvi-merge

1. Set up your <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>. See that page for more information

1. Use <span class="artifact-n">[anvi-profile](/software/anvio/help/7.1/programs/anvi-profile)</span> to create a <span class="artifact-n">[single-profile-db](/software/anvio/help/7.1/artifacts/single-profile-db)</span> for each of your samples (formatted into a <span class="artifact-n">[bam-file](/software/anvio/help/7.1/artifacts/bam-file)</span>) *(Note: for each of these runs, you'll need to use the same <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> and parameters)*

1. Use `anvi-merge` to merge those <span class="artifact-n">[single-profile-db](/software/anvio/help/7.1/artifacts/single-profile-db)</span>s into a single database, called a <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span>. This will look something like the following:

<div class="codeblock" markdown="1">
anvi&#45;merge &#45;c cool_contigs.db \
            Single_profile_db_1 Single_profile_db_2 \
            &#45;o cool_contigs_merge
</div>
                    
This will put all of the output files (the final <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span> as well as a <span class="artifact-n">[misc-data-items-order](/software/anvio/help/7.1/artifacts/misc-data-items-order)</span> which is the result of your hierarchical clustering and describes the order to display your contigs in) into the folder `cool_contigs_merge `.
    

## Other Parameters

You must give `anvi-merge` your contigs database and single profile databases. However, you can also provide more information or give addtional instructions. Use the flag `-h` at any time to display the help menu.

### Hierarchical Clustering 

#### To run or not to run? 
* Use the flag `--skip-hierarchical-clustering` to turn hierarchical clustering off. This will save on computation time, but will skip out on creating the tree of contigs at the center of the interactive interface. If you have more than 25,000 splits in the final profile, this will be set automatically. 
* Use the flag `--enforce-hierarchical-clustering` to turn hierarchical clustering back on. This will take a long time, but will produce a lovely contigs tree for the interactive interface. 

#### Additional parameters
* Provide a custom distance metric for clustering using the flag `--distance.` (The default is euclidean)
* Provide a custom linkage method for clustering using the flag `--linkage.` (The default is ward)

### Providing additional information
* Provide the sample name using the flag `-S`. If you don't, anvi'o will come up with one, but it probably won't be any good. 
* Provide a text file in markdown to describe the project using the flag `--description`. This will show up when you later use the interactive interface to analyze your profile-db. 

### Output Information
* Provide an output destination with the flag `-o`.
* Add the flag `-W` to overwrite existing files in that directory. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-merge.md) to update this information.


## Additional Resources


* [Another description as part of the metagenomic workflow](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-profile)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-merge) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

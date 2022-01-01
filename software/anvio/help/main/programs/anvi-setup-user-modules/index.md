---
layout: page
title: anvi-setup-user-modules [program]
categories: [anvio]
comments: false
redirect_from: /m/anvi-setup-user-modules
image:
  featurerelative: ../../../images/header.png
  display: true
---

Set up user-defined metabolic pathways into an anvi&#x27;o-compatible database.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/ivagljiva.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Iva Veseli</span><div class="page-author-social-box"><a href="mailto:iveseli@uchicago.edu" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://github.com/ivagljiva" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[user-modules-data](../../artifacts/user-modules-data) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[modules-db](../../artifacts/modules-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[user-modules-data](../../artifacts/user-modules-data) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>


## Usage


This program creates a <span class="artifact-n">[modules-db](/software/anvio/help/main/artifacts/modules-db)</span> out of a set of user-defined metabolic modules, for use by <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/main/programs/anvi-estimate-metabolism)</span>.

It takes as input a directory containing module files for each user-defined module, formatted in the same way as KEGG modules are. It parses these modules into the `USER_MODULES.db` database. This directory of user-defined data is referred to as <span class="artifact-n">[user-modules-data](/software/anvio/help/main/artifacts/user-modules-data)</span>, and the help page for that artifact contains a detailed account of how to create your own module definitions and estimate their completeness.

This page will give a few details specific to running <span class="artifact-n">[anvi-setup-user-modules](/software/anvio/help/main/programs/anvi-setup-user-modules)</span>.

### Default Usage

To run this program, you must provide an input directory containing your module definitions:

<div class="codeblock" markdown="1">
anvi&#45;setup&#45;user&#45;modules &#45;&#45;user&#45;modules /path/to/user/data/directory
</div>

This input directory must have a specific format (see section below). The `USER_MODULES.db` will be generated in this directory, so you can use the same path to provide your data to <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/main/programs/anvi-estimate-metabolism)</span> when you want to estimate completeness for these modules.

### Input directory format

The directory you provide to the `--user-modules` parameter must have another folder inside of it, which must be called `modules`. Inside that `modules` folder, you should put text files containing the definitions of your metabolic modules - one file per module. The file should be named according to the identifier you want the module to have, and should not have any extension.

Here is an example schematic of a proper input directory:
```
MY_METABOLISM_DATA_DIR
  |
  |- modules
        |- U00001
        |- U00002
        |- U00003
        |- U00004
```
The `U0000x` files in the schematic above each contains a definition for one module. Running `anvi-setup-user-modules --user-modules MY_METABOLISM_DATA_DIR` will produce a `USER_MODULES.db` file in the `MY_METABOLISM_DATA_DIR` folder which contains 4 modules named U00001, U00002, U00003, and U00004 (assuming those files are formatted correctly).

### How do I format the module files?

We use KEGG's system for describing metabolic modules, so you will need to format your metabolic pathways in the same way. Here is an example, for a module file called `U00002` (like in the schematic above):
```
ENTRY       U00002
NAME        Nitrogen fixation (full Nif gene set)
DEFINITION  K02588+K02586+K02591-K00531 K02587 K02592 K02585
ORTHOLOGY   K02588  NifH
            K02586  NifD
            K02591  NifK
            K00531  anfG
            K02587  NifE
            K02592  NifN
            K02585  NifB
CLASS       User modules; Energy metabolism; Nitrogen metabolism
ANNOTATION_SOURCE  K02588  KOfam
                    K02586  KOfam
                    K02591  KOfam
                    K00531  KOfam
                    K02587  KOfam
                    K02592  KOfam
                    K02585  KOfam
///
```
As you can see, there are different data types in the file, named by the all-capital word at the beginning of the line (we call this the 'data name'). The second column of the file is the value corresponding to that type of information ('data value'). Some data names, like ORTHOLOGY and ANNOTATION_SOURCE, also have a 3rd column further defining the data value (which we call the 'data definition'). Each field in the file should be separated by _at least two spaces_. And the file must end with '///' on the last line (don't ask us why).

The data names you see in the example above are the minimum you should include to define the module. Here is a bit more information about each type of data:
- ENTRY: this is the identifier for the module. It can be anything you want, but should be just one word (underscores and dashes allowed). It should also be the same as the name of the module file. Importantly, this identifier should not be the same as any KEGG module, or you will get an error during setup.
- NAME: this is the name of the metabolic pathway, which can be any arbitrary string (spaces allowed)
- DEFINITION: this is the set of enzymes required for the reactions in the metabolic pathway. The enzymes should be identified by their accession numbers in their respective annotation source - in the example above, these are all KOfams, so the enzyme accessions are KO numbers. However, you can use enzymes from any annotation source you like (COGs, Pfams, custom HMMs, etc), as long as you have a way to annotate them in your contigs database. The rules for defining a metabolic pathway in the KEGG fashion are described in the [technical details section](https://merenlab.org/software/anvio/help/main/programs/anvi-estimate-metabolism/#what-data-is-used-for-estimation) of the help page for <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/main/programs/anvi-estimate-metabolism)</span>, so please read through that section for help with designing your pathway definition.
- ORTHOLOGY: this section maps the identifier for an enzyme (the second column, or data value) to the functional definition of that enzyme (the third column, or data definition). You need one of these lines for every enzyme in your module DEFINITION line.
- CLASS: this line categorizes the module. It must be one string with three sections, separated by semi-colons. The first section is the 'class' of the module, the second section is the 'category' of the module, and the third section is the 'sub-category'. Feel free to use the existing KEGG categories to describe your pathway, or to make up something entirely new. Or, if you don't care about this categorization at all, you could just put random strings in each section (as long as you have two semi-colons in the string, you will be golden).
- ANNOTATION_SOURCE: this section maps the identifier for an enzyme (the second column, or data value) to its annotation source (the third column, or data definition). You need one of these lines for every enzyme in your module DEFINITION line. The annotation source must match the functional annotation source in the contigs database that is associated with the enzyme's annotations. For instance, KOfams annotated with <span class="artifact-n">[anvi-run-kegg-kofams](/software/anvio/help/main/programs/anvi-run-kegg-kofams)</span> have source 'KOfam' (as above), the 2020 COG source from <span class="artifact-n">[anvi-run-ncbi-cogs](/software/anvio/help/main/programs/anvi-run-ncbi-cogs)</span> is 'COG20_FUNCTION', the source for custom HMM profiles given to <span class="artifact-n">[anvi-run-hmms](/software/anvio/help/main/programs/anvi-run-hmms)</span> is the `--hmm-source` directory name, and so on.

You can also define other data names, if you want. Some common ones that can be found in KEGG modules are COMPOUND, REACTION, PATHWAY, COMMENT, REFERENCE, and AUTHORS; but you are not limited by the ones used by KEGG.

Why must we format the module files this way, you ask? Well, to be honest, KEGG modules are formatted like this, and our infrastructure for working with that data has simply been adapted to work with arbitrary, user-defined data. KEGG makes the rules :)

### Specifying KEGG data to be used for sanity checking

If you haven't yet run <span class="artifact-n">[anvi-setup-kegg-kofams](/software/anvio/help/main/programs/anvi-setup-kegg-kofams)</span> on your computer, you will get an error when you try to run this program. This is because KEGG data is always used in addition to user-defined modules, and we need to be aware of which KEGG modules exist so we can make sure none of the user-defined modules have the same identifiers as these.

By default, this program looks for the KEGG data in the default location, so if you have set up KEGG data in a non-default directory, you should specify the path to that directory using the `--kegg-data-dir` parameter:

<div class="codeblock" markdown="1">
anvi&#45;setup&#45;user&#45;modules &#45;&#45;user&#45;modules /path/to/user/data/directory &#45;&#45;kegg&#45;data&#45;dir /path/to/KEGG/data/directory
</div>

If you have multiple KEGG data directories on your computer, you should specify the one that you intend to use (along with this user-defined data) for <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/main/programs/anvi-estimate-metabolism)</span> downstream. It will not make a difference if all of your modules have identifiers unique from KEGG ones, but just in case they overlap, it is better to catch this during setup rather than later during metabolism estimation. :)


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-setup-user-modules.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-setup-user-modules) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

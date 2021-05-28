---
layout: post
title: "Installing third-party software"
excerpt: "Recipes to install various software tools anvi'o uses"
modified: 2016-06-18
categories: [anvio]
comments: true
---

This article describes some recipes for the installation of softwre used by anvi'o but developed by other groups.

{:.warning}
**If you are an end-user and trying to install anvi'o, you shouldn't be on this page**. Instead, please find our installation instructions and follow them.

But the article is limited to those tools that we do not install through our conda recipes. Please see our GitHub repository to see a complete list of third-party software that are [required for anvi'o to install](https://github.com/merenlab/anvio/blob/master/conda-recipe/anvio-minimal/meta.yaml) or [critical but not required for its runtime](https://github.com/merenlab/anvio/blob/master/conda-recipe/anvio/meta.yaml) (i.e., they are critical if you want to do anything comprehensive with anvi'o, but without them you still can install and run anvi'o).

{% include _toc.html %}

{:.notice}
{% include _fixthispage.html source="_posts/anvio/2016-06-18-installing-third-party-software.md" %}

## DSSP

[DSSP](https://swift.cmbi.umcn.nl/gv/dssp/index.html) _(Dictionary of Secondary Structure Prediction) is a program for annotating amino acids in protein structures with structural information, such as secondary structure (alpha, helix, etc.) and solvent accessibility_.

**Citation**: [https://www.ncbi.nlm.nih.gov/pubmed/6667333](https://www.ncbi.nlm.nih.gov/pubmed/6667333)

**Citation**: [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3013697/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3013697/)

Do you already have DSSP? Go to your terminal, and type `mkdssp --version`. If you do not get an error, read ahead to make sure your DSSP is operational. If you get an error, you need to install DSSP.

Use this command line to install DSSP with `conda`:

```
conda install -c salilab dssp
```

If you get a `command not found` error, you have to install `conda`, and easy install instructions can be found on [their website](https://conda.io/docs/user-guide/install/index.html#)

Test your DSSP installation by running the following commands your terminal:

```bash
wget http://files.rcsb.org/view/1H97.cif #download myoglobin structure
mkdssp -i 1H97.cif -o myoglobin_DSSP.txt
```

If no error is produced, DSSP is working (you can open up `myoglobin_DSSP.txt` to see what the output looks like).


## MODELLER

[MODELLER](https://salilab.org/modeller/) is a program for homology or comparative modeling of protein three-dimensional structures using a reference database.

**Citation**: [https://www.ncbi.nlm.nih.gov/pubmed/18428767](https://www.ncbi.nlm.nih.gov/pubmed/18428767)

**Citation**: [https://www.ncbi.nlm.nih.gov/pubmed/10940251](https://www.ncbi.nlm.nih.gov/pubmed/10940251)

**Citation**: [https://www.ncbi.nlm.nih.gov/pubmed/8254673](https://www.ncbi.nlm.nih.gov/pubmed/8254673)

Go to your terminal and type `mod` and then attempt to autocomplete the command by pressing the tab key twice. If something like `mod9.xx` pops up, where `xx` is some number, you have MODELLER. If `xx` is lower than `19`, you should reinstall MODELLER using the following instructions.

If you don't have any version of MODELLER, their website has excellent cross-platform [intallation instructions](https://salilab.org/modeller/release.html#install). If you are using conda, you can just run:

```
conda install -c salilab modeller
```

Regardless of your method of installation, check it worked by typing `mod` into your terminal and seeing that it autocompletes after pressing the tab key twice.

Along with an installation, you will need a license key (which is free for academic use). You can go ahead and [register for a license key](https://salilab.org/modeller/registration.html) on their website and put the key where they tell you, **or** you can just pretend you have one and anvi'o will give you precise instructions on how to get one and what to do with it when the time comes.

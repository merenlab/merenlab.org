---
layout: post
authors: [meren]
title: "Installing the oligotyping pipeline"
excerpt: "Installing the pipeline shall not be a PITA."
modified: 2014-08-16
tags: [pip]
categories: [oligotyping]
comments: true
---


{:.notice}
There is also a [docker image for oligotyping]({% post_url oligotyping/2014-09-02-virtualbox %}) if you are interested in a painless run.

Installing the oligotyping pipeline on a MAC or a Linux computer is as easy as writing this command in the terminal window (assuming you already have [pip](https://pypi.python.org/pypi/pip)):

    sudo pip install oligotyping

However, the installation requires a couple of software to be installed on your system. If you feel really really lazy, you can try the command above, and then take a look at this document if you get errors. In an ideal world, I would prefer to read this document and follow the steps mentioned to make sure things are in order.

The oligotyping pipeline is now very easy to install on a MAC computer (I am not going to explain it for Linux users, assuming they can take a look at this document and figure out what needs to be done). But please do not try to download installers online to install required software if I haven’t specifically requested that.

I will first provide instructions to install the **latest stable pipeline**, then I will add some information for the ones who wish to use the **latest unstable snapshot** that is being developed on [GitHub code repository for oligotyping](https://github.com/meren/oligotyping).

{: .notice}
**Note for MAC users**: I assume that you already have [XCode and XCode Command Line Tools](https://developer.apple.com/xcode/downloads/) are installed on your system. Without them your MAC can’t compile things, so please install them first. I stumbled upon a case where missing XCode or XCode Command Line Tools did cause issues. You should also make sure you accept the XCode license agreement, as is mentioned [here](http://apple.stackexchange.com/questions/175069/how-to-accept-xcode-license/), you can do that by running this in your terminal:
    
    sudo xcodebuild -license

{: .notice}
**Another note for MAC users**: To avoid conflicts during multithreading, you should run this command in your terminal to tell matplotlib to use the Agg backend: ```echo "backend: Agg" >> ~/.matplotlib/matplotlibrc```. You will get an error if you don't have the directory ```~/.matplotlib/``` in place. In that case you should first create it (hint: ```mkdir ~/.matplotlib/```), and then run the same command again.

When the installation is done, you can refer to [this post]({% post_url oligotyping/2012-05-11-oligotyping-pipeline-explained %}) to learn more about the use of the pipeline.

If you have questions, please send them to the [discussion forum](https://groups.google.com/forum/#!forum/oligotyping).


## Installing the latest stable version using conda + pip (easy-peasy)

For this to work, you need [miniconda](https://docs.conda.io/en/latest/miniconda.html) to be installed on your system. If you are not sure whether it is installed or not, open a terminal (hopefully an [iTerm](https://www.iterm2.com/), if you are using Mac) and type `conda`. You should see an output like this instead of a 'command not found' error (your version might be different):

```bash
$ conda --version
conda 4.8.3
```

If you don't have conda installed, then you should first install it through their [installation page](https://docs.conda.io/en/latest/miniconda.html). Once you have confirmed you have conda installed, run this command to make sure you are up-to-date:

``` bash
conda update conda
```

{:.warning}
Please make sure you create a new conda environment oligotyping (you can make sure you are not in a conda environment by opening a new terminal and running `conda deactivate`. You can see which environments exist on your computer by running `conda env list`.

Then, create a new environment for oligotyping:

``` bash
conda create -y --name oligotyping python=3.7
```

Once it is ready, activate your new environment:

``` bash
conda activate oligotyping
```

And finally install the following packages first:

``` bash
conda install -y r-base blast -c bioconda
```

Once everything is done, install the following R libraries:

```
conda install -y r-ggplot2 r-vegan r-gplots r-gtools r-reshape r-optparse r-pheatmap r-rcolorbrewer r-compute.es
```

If you came this far, it means you have installed everything necessary.

You can either execute these last two steps, OR directly jump to the next section to follow the active codebase (which would be my suggestion).

```
pip install oligotyping
```

When you type this command you should get a help menu, instead of a command not found error:

    oligotype --help


## Tracking the active oligotyping codebase without installation (artistic mode)

This is to make sure you are following the very latest state of the codebase.

First, setup your environment for oligotyping:

``` bash
conda update conda
conda create -y --name oligotyping python=3.7
conda activate oligotyping
```

Now it is time to get a copy of the oligotyping codebase. Here I will suggest `~/github/` as the base directory, but you can change if you want to something else (in which case you must remember to apply that change all the following commands, of course):

``` bash
# setup the code directory and get the oligotyping codebase
mkdir -p ~/github && cd ~/github/
git clone https://github.com/merenlab/oligotyping.git
conda activate oligotyping
```

Let's install some helper tools next:

``` bash
conda install -y r-base blast -c bioconda
conda install -y r-ggplot2 r-vegan r-gplots r-gtools r-reshape r-optparse r-pheatmap r-rcolorbrewer r-compute.es -c conda-forge
```

When this is done successfully, you can run the following to install oligotyping Python dependencies:

```
pip install -r ~/github/oligotyping/requirements.txt
```

Now we will setup your conda environment in such a way, every time you activate oligotyping within it, you will get the very latest updates from the main repository. Just copy-paste the entire thing below into your terminal:


```bash
cat <<EOF >${CONDA_PREFIX}/etc/conda/activate.d/oligotyping.sh
# creating an activation script for the the conda environment for anvi'o
# development branch so (1) Python knows where to find anvi'o libraries,
# (2) the shell knows where to find anvi'o programs, and (3) every time
# the environment is activated it synchronizes with the latest code from
# active GitHub repository:
export PYTHONPATH=\$PYTHONPATH:~/github/oligotyping/
export PATH=\$PATH:~/github/oligotyping/bin:~/github/oligotyping/sandbox
echo -e "\033[1;34mUpdating from oligotyping GitHub \033[0;31m(press CTRL+C to cancel)\033[0m ..."
cd ~/github/oligotyping && git pull && cd -
EOF
```

{:.warning}
If you are using zsh by default these may not work. If you run into a trouble here or especially if you figure out a way to make it work both for zsh and bash, please let us know.

If everything worked, you should be able to type the following commands in a new terminal and see similar outputs:

```bash
meren ~ $ conda activate oligotyping
Updating from oligotyping GitHub (press CTRL+C to cancel) ...
Already up to date.
/Users/meren
(oligotyping) meren ~ $

(oligotyping) meren ~ $ which oligotype
/Users/meren/github/oligotyping/bin/oligotype

(oligotyping) meren ~ $ decompose
usage: decompose [-h] [-m FLOAT] [-X] [-d INTEGER] [-A INTEGER] [-M INTEGER]
                 [-V INTEGER] [-t CHAR] [-o OUTPUT_DIRECTORY] [-p STR] [-g]
                 [-S] [-H] [-R] [-F] [-K] [-T] [-N INTEGER] [-E FILEPATH]
                 [--skip-gen-html] [--skip-gen-figures]
                 [--skip-check-input-file] [--skip-gexf-files] [--quick]
                 [--version]
                 FILEPATH
decompose: error: the following arguments are required: FILEPATH
```

If that is the case, you’re all set.

If you followed these instructions, every time you open a terminal you will have to run the following command to activate your oligotyping environment:

```
conda activate oligotyping
```

{:.notive}
We hope the oligotyping serves you well. It's development has stalled since our group's interests have shifted from amplicon data to 'omics. There are many ways to improve especially the minimum entropy decomposition algorithm. Due to the core principle behind it (i.e., the information theory), with only relatively minor improvements this algorithm will almost certainly perform better to recover subtle biological signal in complex data than any other algorithm that relies on pairwise sequence alighments to count mismatches. If you are interested, plese send us an e-mail, and we will give you our best ideas so you can improve the code and make it your own.

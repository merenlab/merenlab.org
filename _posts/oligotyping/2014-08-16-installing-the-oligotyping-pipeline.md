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

{% include _toc.html %}

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


## Installing the latest stable version using conda + pip (suggested method)

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
conda create -y --name oligotyping python=3.6
```

Once it is ready, activate your new environment:

``` bash
conda activate oligotyping
```

And finally install the following packages first:

``` bash
conda install -y git r-base blast pip
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


## Using the oligotyping pipeline directly from the upstream without installation (Meren's way)

This is to make sure you are following the very latest state of the codebase.

First, do this:

``` bash
pip install virtualenv
```

Now it is time to get a copy of the oligotyping codebase. Here I will suggest `~/github/` as the base directory, but you can change if you want to something else (in which case you must remember to apply that change all the following commands, of course):

``` bash
# setup the code directory and get the oligotyping codebase
mkdir -p ~/github && cd ~/github/
git clone https://github.com/merenlab/oligotyping.git
```

Here we will setup a directory to keep the Python virtual environment for OLIGOTYPING (virtual environment within a virtual environment, keep your totem nearby):

``` bash
mkdir -p ~/virtual-envs/

rm -rf ~/virtual-envs/oligotyping

virtualenv ~/virtual-envs/oligotyping --pip 20.0.2
source ~/virtual-envs/oligotyping/bin/activate
```

When this is done successfully, you can run the following two lines to install oligotyping Python dependencies and deactivate the Python environment:

```
pip install -r ~/github/oligotyping/requirements.txt
deactivate
```

Now we will setup your conda environment in such a way, every time you activate OLIGOTYPING within it, you will get the very latest updates from the `master` repository:

``` bash
# updating the activation script for the Python virtual environmnet
# so (1) Python knows where to find OLIGOTYPING libraries, (2) BASH knows
# where to find its programs, and (3) every the environment is activated
# it downloads the latest code from the `master` repository
echo -e "\n# >>> OLIGOTYPING STUFF >>>" >> ~/virtual-envs/oligotyping/bin/activate
echo 'export PYTHONPATH=$PYTHONPATH:~/github/oligotyping/' >> ~/virtual-envs/oligotyping/bin/activate
echo 'export PATH=$PATH:~/github/oligotyping/bin:~/github/oligotyping/sandbox' >> ~/virtual-envs/oligotyping/bin/activate
echo 'cd ~/github/oligotyping && git pull && cd -' >> ~/virtual-envs/oligotyping/bin/activate
echo "# <<< OLIGOTYPING STUFF <<<" >> ~/virtual-envs/oligotyping/bin/activate
```

Finally we define an alias, `oligotyping-activate-master`, so when you are in your conda environment for `oligotyping` you can run it as a command to initiate everything like a pro:

```
echo -e "\n# >>> OLIGOTYPING STUFF >>>" >> ~/.bash_profile
echo 'alias oligotyping-activate-master="source ~/virtual-envs/oligotyping/bin/activate"' >> ~/.bash_profile
echo "# <<< OLIGOTYPING STUFF <<<" >> ~/.bash_profile
source ~/.bash_profile
```

At this stage if you run `oligotyping-activate-master`, you should be able to run this command without any errors:

```
$ oligotype -h

$ which oligotype
~/github/oligotyping/bin/oligotype
```

If that is the case, you're golden.


So, this is the end of setting up the active OLIGOTYPING codebase on your computer so you can follow our daily additions to the code before they appear in stable releases, and use OLIGOTYPING exactly the way we use on our computers for our science on your own risk.

**Please note that given this setup so far, every time you open a terminal you will have to first activate conda, and then the Python environment:**

```
conda activate oligotyping
oligotyping-activate-master
```

You can always use `~/.bashrc` or `~/.bash_profile` files to add aliases to make these steps easier for yourself, or remove them when you are tired.

<details markdown="1"><summary>Show/hide Meren's BASH profile setup</summary>

This is all personal taste and they may need to change from computer to computer, but I added the following lines at the end of my `~/.bash_profile` to easily switch between different versions of OLIGOTYPING on my Mac system:

{:.notice}
If you are using Anaconda rather than miniconda, or you are using Linux and not Mac, you will have to find corresponding paths for lines that start with `/Users` down below :)

``` bash

init_oligotyping_master () {
    {
        deactivate && conda deactivate
    } &> /dev/null

    export PATH="/Users/$USER/miniconda3/bin:$PATH"
    . /Users/$USER/miniconda3/etc/profile.d/conda.sh
    conda activate oligotyping
    oligotyping-activate-master
    export PS1="\[\e[0m\e[40m\e[1;30m\] :: oligotyping master :: \[\e[0m\e[0m \[\e[1;34m\]\]\w\[\e[m\] \[\e[1;31m\]>>>\[\e[m\] \[\e[0m\]"
}

alias om=init_oligotyping_master
```

With this setup, in a new terminal window I can type `om` to activate oligotyping.


**But please note** that both aliases run `deactivate` and `conda deactivate` first, and they may not work for you especially if you have a fancy setup. I'd be very happy to improve these shortcuts.
</details>




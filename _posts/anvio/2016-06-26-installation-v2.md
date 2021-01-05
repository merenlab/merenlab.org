---
layout: post
title: "Installing anvi'o"
excerpt: "Instructions to install the v2 branch of the platform."
modified: 2019-05-14
tags: []
categories: [anvio]
redirect_from:
 - /2015/05/01/installation/
 - /install-anvio/
comments: true
image:
  feature: https://github.com/merenlab/anvio/raw/master/anvio/data/interactive/images/logo.png
---

{% include _toc.html %}

{% include _project-anvio-version.html %}

This article explains basic steps of installing anvi'o using rather conventional methods both for end users and current of future developers.

<details markdown="1"><summary>Show/hide A docker solution for those who are in a hurry</summary>

We do recommend you to install anvi'o on your system as explained below, but **if you just want to run anvi'o without any installation**, you can actually do it within minutes using [docker](https://docs.docker.com/get-docker/).

The docker solution is very simple, guaranteed to work, and very effective to do quick analyses or visualize anvi'o data currencies from others without having to install anything. A more detailed article on how to run anvi'o in docker [is here]({% post_url anvio/2015-08-22-docker-image-for-anvio %}), but here is a brief set of steps.

Assuming you have docker installed and running on your computer, first pull the container:

``` bash
docker pull meren/anvio:7
```

This step will take a few minutes and require about 15Gb disk space. Once it is done, you can run it the following way:

```
docker run --rm -it -v `pwd`:`pwd` -w `pwd` -p 8080:8080 meren/anvio:7
```

And that's it! You are now in a virtual environment that runs anvi'o. You can exit this environment by pressing `CTRL+D`.

{:.warning}
If you wish to do resource demanding analyses, don't forget to increase CPU and memory resources allocated for anvi'o using the docker Preferences menu.

If you at some point want to remove all containers and reclaim all the storage space, you can run this after exiting all containers:

```
docker system prune --force -a
```
</details>

Please consider opening an <a href="https://github.com/meren/anvio/issues">issue</a> for technical problems, or join us on Slack if you need help:

{% include _join-anvio-slack.html %}

{:.notice}
{% include _fixthispage.html source="_posts/anvio/2016-06-26-installation-v2.md" %}

{:.warning}
We thank [Daan Speth](https://twitter.com/daanspeth), [Jarrod Scott](https://orcid.org/0000-0001-9863-1318), [Susheel Bhanu Busi](https://scholar.google.com/citations?user=U0g3IzQAAAAJ&hl=en), and [Mike Lee](https://twitter.com/AstrobioMike), who kindly invested their time to test the installation instructions on this page on different systems and/or made suggestions to the document to ensure a smoother installation experience for everyone.

## (1) Setup conda

This is a very simple and effective way to install anvi'o on your system along with most of its dependencies.

{:.notice}
<div id="windowsreminder" style="min-height: 50px;" markdown="1">
<img src="/images/windows10.png" style="float:left; width:50px; border: None; display: block; margin-right: 5px;">Although these installation instructions primarily target and rigorously tested for Linux and Mac OSX, you will be able to follow them if you are using Microsoft Windows **if and only if you first install the [Linux Subsystem for Windows](https://docs.microsoft.com/en-us/windows/wsl/install-win10)**. Our users have reported success stories with Ubuntu on WSL.
</div>

**For this to work, you need [miniconda](https://docs.conda.io/en/latest/miniconda.html) to be installed on your system.** If you are not sure whether it is installed or not, open a terminal (such as [iTerm](https://www.iterm2.com/), if you are using Mac) and type `conda`. You should see an output like this instead of a 'command not found' error (your version might be different):

```bash
$ conda --version
conda 4.9.2
```

If you don't have conda installed, then you should first install it through their [installation page](https://docs.conda.io/en/latest/miniconda.html). Once you have confirmed you have conda installed, run this command to make sure you are up-to-date:

``` bash
conda update conda
```

Good? Good! You are almost there!

## (2) Setup an anvi'o environment

{:.notice}
It is a good idea to **make sure you are not already in a conda environment** before you run the following steps. Just to be clear, you can indeed install anvi'o in an existing conda environment, but if things go wrong, we kindly ask you to refer to meditation for help, rather than [anvi'o community resources]({% post_url anvio/2019-10-07-getting-help %} since there is no way we can help you if you are installing anvi'o in a different conda environment :) If you want to see what environments do you have on your computer and whether you already are in one of them in your current terminal by running `conda env list`. **If all these are too much for you and all you want to do is to move on with the installation**, simply do this: open a new terminal, and run `conda deactivate`, and continue with the rest of the text.

You have two options here.

### A shortcut worth trying

{:.notice}
If this shortcut doesn't work for you for some reason, you will not worry and try the slower option down below.

Resolving dependencies (especially on Mac systems) can take a very long time for conda (which is a [known problem](https://github.com/conda/conda/issues/7239)), hence, here we will use a serious shortcut to generate an environment for anvi'o `v{% include _project-anvio-version-number.html %}`.


First you will need to get a copy of the following file, but you have two options depending on the operating system you're using.

If you are using **Mac OSX**, use this one:

``` bash
curl https://merenlab.org/files/anvio-conda-environments/anvio-environment-7-MACOS.yaml \
     --output anvio-environment-7.yaml
```

If you are using **Linux/Windows**, use this one:

``` bash
curl https://merenlab.org/files/anvio-conda-environments/anvio-environment-7-LINUX.yaml \
     --output anvio-environment-7.yaml
```

Run this to make sure you don't already have an environment called `anvio-v{% include _project-anvio-version-number.html %}`:

```
conda env remove --name anvio-7
```

Now create a new `anvio-v{% include _project-anvio-version-number.html %}` environment using the file you just downloaded:

```
conda env create -f anvio-environment-7.yaml
```

If this didn't go well, jump to the slower but reliable option. If it did go well, then you should activate that environment,

```
conda activate anvio-7
```

and jump to "[Download and install anvi'o](#3-install-anvio)".

### Slower but reliable option

First create a new conda environment:

``` bash
conda create -y --name anvio-7 python=3.6
```

And activate it:

```
conda activate anvio-7
```

Now you are in a pristine environment, in which you will install all conda packages that anvi'o will need to work properly. This looks scary, but it will work if you just copy paste it and press ENTER:

``` bash
conda install -y -c bioconda "sqlite >=3.31.1"
conda install -y -c bioconda prodigal
conda install -y -c bioconda mcl
conda install -y -c bioconda muscle
conda install -y -c bioconda hmmer
conda install -y -c bioconda diamond
conda install -y -c bioconda blast
conda install -y -c bioconda megahit
conda install -y -c bioconda spades
conda install -y -c bioconda bowtie2
conda install -y -c bioconda bwa
conda install -y -c bioconda samtools
conda install -y -c bioconda centrifuge
conda install -y -c bioconda trimal
conda install -y -c bioconda iqtree
conda install -y -c bioconda trnascan-se
conda install -y -c bioconda r-base
conda install -y -c bioconda r-stringi
conda install -y -c bioconda r-tidyverse
conda install -y -c bioconda r-magrittr
conda install -y -c bioconda r-optparse
conda install -y -c bioconda bioconductor-qvalue
conda install -y -c bioconda fasttree
conda install -y -c conda-forge h5py=2.8.0

# this may cause some issues. if it doesn't install,
# don't worry:
conda install -y -c bioconda fastani
```

Now you can jump to "[Download and install anvi'o](#3-install-anvio)"!


## (3) Install anvi'o

Here you will first download the Python source package for the official anvi'o release:

```
curl -L https://github.com/merenlab/anvio/releases/download/v7/anvio-7.tar.gz \
        --output anvio-7.tar.gz
```

And install it using `pip` like a boss:

```
pip install anvio-7.tar.gz
```

If everything went fine, you can jump to "[Check your anvi'o setup](#4-check-your-installation)" to see if things worked for you, and then you are free to go!


## (4) Check your installation

If you are here, you are ready to check if everything is working on your system.

The easiest way to do it is to run a self-test and see if everything is in order:

``` bash
anvi-self-test
```

{:.notice}
If you don't want anvi'o to show you a browser window at the end and quietly finish testing if everything is OK, add `--no-interactive` flag to the command above.

{:.warning}
It is absolutely normal to see 'warning' messages. In general anvi'o is talkative as it would like to keep you informed. In an ideal world you should keep a careful eye on those warning messages, but in most cases they will not require action.

If everything goes smoothly, your browser should pop-up and show you something like this at the end of `anvi-self-test`:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/misc/mini-test-screenshot.png"><img src="{{ site.url }}/images/anvio/misc/mini-test-screenshot.png" width="50%" /></a>
</div>

{:.notice}
The screenshot above is from 2015 and will be vastly different from the [interactive interface]({% post_url anvio/2016-02-27-the-anvio-interactive-interface %}) you should see in your browser. It is still here so we remember where we came from 😇

All fine? Perfect! Now you have a computer that can run anvi'o!

### Final words

Here are some reminders we hope you can quickly go through before leaving this page:

* **Your tests looked like they worked but your browser didn't pop-up?** OK. There may be to reasons for this. **If you are on your personal computer**, this may mean that the Python setup on your system is unable to find your default browser. Not a biggie. When an interactive interface is initiated successfully, you will see the address through which you can access to the interactive. Copy that address and paste it in your browser's address bar. **If you are connected to a remote computer**, it may take a bit more setup to make to use your local browser to access to an anvi'o interactive interactive interface running remotely. [Read this article]({% post_url anvio/2018-03-07-working-with-remote-interative %}) (or ask your systems administrator to read it) to learn how you can forward displays from servers to your personal computer. Finally, **Windows users may be greeted by a "This site can't be reached" error** even when their browsers show up :/ You can overcome this problem by copy-pastig the address you see in your terminal to your browser address bar (this is what happens when you bring a Microsoft Windows to an open-source fight).

* **When you open a new terminal you get command not found error when you run anvi'o commands?**  Depending on your conda setup, you may need to activate the anvi'o conda environment every time you open a new terminal window. In these cases you will either need to run `source activate anvio-7` or `conda activate anvio-7` (this assumes you named your conda environment for anvio `anvio-7`, but you can always list your conda environments by typing `conda env list`.

* **You are not using [Chrome](https://www.google.com/chrome/) as your default browser and anvi'o complains about it?** We hate the idea of asking you to change your browser preferences for anvi'o :( But currently, Chrome maintains the most efficient SVG engine among all browsers we tested as of 2021. For instance, Safari can run the anvi'o interactive interface, however it takes orders of magnitude more time and memory compared to Chrome. Firefox, on the other hand, doesn't even bother drawing anything at all. Long story short, the anvi'o interactive interface __will not perform optimally__ with anything but Chrome. So you need Chrome. Moreover, if Chrome is not your default browser, every time interactive interface pops up, you will need to copy-paste the address bar into a Chrome window. You can learn what is your default browser by running this command in your terminal:

<div markdown="1" style="margin-top: -32px; margin-left: 80px;">
``` bash
python -c 'import webbrowser as w; w.open_new("http://")'
```
</div>

If you are here, you are done! Congratulations, and thank you very much for your patience!

Now you can take a look up some anvi'o resources [here]({{ site.url }}/software/anvio), or come say hi to us on Slack.

{% include _join-anvio-slack.html %}


## (5) Follow the active development (you're a wizard, arry)

{:.warning}
This section is not meant to be followed by those who would define themselves as *end users* in a conventional sense. But we are not the kinds of people who would dare to tell you what you can and cannot do. FWIW, our experience suggests that if you are doing microbiology, you will do computers no problem if you find this exciting.

If you follow these steps, you will have anvi'o setup on your system in such a way, every time you initialize your anvi'o environment you will get **the very final state of the anvi'o code**. Plus, you can have both the stable and active anvi'o on the same computer.

Nevertheless, it is important to keep in mind that there are multiple advantages and disadvantages to working with the active development branch. Advantages are obvious and include,

* **Full access to all new features and bug fixes in real-time**, without having to wait for stable releases to be announced.

* A working system to **hack anvi'o and/or add new features to the code** (this strategy is exactly how we develop anvi'o and use it for our science at the same time at our lab). 

In contrast, disadvantages include,

* **Unstable intermediate states may frustrate you with bugs, and in extremely rare instances loss of data** (this happened only once so far during the last five years, and required one of our users to re-generate their contigs databases).

* Difficulty to mention the anvi'o version in a paper for reproducibility. Although this can easily be solved by sharing not the version number of anvi'o but the cryptographic hash of the last commit for reproducibility. If you ever struggle with this, please let us know and we will help you.

If you are still here, let's start.

---

First make sure you are not in any environment by running `conda deactivate`. Then, make sure you don't have an environment called `anvio-dev` (as in *anvi'o development*):

```
conda env remove --name anvio-dev
```

Now we can continue with setting up the conda environment.

### Setting up the conda environment

First create a new conda environment:

``` bash
conda create -y --name anvio-dev python=3.6
```

And activate it:

```
conda activate anvio-dev
```

Install necessary packages:

``` bash
conda install -y -c bioconda "sqlite >=3.31.1"
conda install -y -c bioconda prodigal
conda install -y -c bioconda mcl
conda install -y -c bioconda muscle
conda install -y -c bioconda hmmer
conda install -y -c bioconda diamond
conda install -y -c bioconda blast
conda install -y -c bioconda megahit
conda install -y -c bioconda spades
conda install -y -c bioconda bowtie2
conda install -y -c bioconda bwa
conda install -y -c bioconda samtools
conda install -y -c bioconda centrifuge
conda install -y -c bioconda trimal
conda install -y -c bioconda iqtree
conda install -y -c bioconda trnascan-se
conda install -y -c bioconda r-base
conda install -y -c bioconda r-stringi
conda install -y -c bioconda r-tidyverse
conda install -y -c bioconda r-magrittr
conda install -y -c bioconda r-optparse
conda install -y -c bioconda bioconductor-qvalue
conda install -y -c bioconda fasttree
conda install -y -c conda-forge h5py=2.8.0

# this may cause some issues. if it doesn't install,
# don't worry:
conda install -y -c bioconda fastani
```
Now you are ready for the code.

### Setting up the local copy of the anvi'o codebase

If you are here, it means you have a conda environment with everything except anvi'o itself. We will make sure this environment _has_ anvi'o by getting a copy of the anvi'o codebase from GitHub.

Here I will suggest `~/github/` as the base directory to keep the code, but you can change if you want to something else (in which case you must remember to apply that change all the following commands, of course). Setup the code directory:

``` bash
mkdir -p ~/github && cd ~/github/
```

Get the anvi'o code:

{:.warning}
If you only plan to follow the development branch you can skip this message. But if you are not an official anvi'o developer but intend to change anvi'o and send us pull requests to reflect those changes in the official repository, you may want to clone anvi'o from your own fork rather than using the following URL. Thank you very much in advance and we are looking forward to seeing your PR!

```
git clone --recursive https://github.com/merenlab/anvio.git
```

Now it is time to install the Python dependencies of anvi'o:

``` bash
cd ~/github/anvio/
pip install -r requirements.txt
```

Now all dependencies are in place, and you have the code. One more step.

### Linking conda environment and the codebase

Now we have the codebase and we have the conda environment, but they don't know about each other.

Here we will setup your conda environment in such a way that every time you activate it, you will get the very latest updates from the main anvi'o repository. While you are still in anvi'o environment, copy-paste these lines into your terminal:

``` bash
cat <<EOF >${CONDA_PREFIX}/etc/conda/activate.d/anvio.sh
# creating an activation script for the the conda environment for anvi'o
# development branch so (1) Python knows where to find anvi'o libraries,
# (2) the shell knows where to find anvi'o programs, and (3) every time
# the environment is activated it synchronizes with the latest code from
# active GitHub repository:
export PYTHONPATH=\$PYTHONPATH:~/github/anvio/
export PATH=\$PATH:~/github/anvio/bin:~/github/anvio/sandbox
echo -e "\033[1;34mUpdating from anvi'o GitHub \033[0;31m(press CTRL+C to cancel)\033[0m ..."
cd ~/github/anvio && git pull && cd -
EOF
```

{:.warning}
If you are using `zsh` by default these may not work. If you run into a trouble here or especially if you figure out a way to make it work both for `zsh` and `bash`, please let us know.

If everything worked, you should be able to type the following commands in a new terminal and see similar outputs:

```
meren ~ $ conda activate anvio-dev
Updating from anvi'o GitHub (press CTRL+C to cancel) ...

(anvio-dev) meren ~ $ which anvi-self-test
/Users/meren/github/anvio/bin/anvi-self-test

(anvio-dev) meren ~ $ anvi-self-test -v
Anvi'o .......................................: hope (v7-dev)

Profile database .............................: 35
Contigs database .............................: 20
Pan database .................................: 14
Genome data storage ..........................: 7
Auxiliary data storage .......................: 2
Structure database ...........................: 2
Metabolic modules database ...................: 2
tRNA-seq database ............................: 1

(anvio-dev) meren ~ $
```

If that is the case, you're all set.

Every change you will make in anvi'o codebase will immediately be reflected when you run anvi'o tools (but if you change the code and do not revert back, git will stop updating your branch from the upstream). 

If you followed these instructions, every time you open a terminal you will have to run the following command to activate your anvi'o environment:

```
conda activate anvio-dev
```

If you are here, you can now jump to "[Check your anvi'o setup](#4-check-your-installation)" to see if things worked for you using `anvi-self-test`.


## Bonus: An alternative BASH profile setup

{:.notice}
This section is written by Meren and reflects his setup on a Mac system that runs miniconda where `bash` is [setup as the default shell](https://itnext.io/upgrading-bash-on-macos-7138bd1066ba). If you are using another shell and if you would like to share your solution, please send a PR!

This is all personal taste and they may need to change from computer to computer, but I added the following lines at the end of my `~/.bash_profile` to easily switch between different versions of anvi'o on my Mac system:


``` bash
# This is where my miniconda base is, you can find out
# where is yours by running this in your terminal:
#
#    conda env list | grep base
#
export MY_MINICONDA_BASE="/Users/$USER/miniconda3"

init_anvio_7 () {
    deactivate &> /dev/null
    conda deactivate &> /dev/null
    export PATH="$MY_MINICONDA_BASE/bin:$PATH"
    . $MY_MINICONDA_BASE/etc/profile.d/conda.sh
    conda activate anvio-7
    export PS1="\[\e[0m\e[47m\e[1;30m\] :: anvi'o v7 :: \[\e[0m\e[0m \[\e[1;32m\]\]\w\[\e[m\] \[\e[1;31m\]>>>\[\e[m\] \[\e[0m\]"
}


init_anvio_dev () {
    deactivate &> /dev/null
    conda deactivate &> /dev/null
    export PATH="$MY_MINICONDA_BASE/bin:$PATH"
    . $MY_MINICONDA_BASE/etc/profile.d/conda.sh
    conda activate anvio-dev
    export PS1="\[\e[0m\e[40m\e[1;30m\] :: anvi'o v7 dev :: \[\e[0m\e[0m \[\e[1;34m\]\]\w\[\e[m\] \[\e[1;31m\]>>>\[\e[m\] \[\e[0m\]"
}

alias anvio-7=init_anvio_7
alias anvio-dev=init_anvio_dev
```

You can either open a new terminal window or run `source ~/.bash_profile` to make sure these changes take effect. Now you should be able to type `anvio-7` to initialize the stable anvi'o, and `anvio-dev` to initialize the development branch of the codebase.

Here is what I see in my terminal for `anvio-7`:

```
meren ~ $ anvi-self-test -v
-bash: anvi-self-test: command not found

meren ~ $ anvio-7

:: anvi'o v7 :: ~ >>>

:: anvi'o v7 :: ~ >>> anvi-self-test -v
Anvi'o .......................................: hope (v7)

Profile database .............................: 35
Contigs database .............................: 20
Pan database .................................: 14
Genome data storage ..........................: 7
Auxiliary data storage .......................: 2
Structure database ...........................: 2
Metabolic modules database ...................: 2
tRNA-seq database ............................: 1
```

Or for `anvio-dev`:

```
meren ~ $ anvi-self-test -v
-bash: anvi-self-test: command not found

:: anvi'o v7 :: ~ >>> anvio-dev

:: anvi'o v7 dev :: ~ >>>

:: anvi'o v7 dev :: ~ >>> anvi-self-test -v
Anvi'o .......................................: hope (v7-dev)

Profile database .............................: 35
Contigs database .............................: 20
Pan database .................................: 14
Genome data storage ..........................: 7
Auxiliary data storage .......................: 2
Structure database ...........................: 2
Metabolic modules database ...................: 2
tRNA-seq database ............................: 1
```

**But please note** that both aliases run `deactivate` and `conda deactivate` first, and they may not work for you especially if you have a fancy setup.


## Other installation options

You will always find the official archives of anvi'o code as at the bottom of our GitHub releases as `anvio-X.tar.gz`:

[https://github.com/merenlab/anvio/releases/latest](https://github.com/merenlab/anvio/releases/latest)

The best way to see what additional software you will need running on your computer for anvi'o to be happy is to take a look at the contents of [this conda recipe](https://github.com/merenlab/anvio/blob/master/conda-recipe/anvio/meta.yaml) (which is a conda build recipe, but it will give you the idea (ignore anvio-minimal, you basically have that one taken care of when you have anvi'o installed)).

Don't be a stranger, and let us know if you need help.

{% include _join-anvio-slack.html %}

---

{:.notice}
{% include _fixthispage.html source="_posts/anvio/2016-06-26-installation-v2.md" %}

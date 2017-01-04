---
layout: post
title: "Installing third-party software"
excerpt: "Recipes to install various software tools anvi'o uses"
modified: 2016-06-18
categories: [anvio]
comments: true
---

This article describes some recipes to install software written by other groups.

Most software we rely on to enhance anvi'o's abilities *do* allow us to re-distribute their code, or have them pre-installed, however, we do not want to follow that route. Although doing that would have made your life much easier, internalizing third-party software from within other platforms directly makes users unable to appreciate other groups' efforts.

As an apology, we will do our best to keep this article up-to-date, so installing third-party software anvi'o uses will not be a big hassle for you. Thank you for your understanding, and your patience in advance.

{% include _toc.html %}

{:.notice}
{% include _fixthispage.html source="anvio/2016-06-18-installing-third-party-software.md" %}

## Prodigal

[Prodigal](http://prodigal.ornl.gov/) is a bacterial and archaeal gene finding program developed at Oak Ridge National Laboratory and the University of Tennessee. Everytime you create a contigs database in anvi'o with `anvi-gen-profile-database`, you use it.

**Citation**: [http://www.biomedcentral.com/1471-2105/11/119](http://www.biomedcentral.com/1471-2105/11/119)

Go to your terminal, and type `prodigal -v` if you get an error, you need to install it, __if the version number is smaller than 2.6.2__, you need to update it.

Here is how to install v2.6.2 (the first line will not work if you don't have wget, but you can get wget installed esily typing `sudo port install wget` if you are using MacPorts system on your Mac computer):

{% highlight bash %}
wget https://github.com/hyattpd/Prodigal/archive/v2.6.2.tar.gz
tar -zxvf v2.6.2.tar.gz && cd Prodigal-2.6.2/ && make
sudo cp prodigal /usr/local/bin/
{% endhighlight %}

Type `prodigal -v` again to make sure everything is alright, and you get the proper version number.

## HMMER

[HMMER](http://hmmer.org/) uses hidden Markov models to perform sequence search and alignments. Everytime you run `anvi-run-hmmss` program, you use it.

**Citation**: [http://hmmer.org/](http://hmmer.org/)

Go to your terminal, and type `hmmscan -h`, if you get an error, you need to install HMMER, if the version number is less than 3.1, you need to update it.

Here is how to install v3.1b2:

{% highlight bash %}
wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2.tar.gz
tar -zxvf hmmer-3.1b2.tar.gz
cd hmmer-3.1b2
./configure && make && sudo make install
cd easel && make check && sudo make install
{% endhighlight %}

Type `hmmscan -h` again to make sure everything is alright, and you get the proper version number.

## SQLite

[SQLite](http://www.tutorialspoint.com/sqlite/sqlite_installation.htm) is a software library that implements a self-contained, serverless, zero-configuration, transactional SQL database engine. Anvi'o uses SQLite pretty much all the time.

**Citation**: [https://www.sqlite.org/](https://www.sqlite.org/)

Go to your terminal, type `sqlite3 --version`, if you get an error, you need to install it. Extensive installation instructions are [here](http://www.tutorialspoint.com/sqlite/sqlite_installation.htm). Or you can install it by typing `sudo port install sqlite3` if you are using the port system on a Mac OSX computer.

{:.notice}
Note: Although this is completely optional, you may also want to consider installing [DB Browser for SQLite](http://sqlitebrowser.org/). It is a lightweight, open-source database browser a nice graphical interface that is very easy-to-install. You probably will never need it or use it, but it may be handy at some point.


## GNU Scientific Library

[GSL](http://www.gnu.org/software/gsl/) is a widely used C library for scientific computation. The only thing depends on GSL is the CONCOCT extension in the codebase. The installation is quite straightforward on most systems. If you are using MacPorts, you can type this on your terminal: `port install gsl gsl-devel py27-gsl` (Rika tells me homebrew on Mac works, too). Otherwise, try these commands and you should be OK:

{% highlight bash %}
wget ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz
tar -zxvf gsl-latest.tar.gz
cd gsl-*
./configure && make && sudo make install
{% endhighlight %}

## NumPY

[NumPY](http://www.numpy.org/) is the fundamental package for scientific computing with Python. Anvi'o uses numpy quite often, and probably not in the best way possible.

**Citation**: [https://arxiv.org/abs/1102.1523](https://arxiv.org/abs/1102.1523)

You don't need to install numpy if you get no complaints back when you type `python -c "import numpy"` in your terminal. If you do get an import error, then you need to install numpy. You can try this:

{% highlight bash %}
sudo pip install numpy
{% endhighlight %}

## Cython

[Cython](http://cython.org/) is "*an optimising static compiler for both the Python programming language and the extended Cython programming language*". If `python -c "import Cython"` in your terminal does not complain, you are golden. Otherwise, install it by running this:

{% highlight bash %}
sudo pip install Cython
{% endhighlight %}

## HDF5

[HDF5](https://www.hdfgroup.org/HDF5/) is "*a data model, library, and file format for storing and managing data*". If you are not sure what it is, you probably don't have it, but we are a big fan of HD5 here in anvi'o development side. If you are using macports on your Mac, you can get away with `sudo port install hdf5`, otherwise you can run these commands on your terminal (these are for version 1.8.16, feel free to check whether there is a newer release of HDF5 from [here](https://www.hdfgroup.org/ftp/HDF5/current/src/), and install the most curent tar.bz2 file in that directory instead):

{% highlight bash %}
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.17/src/hdf5-1.8.17.tar.bz2 
tar -jxvf hdf5-1.8.17.tar.bz2
cd hdf5-1.8.17
./configure && make && sudo make install
{% endhighlight %}

{:.notice}
Depending on your operating system and version, you may need to install `libhdf5-dev` package separately to avoid fatal `No such file or directory` errors for various header files (we heard complaints from Debian and Ubuntu users).
 
## Centrifuge

[Centrifuge](http://www.ccb.jhu.edu/software/centrifuge/) is a "*classification engine that enables rapid, accurate and sensitive labeling of reads and quantification of species on desktop computers*".

**Citation**: [http://biorxiv.org/content/early/2016/05/25/054965](http://biorxiv.org/content/early/2016/05/25/054965)

To install centrifuge, you need to first decide where you want to put all its files on your disk. It could be a directory under `/opt`, or `/usr/local`, or somewhere under your user directory, in case you don't have superuser access on the machine you are working on. Once you know where, open a terminal and set an environment variable to point the base directory you want to keep all centrifuge files:

{% highlight bash %}
$ export CENTRIFUGE_BASE="/path/to/a/directory"
{% endhighlight %}

Do not forget to make sure your version of `/path/to/a/directory` is a full path, and starts with a `/` character.

{:.notice}
More on the "**full path**" thingy: Let's say I want to put all centrifuge related stuff in a directory called `CENTRIFUGE` in my home. Here is what I do: First, in my terminal I type `cd` to makes sure I am in my home directory. Then I type `mkdir -p CENTRIFUGE` to make sure the directory `CENTRIFUGE` exists in my home. Then I type `cd CENTRIFUGE` to go into it. Finally I type `pwd` to get the full path, and replace that entire string with `/path/to/a/directory` in the command above (still keeping it in double quotes) before running the export command.

Then you will get the code, and compile it:

{% highlight bash %}
cd $CENTRIFUGE_BASE
git clone https://github.com/infphilo/centrifuge
cd centrifuge
git checkout 30e3f06ec35bc83e430b49a052f551a1e3edef42
make
{% endhighlight %}

This compiles everything, but does not install anything. To make sure binary files are available directly, you can run this:

{% highlight bash %}
$ export PATH=$PATH:$CENTRIFUGE_BASE/centrifuge
{% endhighlight %}

If everything is alright so far, this is what you should see if you run the following command:

{% highlight bash %}
$ centrifuge --version | head -n 1
centrifuge-class version v1.0.1-beta-27-g30e3f06ec3
{% endhighlight %}

Good? Good. If it does not work, it means you made a mistake with your path variables. If it worked, it means you are golden, and now you should add those two lines in your `~/.bashrc` or `~/.bash_profile` file (whichever one is being used on your system, most likely `~/.bash_profile` will work) to make sure it is set in your environment every time you start a new terminal (clearly with the right full path):

{% highlight bash %}
export CENTRIFUGE_BASE="/path/to/a/directory"
export PATH=$PATH:$CENTRIFUGE_BASE/centrifuge
{% endhighlight %}

You can test whether you managed to do this right by opening a new terminal, and typing `centrifuge --version`. Did it work? Good. Then you set your environment variables right.

Now you have a working centrifuge installation. But not databases to do anything with. For that, you will need to download pre-computed indexes (unless you want to go full Voldemort and compile your own indexes). The compressed indexes for Bacteria, Viruses, Human genome is 6.3 Gb, and it will take about 9 Gb on your disk uncompressed. You will download this data and unpack it only for once:

{% highlight bash %}
$ cd $CENTRIFUGE_BASE
$ wget ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/b+h+v.tar.gz
$ tar -zxvf b+h+v.tar.gz && rm -rf b+h+v.tar.gz
{% endhighlight %}

If everything went alright, you should see something similar to this when you run the following command:

{% highlight bash %}
$ ls -lh $CENTRIFUGE_BASE/b+h+v/*cf
-rw-r--r--   6.5G Feb 15 13:18 $CENTRIFUGE_BASE/b+h+v/b+h+v.1.cf
-rw-r--r--   2.3G Feb 15 13:18 $CENTRIFUGE_BASE/b+h+v/b+h+v.2.cf
-rw-r--r--   1.4M Feb 15 13:18 $CENTRIFUGE_BASE/b+h+v/b+h+v.3.cf
{% endhighlight %}

Good? Good! See? You are totally doing this!

## MCL

[MCL](http://www.micans.org/mcl/index.html?sec_software) is "*a fast and scalable unsupervised cluster algorithm for graphs based on simulation of (stochastic) flow in graphs*", developed by [Stijn van Dongen](http://micans.org/stijn/). If when you type `mcl --version` in your terminal, if you are seeing `mcl 14-137` as an output, you are golden. Otherwise you can install it the following way:

{% highlight bash %}
 $ wget http://www.micans.org/mcl/src/mcl-14-137.tar.gz
 $ tar -zxvf mcl-14-137.tar.gz && cd mcl-14-137
 $ ./configure && make && sudo make install
{% endhighlight %}

Once you are done, you should get a simple usage statement instead of a command not found error when you type `mcl` in your terminal. If that is the case, you are done.


## egnogg-mapper

[eggnog-mapper](https://github.com/jhcepas/eggnog-mapper) _is a tool for fast functional annotation of novel sequences (genes or proteins) using precomputed eggNOG-based orthology assignments_.

**Citation**: [http://biorxiv.org/content/early/2016/09/22/076331](http://biorxiv.org/content/early/2016/09/22/076331)

The official codebase for `eggnog-mapper` is [here](https://github.com/jhcepas/eggnog-mapper), and a pre-print by [Jaime Huerta-Cepas](http://big.crg.cat/people/jaime_huerta_cepas) and his colleagues describing the work is [here](http://biorxiv.org/content/early/2016/09/22/076331). If you follow this recipe, you should remember that you will be using `eggNOG` databases with `eggnog-mapper`, and in your writings you should cite the `eggNOG` release, too:

**Citation**: [https://www.ncbi.nlm.nih.gov/pubmed/26582926](https://www.ncbi.nlm.nih.gov/pubmed/26582926)

{:.notice}
`eggnogg-mapper` has online [documentation](https://github.com/jhcepas/eggnog-mapper/wiki) for you to read and set it up on your system yourself, and learn about the details of working with it. This is a recipe for the lazy. If you have a systems administrator, it may be better for them to set it up as a module for everyone. Otherwise, this recipe will tell you how you can you do it within your own space (note that you will need lots of disk space depending on databases you want to download).

To install `eggnog-mapper` you first need to get the source code, and then you will need to collect the precomputed database files. 

First, you need to decide where do you want to put `eggnog-mapper` and its databases (you will need to change that `/path/to/a/directory` line to wherever you want on your disk):

{% highlight bash %}
$ export EGGNOG_MAPPER_BASE="/path/to/a/directory"
{% endhighlight %}

Here is how you get the code:

{% highlight bash %}
$ cd $EGGNOG_MAPPER_BASE
$ git clone https://github.com/jhcepas/eggnog-mapper.git
$ cd eggnog-mapper/
$ git checkout tags/0.12.6
$ export PATH=$PATH:$EGGNOG_MAPPER_BASE/eggnog-mapper
{% endhighlight %}

At this point if you run this command, you should get the following output:

{% highlight bash %}
$ emapper.py --version
emapper-0.12.6
{% endhighlight %}

If all is good, now you can download the databases. Which databases you are going to be downloading is up to you (which will not only affect the disk space you need, but also the runtime to screen your genes). Here I will download everything (because I have time _and_ space):

{% highlight bash %}
$ download_eggnog_data.py euk bact arch viruses -y
{% endhighlight %}

This will take a long *very* long time mostly due to large I/O overhead to decompress some of the databases with large numbers of smaller files (so do not forget to start the process in a `screen`), but fortunately you will not do it again.

If you are here, you have the basic setup done. Congratulations.

As a very final step, you should add these two lines in your `~/.bashrc` or `~/.bash_profile` file (whichever one is being used on your system, most likely `~/.bash_profile` will work) to make sure they are set in your environment every time you start a new terminal (don't forget to update the directory name):

{% highlight bash %}
export EGGNOG_MAPPER_BASE="/path/to/a/directory"
export PATH=$PATH:$EGGNOG_MAPPER_BASE/eggnog-mapper
{% endhighlight %}

## muscle

[muscle](http://drive5.com/muscle/) _is one of the best-performing multiple alignment programs_.

**Citation**: [https://www.ncbi.nlm.nih.gov/pubmed/15034147](https://www.ncbi.nlm.nih.gov/pubmed/15034147)

Anvi'o uses muscle to align amino acid sequences within each protein cluster while running the pangenomic workflow in `v2.1.0` and later versions of anvi'o. Installation is rather easy: go to the [downloads page](http://www.drive5.com/muscle/downloads.htm) for `muscle`, grab the one that matches to your operating system, rename the unzipped binary to 'muscle', and move it into `/usr/local/bin` or whichever directory seems to be working.

If you were successful, this is what you should see when you type `muscle` in your terminal:

{% highlight bash %}
$ muscle -version
MUSCLE v3.8.31 by Robert C. Edgar
{% endhighlight %}

---
layout: post
title: "Installing anvi'o"
excerpt: "Instructions to install the platform."
modified: 2015-05-01 
tags: []
categories: [anvio]
comments: true
---

{% include _toc.html %}

This article describes basics steps of installing Anvi'o. If you run into any issues, please post a comment down below, or open an <a href="https://github.com/meren/anvio/issues">issue</a>.

Here I would like to thank Inés Martínez, Rika Anderson, and Sharon Grim for helping me a lot by volunteering themselves to test the installation.

{: .notice}
Are you using Mac OS X? Then you really should see this one first: <a href="../../../../2015/08/20/installation-on-mac/">OS X installer for anvi'o</a>

## Dependencies

anvi'o has some dependencies, some of which will be taken care of the installer. However, you first need to make sure your system does have all the following software (it is not as scary as it looks, you will be fine):

* [Chrome Web Browser](https://www.google.com/chrome/browser/desktop/). Chrome should be not only installed on your system, but unfortunately it should also be the default browser (otherwise everytime interactive interface pops up, you will need to copy-paste the address to a Chrome window). anvi'o does not support any other browser, and __it will not perform optimally__ on others.

* [Prodigal](http://prodigal.ornl.gov/) (go to your terminal, type `prodigal -v` if you get an error, you need to install it, __if the version number is smaller than 2.6.2__, you need to update it). Here is a quick way to install the 2.6.2 version, but feel free to do visit the [latest releases page](https://github.com/hyattpd/prodigal/releases/) and do it yourself (the first line will not work if you don't have wget, but you can get wget installed esily typing `sudo port install wget` if you are using MacPorts system on your Mac computer):

<div style="padding-left:30px">
<pre>
wget https://github.com/hyattpd/Prodigal/archive/v2.6.2.tar.gz
tar -zxvf v2.6.2.tar.gz && cd Prodigal-2.6.2/ && make
sudo cp prodigal /usr/bin/
</pre>
</div>

* [HMMer](http://hmmer.janelia.org/) (go to your terminal, type `hmmscan`, if you get an error, you need to install this one). It is quite easy to install. Go to the downloads page, [http://hmmer.janelia.org/download.html](http://hmmer.janelia.org/download.html), download the *source* for the latest version, unpack it, go into that directory with your terminal, and then type (don't forget to check whether it is installed properly by typing `hmmscan` again after installation):

<div style="padding-left:30px">
<pre>
./configure && make && sudo make install
</pre>
</div>

* [SQLite](http://www.tutorialspoint.com/sqlite/sqlite_installation.htm) (go to your terminal, type sqlite3, if it gives you an error, you need to install this one). There is installation instructions on the web page if you follow the link. Or you can install it by typing `sudo port install sqlite3` if you are using the port system on your Mac.

* [MyRAST](http://blog.theseed.org/servers/) (this one is optional, but it is very useful, if you type `svr_assign_to_dna_using_figfams` and if it doesn't give you an error, you have it (press `CTRL+C` to quit)).

* [GSL](http://www.gnu.org/software/gsl/), which is the GNU Scientific Library. It is quite straightforward. If you are using MacPorts, you can install `gsl`, `gsl-devel`, and `py27-gsl` packages from the terminal using `port install` (Rika tells me homebrew on Mac works, too). Otherwise, try these commands and you should be OK:

<div style="padding-left:30px">
<pre>
wget ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz
tar -zxvf gsl-latest.tar.gz
cd gsl-*
./configure && make && sudo make install
</pre>
</div>

* [numpy](http://www.numpy.org/). You don't need to install numpy if you get no complaints back when you type `python -c "import numpy"` in your terminal. If you do get an import error, then you need to install numpy. You can try this:

<div style="padding-left:30px">
<pre>
sudo pip install numpy
</pre>
</div>

* [Cython](http://cython.org/). Similar to the numpy case. Test whether you are good to go by typing `python -c "import Cython"` in your terminal, and install it by running this if there is an import error:

<div style="padding-left:30px">
<pre>
sudo pip install numpy
</pre>
</div>

### Other installation notes

* Anvi'o uses SQLite to create self-contained databases to store information. There are many bindings for many programming languages to access to these database files and explore them, and it is also possible to use `sqlite3` program from the terminal to play witht them. [DB Browser for SQLite](http://sqlitebrowser.org/) is a very easy-to-install open-source software that does what `sqlite3` is doing with a nice graphical interface. I urge you to stick with the terminal as much as possible, but in case you insist on using a GUI, this is the one you should give a chance.

## Installation

You can either install a stable release of anvi'o, or you can get a copy of the latest snapshot from the repository (clearly, it is always safer to go with the stable release).

### Installing the latest stable release (safe mode)

The best way to install stable release is to do it through pip. If you are have your dependencies sorted, tyr running this on your terminal:

    sudo pip install anvio

If there are no errors, you are golden.

If you want to make sure everything is working, you can run the "mini test" (explained down below).


### Installing or updating from the current codebase (semi-pro)

If this is your first time with the codebase, get a fresh copy (with all the submodules necessary):

    git clone --recursive https://github.com/meren/anvio.git

Then go into the anvio directory, and then run the installation:

    cd anvio
    sudo python setup.py install

If you already have the codebase, and if your purpose is to _update_ your already existing installation, you are going to need to run these commands from within the anvio directory instead of the ones above:

    git pull
    git submodule update --init --recursive
    sudo python setup.py install

No errors? Perfect. Run the mini test!

### Installation for developers (you're a wizard, ary)

If you are planning to do this you need no introduction, but I will give one anyway. Clone the codebase into a `$DIR` you like:


    cd $DIR
    git clone --recursive https://github.com/meren/anvio.git

Then edit your `~/.bashrc` or `~/.bash_profile` files depending on your system configuration to update your `PATH` and `PYTHONPATH` environment variables:

    export PYTHONPATH=$PYTHONPATH:$DIR/anvio/
    export PATH=$PATH:$DIR/anvio/bin:$DIR/anvio/sandbox

Because you didn't run "build", your C extensions will not be compiled. This is the fast way to get them compiled and put in the right place:

    cd $DIR/anvio
    python setup.py build
    cp build/lib.*/anvio/*so anvio/

That's it.

After sourcing your `.bashrc` (or `.bash_profile`) and get the new environment variables set, you should be able run the 'mini test'.

Now you can edit the codebase, and test it, without re-installing anvi'o over and over again.

## Running the "Mini Test"

"Mini test" is [a tiny shell script](https://github.com/meren/anvio/blob/master/tests/run_mini_test.sh) that runs almost everything implemented in anvi'o on a very small dtaset.

If you have a proper installation, you shouldn't get any errors when you run it.

If you installed anvi'o via pip (following the safe mode), you can go to an empty directory, and type these:

    git clone --recursive https://github.com/meren/anvio.git
    cd anvio/tests
    ./run_mini_test.sh

If you installed anvi'o from the codebase, you can go to the anvi'o source code directory, and simply type these:

    cd tests
    ./run_mini_test.sh

Upon the successful completion of the mini test run, your browser should popup to take you to the interactive interface. When you click that 'Draw' button, you should see something like this (this is the old version of the anvi'o interface, and it shall stay here so we remember where we came from):

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/misc/mini-test-screenshot.png"><img src="{{ site.url }}/images/anvio/misc/mini-test-screenshot.png" width="50%" /></a>
</div>

All fine? Perfect! You have a running installation of anvi'o.

Maybe you would like to continue with [other posts on the platform]({{ site.url }}/projects/anvio).
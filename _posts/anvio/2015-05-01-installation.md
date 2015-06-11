---
layout: post
title: "Installing Anvi'o"
excerpt: "Instructions to install the platform."
modified: 2015-05-01 
tags: []
categories: [anvio]
comments: true
---

{% include _toc.html %}

This article describes basics steps of installing Anvi'o. If you run into any issues, please post a comment down below, or open an <a href="https://github.com/meren/anvio/issues">issue</a>.

Here I would like to thank Inés Martínez, Rika Anderson, and Sharon Grim for helping me a lot by volunteering themselves to test the installation.

## Dependencies

anvi'o has some dependencies, some of which will be taken care of the installer. However, you first need to make sure your system does have the following software (it is not as scary as it looks, you will be fine):

* [Chrome Web Browser](https://www.google.com/chrome/browser/desktop/). Chrome should be not only installed on your system, but unfortunately it should also be the default browser (otherwise everytime interactive interface pops up, you will need to copy-paste the address to a Chrome window). anvi'o does not support any other browser, and it will not perform optimally on others.

* [Prodigal](http://prodigal.ornl.gov/) (go to your terminal, type `prodigal` if you get an error, you need to install it). Here is a quick way to install it (the first line will not work if you don't have wget, but you can get wget installed esily typing `sudo port install wget` if you are using MacPorts system on your Mac computer):

<div style="padding-left:30px">
<pre>
wget http://prodigal.ornl.gov/zips/prodigal.v2_50.tar.gz
tar -zxvf prodigal.v2_50.tar.gz && cd prodigal.v2_50 && make
sudo cp prodigal /usr/bin/
</pre>
</div>

* [HMMer](http://hmmer.janelia.org/) (go to your terminal, type `hmmscan`, if you get an error, you need to install this one). There is an easy installer there if you follow the link. Alternatively, if you are using the port system on your Mac, you can simply type `sudo port install hmmer`.

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

### Other installation notes

* Please make sure your pysam version is 0.7.7, and nothing else. Although the installer is supposed to make sure you have the right version, more than once people run into trouble with pysam. You can learn your pysam version by typing `python` in your terminal, and then copy-pasting this line: `import pysam; pysam.__version__`. Does it say `0.7.7`? If it does, you are OK. If it doesn't, you need to install the right version. You can run these two commands to do that (after the intallation don't forget to chek whether you see the right version in your python terminal):

<div style="padding-left:30px">
<pre>
sudo pip uninstall pysam
pip install pysam==0.7.7
</pre>
</div>

* Anvi'o uses SQLite to create self-contained databases to store information. There are many bindings for many programming languages to access to these database files and explore them, and it is also possible to use `sqlite3` program from the terminal to play witht them. [DB Browser for SQLite](http://sqlitebrowser.org/) is a very easy-to-install open-source software that does what `sqlite3` is doing with a nice graphical interface. I urge you to stick with the terminal as much as possible, but in case you insist on using a GUI, this is the one you should give a chance.

## Installation

You can either install a stable release of anvi'o, or you can get a copy of the latest snapshot from the repository (it is always safer to with the stable release).

### Installing the latest stable release

Go here, and download the latest release:

<p style="padding-left: 30px"><a href="http://merenlab.org/projects/anvio/downloads/" target="_blank">http://merenlab.org/projects/anvio/downloads/</a></p>

Install it by typing these commands:

    tar -zxvf anvio-0.8.5.tar.gz
    cd anvio-0.8.5
    sudo python setup.py install

If there are no errors, you are golden. Do not forget to run the mini test.


### Installing or updating from the current codebase

If this is your first time with the codebase, get a fresh copy:

    git clone https://github.com/meren/anvio.git

Then go into the anvio directory, and then run the installation:

    cd anvio
    sudo python setup.py install

If you already have the codebase, and if your purpose is to _update_ your already existing installation, you are going to need to run these commands from within the anvio directory instead of the ones above:

    git pull
    sudo python setup.py install

No errors? Perfect. Run the mini test!

### Pro installation for developers

If you are planning to do this you need no introduction, but I will give one anyway. Clone the codebase into a `$DIR` you like:


    cd $DIR
    git clone https://github.com/meren/anvio.git

Then edit your `~/.bashrc` or `~/.bash_profile` files depending on your system configuration to update your `PATH` and `PYTHONPATH` environment variables:

    export PYTHONPATH=$PYTHONPATH:$DIR/anvio/
    export PATH=$PATH:$DIR/anvio/bin:$DIR/anvio/sandbox

This is not enough, though. Because you didn't run build, your C extensions will not be compiled. Here is what you will do put them in the right place:

    cd $DIR/anvio
    python setup.py build
    cp build/lib.*/anvio/*so anvio/

That's it. After sourcing your .bashrc (or .bash_profile) and get the new environment variables set, you should be able run the 'mini test'. Now you can edit the codebase without re-installing anvi'o over and over again.

## Running the "Mini Test"

"Mini test" is a minimum test set that runs almost everything in the codebase. If you have a proper installation, you shouldn't get any errors from running this test. To run it go to your anvio installation directory, and type these:

    cd tests
    ./run_mini_test.sh

If you are not on a server where there is no access to a GUI, your browser should popup upon the successful completion of `run_mini_test.sh`, and when you click that 'Draw' button, you should see this:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/misc/mini-test-screenshot.png"><img src="{{ site.url }}/images/anvio/misc/mini-test-screenshot.png" width="50%" /></a>
</div>

All fine? Perfect!


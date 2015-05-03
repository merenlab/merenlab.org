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

### Dependencies

anvi'o has some dependencies, some of which will be taken care of the installer. There are two you need to make sure you have installed:

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

### Installation

First get the codebase:

    git clone https://github.com/meren/anvio.git

Then go into the anvio directory, and then run the installation:

    cd anvio
    sudo python setup.py install

If you want to update your installation, you are going to need to run these commands from within the anvio directory:

    git pull
    sudo python setup.py install


### Running the "Mini Test"

"Mini test" is a minimum test set that runs almost everything in the codebase. If you have a proper installation, you shouldn't get any errors from running this test. To run it go to your anvio installation directory, and type these:

    cd tests
    ./run_mini_test.sh

All fine? Perfect!


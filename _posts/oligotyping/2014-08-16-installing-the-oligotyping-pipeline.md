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


## Installing the latest stable version using pip (suggested method)

Install [macports](https://www.macports.org/) if you haven’t. You can download it from [here](https://www.macports.org/install.php) (this is the only installer you will download). Macports provides an easy access to thousands of open-source software. Once you have it installed, open a terminal and type this:

    sudo port selfupdate

Next, type this command:

    sudo port install python27

This will install Python 2.7 on your system. I know that you already have Python installed, but you don’t want to use that one.

Next, type this command:

    sudo port select --set python python27

This will make sure you use the correct Python version from the correct place. When you type the following command, you should see a response that starts with `/opt/local/...`:

    which python

Next, type these commands:

    sudo port install py27-pip py27-scipy py27-matplotlib py27-biopython py27-django git
    sudo port select --set pip pip27

The first one is to install some dependencies. The second one is to make sure pip27 is the default, and you can simply access to it by typing “pip” instead of “pip27″.

If you don’t have R installed, type this (if you have an older version of R installed, you may have issues with the next step, in that case please consider doing what this step suggests):

    sudo port install R

In the terminal, type R, and enter these lines (if everything goes all right with the first line, you are going to need to press ‘n’ after the second line):

    install.packages(c('vegan', 'ggplot2', 'gplots', 'gtools', 'reshape', 'optparse', 'pheatmap', 'RColorBrewer', 'compute.es'))
    quit()

Install NCBI+ executables (especially blastn (v 2.2.*), they are avilable from this address (the easiest way to install NCBI tools for MAC users is to download the .dmg file):

[ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

After the installation, you need to make sure when you open a terminal and type `blastn -help`, you get the help menu, instead of a command not found error.

If you came this far, it means you have installed everything necessary. Finally, you can run this command to install the oligotyping pipeline:

    sudo pip install oligotyping

When you type this command you should get a help menu, instead of a command not found error:

    oligotype --help

**Note for MAC users:** It seems running this command sometimes results in a `command not found` error, although all files are in place. In that case you may need to update your `$PATH` variable (so the shell knows where to look for Python programs). To do that, first make sure that you see the `oligotype` binary in this output:

    ls /opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin/oligotype

If it is there, update your `PATH` variable permanently by running these commands:

{% gist a9bc1bf856adb3525a95 %}


## Installing the latest stable version from the source code (OK method)


Pip installation makes things easier, but alternatively you can download the source package of the stable version and install it manually. For this, you still need to complete the first 7 steps in the previous recipe, and then instead of running 8th step, you can follow these steps.

Download the source code from this address:

[http://oligotyping.org/files/downloads/](http://oligotyping.org/files/downloads/)

You should download the latest version available, but I will assume you downloaded `oligotyping-1.0.tar.gz` in this example.

Go to the directory where you downloaded this file in terminal, and type these commands:

    tar -zxvf oligotyping-1.0.tar.gz
    cd oligotyping-1.0
    sudo python setup.py install

If everything went alright, when you type this command you should get a help menu, instead of a command not found error:

    oligotype --help


## Installing the latest snapshot from the GitHub repository (Hacker mode [ON])

The repository contains the latest codebase. It may be unstable, but it may also have features that haven’t appeared in a stable version of the oligotyping pipeline. This is what you need to do to install oligotyping from the GitHub repository:

From the terminal, type this:

    git clone git://github.com/meren/oligotyping.git

This will create a directory called ‘oligotyping’.

Then type these commands to install the pipeline:

    cd oligotyping
    sudo python setup.py install


## Using the oligtyping pipeline from directly the master repo without installation (Meren's way)

This is for command line gurus. See the next section if you don't feel comfortable with this one. Keep in mind that you will need to change the directory names to adapt this recipe for your system:

``` bash
# setup a virtual environment with Python 2.7, and activate it
mkdir -p ~/virtual-envs/
virtualenv -p python2.7 ~/virtual-envs/oligotyping-master
source ~/virtual-envs/oligotyping-master/bin/activate

# get the source code, and install requirements
mkdir -p ~/github/
cd ~/github/
git clone git@github.com:merenlab/oligotyping.git
cd oligotyping
pip install -p requirements.txt
deactivate

# setup init scripts for easy use
echo 'export PYTHONPATH=$PYTHONPATH:~/github/oligotyping' >> ~/virtual-envs/oligotyping-master/bin/activate
echo 'export PATH=$PATH:~/github/oligotyping/bin' >> ~/virtual-envs/oligotyping-master/bin/activate
echo 'alias oligotyping-activate-master="source ~/virtual-envs/oligotyping-master/bin/activate"' >> ~/.bash_profile
echo 'cd ~/github/oligotyping/ && git pull && cd -' >> ~/virtual-envs/oligotyping-master/bin/activate
source ~/.bash_profile

# setup matplotlib so it works from within the virtualenv
mkdir -p ~/.matplotlib
echo 'backend: Agg' >> ~/.matplotlib/matplotlibrc
```

Now you can activate the oligotyping pipeline with the latest additions by typing this in your command line:

``` bash
oligotyping-activate-master
```

You will still need to setup other dependencies such as R, for which you can see the instructions in the 'suggested method' section.

---


I thank Les Dethlefsen very much for sharing his experience with the installation and helping me improve the document.

Please don’t hesitate to ask questions about the installation at the [here](https://groups.google.com/forum/#!forum/oligotyping).



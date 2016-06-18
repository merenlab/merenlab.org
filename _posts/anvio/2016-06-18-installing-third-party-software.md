---
layout: post
title: "Installing third-party software"
excerpt: "Recipes to install varius software tools anvi'o uses"
modified: 2016-06-18
categories: [anvio]
comments: true
---

This article describes some recipes to install software written by other groups.

Most software we rely on to enhance anvi'o's abilities *do* allow us to re-distribute their code, or have them pre-installed, however, we do not want to follow that route. Although doing that would have made your life much easier, internalizing third-party software from within other platforms directly makes users unable to appreciate other groups' efforts.

As an apology, we will do our best to keep this article up-to-date, so installing third-party software anvi'o uses will not be a big hassle for you. Thank you for your understanding, and your patience in advance.

{% include _toc.html %}

## Centrifuge

[Centrifuge](http://www.ccb.jhu.edu/software/centrifuge/) is a "*classification engine that enables rapid, accurate and sensitive labeling of reads and quantification of species on desktop computers*".

**Citation**: [http://biorxiv.org/content/early/2016/05/25/054965](http://biorxiv.org/content/early/2016/05/25/054965)

To install centrifuge, you need to first decide where you want to put all its files on your disk. It could be a directory under `/opt`, or `/usr/local`, or somewhere under your user directory, in case you don't have superuser access on the machine you are working on. Once you know where, open a terminal and set an environment variable to point the base directory you want to keep all centrifuge files:

{% highlight bash %}
$ export CENTRIFUGE_BASE="/path/to/a/directory"
{% endhighlight %}

Do not forget to make sure your version of `/path/to/a/directory` is a full path, and starts with a `/` character.

{:.notice}
More on the "**full path**" thingy: Let's say I want to put all centrifuge related stuff in a directory called `CENTRIFUGE` in my home. Here is what I do: First, in my terminal I type `cd` to makes ure I am in my home directory. Then I type `mkdir -p CENTRIFUGE` to make sure the directory `CENTRIFUGE` exists in my home. Then I type `cd CENTRIFUGE` to go into it. Finally I type `pwd` to get the full path, and replace that entire string with `/path/to/a/directory` in the command above (still keeping it in double quotes) before running the export command.

Then you will get the code, and compile it:

{% highlight bash %}
$ cd $CENTRIFUGE_BASE
$ git clone https://github.com/infphilo/centrifuge
$ cd centrifuge && make
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

Good? Good. If it does not work, it means you made a mistake with your path variables. If it worked, it means you are golden, and now you should add those two lines in your `~/.bashrc` or `~/.bash_profile` file (whichever one is being used on your system, most likely `~/.bash_profile` will work) to make sure it is set in your environment everytime you start a new terminal (clearly with the right full path):

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

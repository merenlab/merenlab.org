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

[Centrifuge](http://www.ccb.jhu.edu/software/centrifuge/) is a classification engine that enables rapid, accurate and sensitive labeling of reads and quantification of species on desktop computers.

**Citation**: [http://biorxiv.org/content/early/2016/05/25/054965](http://biorxiv.org/content/early/2016/05/25/054965)

To install centrifuge, you need to first decide where you want to put all its files on your disk. It could be a directory under `/opt`, or `/usr/local`, or somewhere under your user directory if you don't have superuser access on the machine you are working on. Once you know where to keep it, open a terminal and start with this, after replacing the path with the correct one:

{% highlight bash %}
$ export CENTRIFUGE="/path/to/a/directory"
{% endhighlight %}

{:.notice}
Do not forget to make sure your version of `/path/to/a/directory` is a full path, and starts with a `/` character.

You *should* put this line in your `~/.bashrc` or `~/.bash_profile` (whichever one is being used on your system) to make sure it is set in your environment everytime you start a new terminal.

Then you will get the code, and compile it:

{% highlight bash %}
$ cd $CENTRIFUGE
$ git clone https://github.com/infphilo/centrifuge
$ cd centrifuge && make
{% endhighlight %}

This compiles everything, but does not install anything. To make sure binary files are available directly, you can run this (and again, you should add this line into your `~/.bashrc` or `~/.bash_profile` to make sure every new terminal session remembers where they are):

{% highlight bash %}
$ export PATH=$PATH:$CENTRIFUGE/centrifuge
{% endhighlight %}

If everything is alright so far, this is what you should see if you run the following command:

{% highlight bash %}
$ centrifuge --version | head -n 1
centrifuge-class version v1.0.1-beta-27-g30e3f06ec3
{% endhighlight %}

Good? Good.

Next, you will need to download pre-computed indexes (unless you want to go Voldemort and compile your own indexes). The compressed indexes for Bacteria, Viruses, Human genome is 6.3 Gb, and it will take about 9 Gb on your disk uncompressed. You will download these only for once:

{% highlight bash %}
$ cd $CENTRIFUGE
$ wget ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/b+h+v.tar.gz
$ tar -zxvf b+h+v.tar.gz && rm -rf b+h+v.tar.gz
$ centrifuge-inspect --name-table $CENTRIFUGE/b+h+v/b+h+v > b+h+v-names-table.txt
{% endhighlight %}

The following command will make sure you have the necessary conversion table ready for anvi'o when the time comes. This is also something you will have to do only once:

{% highlight bash %}
$ centrifuge-inspect --name-table $CENTRIFUGE/b+h+v/b+h+v > $CENTRIFUGE/b+h+v-names-table.txt
{% endhighlight %}

If everything went alright, you should see something like this when you run the following command:

{% highlight bash %}
$ wc -l $CENTRIFUGE/b+h+v-names-table.txt
   13455 $CENTRIFUGE/b+h+v-names-table.txt
{% endhighlight %}

See? You are totally doing this!


---
layout: post
title: "Installing anvi'o"
excerpt: "Instructions to install the v2 branch of the platform."
modified: 2016-06-26 
tags: []
categories: [anvio]
redirect_from: /2015/05/01/installation/
comments: true
---

{% include _toc.html %}

This article explains basic steps of installing anvi'o using rather conventional methods.

{:.notice}
There are other installation options for [OS X]({% post_url anvio/2015-08-20-installation-on-mac %}), [Ubuntu Linux](https://gist.github.com/meren/f1e674f774a30a209d82), or [Arch Linux](https://gist.github.com/meren/f2f247f78a6402e3370ea929de1d3510) users, and [Docker]({% post_url anvio/2015-08-22-docker-image-for-anvio %}) enthusiasts.

Please post a comment down below if you have any questions about the installation. You may want to consider opening an <a href="https://github.com/meren/anvio/issues">issue</a> for more technical problems.


## An overview of anvi'o dependencies

You need to make sure your system does have all the following software if you are going to follow the installation instructions on this page. It is not as scary as it looks. If you just follow the simple installation instructions behind these links for each of them, you will most probably be fine:

* [Prodigal]({% post_url anvio/2016-06-18-installing-third-party-software %}#prodigal){:target="_blank"}
* [HMMER]({% post_url anvio/2016-06-18-installing-third-party-software %}#hmmer){:target="_blank"}
* [SQLite]({% post_url anvio/2016-06-18-installing-third-party-software %}#sqlite){:target="_blank"}
* [GSL]({% post_url anvio/2016-06-18-installing-third-party-software %}#gnu-scientific-library){:target="_blank"}
* [NumPY]({% post_url anvio/2016-06-18-installing-third-party-software %}#numpy){:target="_blank"}
* [Cython]({% post_url anvio/2016-06-18-installing-third-party-software %}#cython){:target="_blank"}
* [HDF5]({% post_url anvio/2016-06-18-installing-third-party-software %}#hdf5){:target="_blank"}

Here is a very final note about the [Chrome Web Browser](https://www.google.com/chrome/browser/desktop/).  Currenlty, Chrome has the most efficient SVG engine among all browsers. For instance, Safari can run the anvi'o interactive interface, however it takes orders of magnitude more time and memory compared to Chrome. Firefox, on the other hand, doesn't even bother drawing anything at all. Long story short, the anvi'o interactive interface __will not perform optimally__ with anything but Chrome. So you need Chrome. Moreover, if Chrome is not your default browser, everytime interactive interface pops up, you will need to copy-paste the address bar into a Chrome window. You can learn what is your default browser by running this command in your terminal:

<div style="padding-left:50px">
{% highlight bash %}
python -c 'import webbrowser as w; w.open_new("http://")'
{% endhighlight %}
</div>

## Installation

You can either install a stable release of anvi'o, or you can get a copy of the latest snapshot from the repository.

### Installing the latest stable release (safe mode)

The best way to install stable release is to do it through `pip`. If you have your dependencies sorted, try running this on your terminal:

{% highlight bash %}
sudo pip install anvio
{% endhighlight %}

Pleaes keep an eye on the output to make sure there are no errors.

If the installation finishes properly, run this in your terminal:

{% highlight bash %}
anvi-profile --version
{% endhighlight %}

If everything looks alright, you are golden. Now you can run `anvi-self-test` to make sure everything else works as expected.

If you get a '*command not found*' error despite a successful installation, this would mean that the directory where your anvi'o programs live is not listed in your `PATH` variable. Please [follow this post]({{ site_url }}/software/anvio/issues/missing-binaries.html) to fix that first.


### Installing or updating from the current codebase (semi-pro)

If this is your first time with the codebase, get a fresh copy (with all the submodules necessary):

{% highlight bash %}
git clone --recursive https://github.com/meren/anvio.git
{% endhighlight %}

Then go into the anvio directory, and then run the installation:

{% highlight bash %}
cd anvio
sudo python setup.py install
{% endhighlight %}

If you already have the codebase, and if your purpose is to _update_ your already existing installation, you are going to need to run these commands from within the anvio directory instead of the ones above:

{% highlight bash %}
git pull
git submodule update --init --recursive
sudo python setup.py install
{% endhighlight %}

No errors? Perfect. Run `anvi-self-test`!

### Installation for developers (you're a wizard, arry)

If you are planning to do this you need no introduction, but I will give you one anyway. Clone the codebase into a `$DIR` you like:

{% highlight bash %}
cd $DIR
git clone --recursive https://github.com/meren/anvio.git
{% endhighlight %}

Then edit your `~/.bashrc` or `~/.bash_profile` files depending on your system configuration to update your `PATH` and `PYTHONPATH` environment variables:

{% highlight bash %}
export PYTHONPATH=$PYTHONPATH:$DIR/anvio/
export PATH=$PATH:$DIR/anvio/bin:$DIR/anvio/sandbox
{% endhighlight %}

{:.notice}
Of course this approach will not help you install all dependencies through pip. To avoid further headaches, I suggest running `sudo pip install anvio` once, so all necessary python modules are automatically installed, and then removing the anvio installation by running `sudo pip uninstall anvio`.

Finally, because you didn't actually run the proper `build`, your C extensions will not be compiled. This is the fast way to get them compiled and put in the right place:

{% highlight bash %}
cd $DIR/anvio
python setup.py build
cp build/lib.*/anvio/*so anvio/
{% endhighlight %}

That's it.

After sourcing your `.bashrc` (or `.bash_profile`) and get the new environment variables set, you should be able run `anvi-self-test`.

Now you can edit the codebase, and test it without re-installing anvi'o over and over again. Plus, everytime you run `git pull` you can get the very latest improvements in the `master` repository.

## Running the "Mini Test"

You can make anvi'o test itself by running the program `anvi-self-test`.

Upon the successful completion of all the tests, your browser should popup to take you to the interactive interface. When you click that 'Draw' button, you should see something like this (this is one of the older version of [the anvi'o interactive interface]({% post_url anvio/2016-02-27-the-anvio-interactive-interface %}), and it shall stay here so we remember where we came from):

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/misc/mini-test-screenshot.png"><img src="{{ site.url }}/images/anvio/misc/mini-test-screenshot.png" width="50%" /></a>
</div>

All fine? Perfect! Now you have a running installation of anvi'o!

It is time to go through the [tutorial]({% post_url anvio/2016-06-22-anvio-tutorial-v2 %}), or take a look at [all the other posts on the platform]({{ site.url }}/software/anvio).

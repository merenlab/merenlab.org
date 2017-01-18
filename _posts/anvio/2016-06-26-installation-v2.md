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

{% include _project-anvio-version.html %}

This article explains basic steps of installing anvi'o using rather conventional methods.

Please post a comment down below if you have any questions about the installation. You may want to consider opening an <a href="https://github.com/meren/anvio/issues">issue</a> for more technical problems.

<div class="extra-info" markdown="1">

<span class="extra-info-header">A note on Chrome</span>

Currently, the [Chrome Web Browser](https://www.google.com/chrome/browser/desktop/) has the most efficient SVG engine among all browsers we tested. For instance, Safari can run the anvi'o interactive interface, however it takes orders of magnitude more time and memory compared to Chrome. Firefox, on the other hand, doesn't even bother drawing anything at all. Long story short, the anvi'o interactive interface __will not perform optimally__ with anything but Chrome. So you need Chrome. Moreover, if Chrome is not your default browser, every time interactive interface pops up, you will need to copy-paste the address bar into a Chrome window. You can learn what is your default browser by running this command in your terminal:

``` bash
python -c 'import webbrowser as w; w.open_new("http://")'
```

</div>

## Painless installation with Homebrew

If you are using Mac and have [Homebrew](http://brew.sh/) installed on your computer, all you need to do is to run this to have anvi'o installed on your system, and skip the rest of this page:

``` bash
brew install anvio
```

*Well, the only problem is that it will not work this time. We are on it. Give us a couple of days.*

If you are installing anvi'o on another operating system (or if you an OS X user who does not like Homebrew for some reason) please continue reading.


## Installation (with varying levels of pain)

{:.notice}
We also other installation options for [Arch Linux](https://gist.github.com/meren/f2f247f78a6402e3370ea929de1d3510) users, and [Docker]({% post_url anvio/2015-08-22-docker-image-for-anvio %}) enthusiasts, too.

First things first. You need to make sure your system does have all the following software if you are going to follow any of the following installation instructions. It is not as scary as it looks. If you just follow these links, you will most probably be golden:

* [Prodigal]({% post_url anvio/2016-06-18-installing-third-party-software %}#prodigal){:target="_blank"}
* [HMMER]({% post_url anvio/2016-06-18-installing-third-party-software %}#hmmer){:target="_blank"}
* [SQLite]({% post_url anvio/2016-06-18-installing-third-party-software %}#sqlite){:target="_blank"}
* [GSL]({% post_url anvio/2016-06-18-installing-third-party-software %}#gnu-scientific-library){:target="_blank"}
* [HDF5]({% post_url anvio/2016-06-18-installing-third-party-software %}#hdf5){:target="_blank"}

Finally you will need `virtualenv`. This should work for most:

``` bash
pip install virtualenv
```

If you don't have `pip`, you will need to visit [this web page](https://pip.pypa.io/en/stable/installing/) to get it installed. Don't forget, anvi'o still uses Python 2.7. If you run into any trouble, send an e-mail to [Google Groups for anvi'o](https://groups.google.com/forum/#!forum/anvio).

OK. If you are still here, you may have gone through the most painful part already. Anvi'o developers are very proud of you.


### Installing the latest stable release (safe mode)

This is the best way to install the stable release. You will do everything in a Python virtual environment. If you are not experienced with computer stuff, do not worry. If you have taken care of your dependencies mentioned above, the rest is very simple.

We first need to create a new virtual environment for anvi'o. Since it is easier to keep all virtual environments in one place, I will first create a directory in my home:

``` bash
mkdir ~/virtual-envs/
```

Then create a new virtual environment for anvi'o under that directory, and activate it:

``` bash
virtualenv ~/virtual-envs/anvio
source ~/virtual-envs/anvio/bin/activate
```

Make sure your paths look alright. Yours should look similar to this:

``` bash
(anvio) meren ~ $ which pip
/Users/meren/virtual-envs/anvio/bin/pip
```

Now you can do the installation:

``` bash
pip install numpy
pip install scipy
pip install cython
pip install anvio
```

If all looks good, now you should be able to run `anvi-self-test`:

``` bash
anvi-self-test --suite mini
```

If this runs successfully, a browser window will popup. Don't forget to go back to your terminal and press `CTRL+C` to kill the server. To leave the virtual environment, you can run the command `deactivate`.

Now every time you want to use anvi'o, you will need to activate the virtual environment. If you like things to be convenient as much as we do, you may want to run the following command so you have a new command, `anvi-activate` that activates your anvi'o installation:

``` bash
echo 'alias anvi-activate="source ~/virtual-envs/anvio/bin/activate"' >> ~/.bash_profile
```

When I open a new terminal, things look like this:

``` bash
meren ~ $ anvi-interactive -v
-bash: anvi-interactive: command not found
meren ~ $ anvi-activate
(anvio) meren ~ $ anvi-interactive -v
Anvio version ................................: 2.1.0
Profile DB version ...........................: 17
Contigs DB version ...........................: 8
Pan DB version ...............................: 4
Samples information DB version ...............: 2
Genome data storage version ..................: 1
Auxiliary data storage version ...............: 3
Anvi server users data storage version .......: 1
(anvio) meren ~ $ 
```

### Installing or updating from the active codebase (because why not)

This will allow you to go beyond the stable version and follow the very current version of the codebase (we assume you already have taken of your dependencies).

Let's setup a new virtual environment and activate it:

``` bash
virtualenv ~/virtual-envs/anvio-dev
source ~/virtual-envs/anvio-dev/bin/activate
```

### I need to get the codebase

So this is your first time with the codebase. Get a fresh copy (with all the submodules necessary):

``` bash
cd
git clone --recursive https://github.com/meren/anvio.git
```

Then go into the `anvio` directory, and then run the installation:

``` bash
cd anvio
source ~/virtual-envs/anvio-dev/bin/activate
pip install -r requirements.txt
python setup.py install
```

### I already have the codebase

So you want to _update_ your already existing installation. Follow these steps:

{% highlight bash %}
cd
cd anvio
git pull
git submodule update --init --recursive
source ~/virtual-envs/anvio-dev/bin/activate
pip install -r requirements.txt
python setup.py install
{% endhighlight %}

### What now?

Now it is time to run `anvi-self-test --suite mini`, of course.

If you want to make things simpler, you can add an alias to your `~/.bash_profile` to easily switch to this environment:

``` bash
echo 'alias anvi-activate-dev="source ~/virtual-envs/anvio-dev/bin/activate"' >> ~/.bash_profile
```

### Installation for developers (you're a wizard, arry)

If you are planning to do this, you really need no introductions, but I will give you one anyway. Clone the codebase into a `$DIR` you like:

{% highlight bash %}
cd $DIR
git clone --recursive https://github.com/meren/anvio.git
{% endhighlight %}

Create a virtual environment (`master` to remind you that you are following the GitHub `master`), and do the initial setup, and leave it:

``` bash
virtualenv ~/virtual-envs/anvio-master
source ~/virtual-envs/anvio-master/bin/activate
cd $DIR/anvio # don't forget to update the $DIR with the real path
pip install -r requirements.txt
python setup.py build
cp build/lib.*/anvio/*so anvio/
rm -rf anvio.egg-info build dist
deactivate
```

Then update your activation batch to add necessary environment variables (keep in mind that you need to update the `$DIR` variable with whatever it shows in your system):


``` bash
echo 'export PYTHONPATH=$PYTHONPATH:$DIR/anvio/' >> ~/virtual-envs/anvio-master/bin/activate
echo 'export PATH=$PATH:$DIR/anvio/bin:$DIR/anvio/sandbox' >> ~/virtual-envs/anvio-master/bin/activate
```

That's it. If you like, add an alias to your `~/.bash_profile` to activate this quickly:

``` bash
echo 'alias anvi-activate-master="source ~/virtual-envs/anvio-master/bin/activate"' >> ~/.bash_profile
source ~/.bash_profile
```

Finally, if you would like to pull the latest commits from GitHub every time you switch to the `master`, add these to your activation batch (you will need to update `$DIR` once again):

``` bash
echo 'cd $DIR && git pull && cd -' >> ~/virtual-envs/anvio-master/bin/activate
```

You are golden.


## Running the "Mini Test"

You can make anvi'o test itself by running the program `anvi-self-test`.

Upon the successful completion of all the tests, your browser should popup to take you to the interactive interface. When you click that 'Draw' button, you should see something like this (this is one of the older version of [the anvi'o interactive interface]({% post_url anvio/2016-02-27-the-anvio-interactive-interface %}), and it shall stay here so we remember where we came from):

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/misc/mini-test-screenshot.png"><img src="{{ site.url }}/images/anvio/misc/mini-test-screenshot.png" width="50%" /></a>
</div>

All fine? Perfect! Now you have a running installation of anvi'o!

It is time to go through the [tutorial]({% post_url anvio/2016-06-22-anvio-tutorial-v2 %}), or take a look at [all the other posts on the platform]({{ site.url }}/software/anvio).

---
layout: post
title: "OS X installer for anvi'o"
excerpt: "Do you want to try anvi'o? Do you have a MAC computer? We got you covered."
modified: 2015-08-20
tags: []
categories: [anvio]
comments: true
authors: [meren]
---

{:.warning}
This is petty old now, and we will not add new OS X installers. Please do not use this, and try Homebrew installation instructions you can find in [this post]({% post_url anvio/2016-06-26-installation-v2 %}).

{% include _project-anvio-version.html %}


One of our highest priorities with [anvi'o]({{ site.url }}/software/anvio) has always been the ease of use. Pain-free installation is clearly an important requirement to achieve that. [Installing anvi'o]({% post_url anvio/2016-06-26-installation-v2 %}) is not terribly hard, however, there is still room for improvement.

We now have a new installer, and if you are using Mac OS X (version 10.9+, i.e., Maverick or Yosemite), you are pretty much covered!

The **only** thing you need to do is to download the installer on your MAC computer, and install anvi'o with one double-click and one drag gesture with your mouse.


## Download the installer (v2.0.2)

You can download the latest version of the installer using the download link (it is on the bottom right corner):

<iframe src="https://widgets.figshare.com/articles/3505763/embed?show_title=1" width="100%" height="140" frameborder="0"></iframe>


## Install anvi'o

Double-click on the installer. You should be welcomed by this screen:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2015-08-20-installation-on-mac/welcome-screen.png"><img src="{{ site.url }}/images/anvio/2015-08-20-installation-on-mac/welcome-screen.png" /></a>
</div>

Drag the anvi'o logo into the Applications folder.

## Run the test.

{: .notice}
Before I continue, I would first like to invite you to install [iTerm 2](https://www.iterm2.com/) on your Mac if you are not already using it. It is a Terminal replacement, and it is much more talented than the default terminal application on OS X systems. If iTerm 2 is installed on your computer, anvi'o will happily start in an iTerm 2 window (otherwise it will fall back to the ugly Terminal). Please see FAQ down below if you have a different terminal application.

Once the installation is completed, click the anvi'o logo in your Applications folder. This should open a terminal window.

{:.notice}
If your MAC does not open the application due to your security settings, you may need to do the following just for once: open your Applications folder, find anvi'o icon (right-click on it), and select `Open` from the right-click menu. Once you have done it, you can start anvi'o just like any other application.

To run the default test, simply type this command:

{% highlight bash %}
anvi-self-test
{% endhighlight %}

This should do a bunch of things, and if everything is in order, your browser should pop-up, and you should be able to see this display when you click the 'Draw' button:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2015-08-20-installation-on-mac/test-screen.png"><img src="{{ site.url }}/images/anvio/2015-08-20-installation-on-mac/test-screen.png" style="width: 80%;" /></a>
</div>

Do you see it? Now you can follow the [user tutorial]({% post_url anvio/2016-06-22-anvio-tutorial-v2 %}) if you like.

Thank you for giving it a try.

## Older versions


<iframe src="https://widgets.figshare.com/articles/3471440/embed?show_title=1" width="100%" height="140" frameborder="0"></iframe>

---

<iframe src="http://wl.figshare.com/articles/1613293/embed?show_title=1" width="100%" height="140" frameborder="0"></iframe>

---

<iframe src="http://wl.figshare.com/articles/1574141/embed?show_title=1" width="100%" height="140" frameborder="0"></iframe>

---

<iframe src="http://wl.figshare.com/articles/1535533/embed?show_title=1" width="100%" height="140" frameborder="0"></iframe>

---

<iframe src="http://wl.figshare.com/articles/1513854/embed?show_title=1" width="100%" height="140" frameborder="0"></iframe>

## FAQ

### Is there a way to start anvi'o in any terminal window?

Yes! Just copy-paste this line into any terminal window you are working in:

{% highlight bash %}
source /Applications/Anvio.app/Contents/Resources/Scripts/activate.sh
{% endhighlight %}

Even better, you can add an `anvio` alias in your shell, so you can invoke the platform anytime from your terminal by typing `anvio`. Run this command line on your terminal, and start a new terminal to test whether `anvio` alias is working for you:

{% highlight bash %}
echo 'alias anvio="source /Applications/Anvio.app/Contents/Resources/Scripts/activate.sh"' >> ~/.bash_profile
{% endhighlight %}

### Now it is installed on my MAC, am I good to go with my analyses?

Yes you can do that, but I would urge you to get your system administrator to have it installed on your server as well --everything they may need to install anvi'o on a server system is explained [here]({% post_url anvio/2016-06-26-installation-v2 %}), and we would be very interested to answer any further questions. If you have anvi'o installed on your personal computer as well as on your server computer, you can do compute-intensive tasks (such as profiling, merging multiple profiles, or automatic genome binning) on the server, and download the results and visualize them, or perform human-guided genome binning on your computer using the interactive interface.

### Who is responsible for this installer?

OK. I wasn't going to mention it to not embarrass him, but since you asked, I have no other option but to tell you that [Özcan Esen](https://twitter.com/ozcanesen) spent a couple of sleepless nights to create the initial version of this installer.

---
layout: page
title: Missing binaries...
modified: 2015-10-13
excerpt: "You can't find your binaries after pip installation?"
comments: true
---

{% include _toc.html %}

So you had a successfull pip installation, yet when you type this, you get a 'command not found' error, instead of a nice looking output:

    anvi-profile -v

Unfortunately this is a common problem, because .. well, let's just say 'reasons'.

What you will do to solve it is to update your `PATH` variable, so your shell knows where did pip put all the anvi'o programs.

{: .notice}
A bit more info for the curious: PATH is an environmental variable and contains all the places to look for on your computer when you enter a command. For instance, if you type `which ls`, you will get a location, which will show where the `ls` program lives. Every time you type ls, your shell knows where to find it, because the directory that contains `ls` (and all the other programs) is in your PATH variable. You can check it yourself. Type `echo $PATH`, and see all directories that are known to your system as places where programs can live and be found.

I could give you a single line of command to solve this issue, however, the directory path on my computer may not be identical to yours (and it would require just too much heuristics to find what is the actual directory on your's). Therefore it is important you do these steps on your computer. I am sorry it is not always easy and straightforward to install open source software, and thank you very much for your patience.

## Finding where the anvi'o programs are

In order to add the proper directory into your `PATH` variable, first we need to find out where anvi'o programs are. We will ask pip.

If you type this, you will get a list of file paths for each file that is installed on your system for anvi'o:

{% highlight bash %}
$ pip show -f anvio
{% endhighlight %}

All these are relative paths. But relative to what? The root directory is at the very top of your output, which is very hard to get. But if you type this, you will get that particular line where pip installed everything anvi'o:

{% highlight bash %}
$ pip show -f anvio > /tmp/_output && head -n 20 /tmp/_output | grep 'Location'
Location: /opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages
{% endhighlight %}

From this output I learn that my anvio installation is under `/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages`. But this doesn't really help. Because I still don't know where the binaries went. Let's see what is the relative path of `anvi-profile` among the installed files:

{% highlight bash %}
$ pip show -f anvio > /tmp/_output && grep 'anvi-profile' /tmp/_output
  ../../../bin/anvi-profile
{% endhighlight %}

OK. So three directories up from the installation directory, there is a directory called `bin`, and anvi'o programs went in there. Let's see:

{% highlight bash %}
$ cd /opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages
$ cd ../../../bin
$ ls anvi-profile
anvi-profile
{% endhighlight %}

Yes. Now if you type this while in that directory, you will see that everything works magically:

{% highlight bash %}
$ ./anvi-profile -v
{% endhighlight %}

So this is the directory where anvi'o programs live on your computer:

{% highlight bash %}
$ pwd
/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin
{% endhighlight %}

Please note that the directory names on my computer and yours may differ.

## Updating the PATH variable

Now you know which directory is missing from your `PATH` variable, it is time to make sure it will always be there. The best way to do it is to add this line into your `~/.bash_profile`, a special file that is read to update all environmental variables every time you start a new terminal:

    export PATH="/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin:${PATH}"

If you don't want to try to find and open that file, you can run this command, and it will add it at the end of the file:

{% highlight bash %}
$ echo 'PATH="/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin:${PATH}"' >> ~/.bash_profile
{% endhighlight %}

## Testing if everything works

Open a new terminal and type this to make sure the directory you just added is somewhere in the output:

{% highlight bash %}
$ echo $PATH
{% endhighlight %}

If the output looks alright, then this should work this time:

{% highlight bash %}
$ anvi-profile -v
{% endhighlight %}

If this doesn't solve your problem please post a comment or send me an e-mail, and I will help you sort it out.



<div style="display: block; height: 200px;">&nbsp;</div>

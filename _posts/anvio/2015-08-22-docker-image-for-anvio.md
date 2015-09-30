---
layout: post
title: "Docker image for anvi'o"
excerpt: "Fresh anvi'o builds for the lazy."
modified: 2015-08-22
tags: []
categories: [anvio]
comments: true
author: meren
---

{% include _toc.html %}

If you are reading this, you probably already know what [Docker](https://www.docker.com/https://www.docker.com/) is. If you don't know, here is a quote from their web site:

<blockquote>
An open platform for distributed applications for developers and sysadmins.
</blockquote>

I have known about Docker for a long time, however I never really had a chance to use it. As so many of you already know, being a researcher and not keeping up with the latest developments outside of your immediate field is _disturbingly_ easy.

After reading our [oldschool installation manual]({% post_url anvio/2015-05-01-installation %}), my good friend [Çağlar Onur](https://twitter.com/caglar10ur) sent me [a pull request](https://github.com/meren/anvio/pull/191) with an initial recepie to build Docker images for anvi'o. Following his contribution, we will maintain [anvi'o docker images online](https://hub.docker.com/r/meren/anvio/tags/).

I have to admit that my introduction to Docker was similar to meeting someone you know you will not get along well simply because (1) everyone talks about them non-stop, and (2) you are a person who learned to not trust what is popular. But then, after playing with it for just two days I can see that Docker does deserve to be talked about, plus, beyond the hype that surrounds the technology, it seems it is the perfect solution to get anvi'o running on a server system with 0 pain.

## Running the anvi'o docker image

I assume the server you are working on has Docker installed, and the docker service is up and running. I run these commands on my Ubuntu Linux box to install Docker, and running it: 

{% highlight bash %}
$ sudo apt-get install docker.io
$ sudo service docker.io start
{% endhighlight %}

Now you can simply type this to get all the anvi'o tags available to you:

{% highlight bash %}
$ sudo docker pull meren/anvio
{% endhighlight %}

Here are the available tags as of today:

{% highlight bash %}
$ sudo docker images               
REPOSITORY          TAG                 IMAGE ID            CREATED             VIRTUAL SIZE
meren/anvio         1.1.0               71a85f4a8a68        3 weeks ago         482.6 MB
meren/anvio         1.0.2               f3c1998acc63        5 weeks ago         481.3 MB
{% endhighlight %}

{:.notice}
Although we intent to do our best to keep this post up-to-date, please remember that there may be a newer docker image at any given time. In that case, please use the most recent docker image and make sure the user tutorial matches to the version you are using.

Now you can simply type this to get the version you are interested in running:

{% highlight bash %}
$ sudo docker run --rm -it meren/anvio:1.1.0
{% endhighlight %}

This will give you a new prompt in your terminal that will look like this:

{% highlight bash %}
:: anvi`o ::  / >>>
{% endhighlight %}

If you would like to run the mini test, you can type these commands:

{% highlight bash %}
:: anvi`o ::  / >>> cd anvio-tests/
:: anvi`o ::  /anvio-tests >>> ./run_mini_test.sh 
{% endhighlight %}

You can press CTRL+D to go back to your host terminal.

### How about my files?

You will soon realize that this new virtual environment in which you run anvi'o does not have access to the rest of your disk. Let's assume the data you would like to analyze using anvi'o is in a directory at `/home/meren/my_data`. To have access to this directory from within the docker image, you can _bind_ it to a directory in the virtual environment. Let's assume, we want to have access to the content of `/home/meren/my_data`, from `/my_data` directory from within the docker image:

{% highlight bash %}
$ sudo docker run --rm -v /home/merenbey/my_data:/my_data -it meren/anvio:1.1.0
{% endhighlight %}

That's it. Here is a longer story:

{% highlight bash %}
$ pwd
/home/meren
$ ls my_data/
X.bam  Y.bam  contigs.fa
$ sudo docker run --rm -v /home/merenbey/my_data:/my_data -it meren/anvio:1.1.0                                      
:: anvi`o ::  / >>> ls /my_data/
X.bam  Y.bam  contigs.fa
{% endhighlight %}

## Rebuilding the anvi'o docker image

If you would like to rebuild the Docker image for anvi'o on your own server, you can simply clone the anvi'o repository:

{% highlight bash %}
$ git clone --recursive https://github.com/meren/anvio.git
{% endhighlight %}

Switch to a particular tag (here is [a complete list of available releases and tags](https://github.com/meren/anvio/releases)):

{% highlight bash %}
git checkout tags/v1.1.0
{% endhighlight %}

Add/remove things you want, do your changes in the code, and build the new docker image (replace the username and tag with your preferences):

{% highlight bash %}
docker build -t meren/anvio:1.1.0 .
{% endhighlight %}

And optionally, push it to your account on the hub to allow other people run it easily:

{% highlight bash %}
docker push meren/anvio
{% endhighlight %}

Don't hesitate to get in touch if you have any questions, and check out the [anvi'o project page]({{ site.url }}/projects/anvio).

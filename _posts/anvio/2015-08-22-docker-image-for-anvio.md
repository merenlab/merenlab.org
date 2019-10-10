---
layout: post
title: "Docker image for anvi'o"
excerpt: "Fresh anvi'o builds for the lazy."
modified: 2015-08-22
tags: []
categories: [anvio]
comments: true
authors: [meren]
---

{% include _project-anvio-version.html %}

{% include _toc.html %}

{% include _join-anvio-slack.html %}

If you want to run anvi'o quickly on a server or a desktop system with zero pain and no installation, this could be a good solution for you (otherwise try the [standard way of installing anvi'o]({% post_url anvio/2016-06-26-installation-v2 %}). 

{:.notice}
Special thanks go to to my good friend [Çağlar Onur](https://twitter.com/caglar10ur) for forcing us to [have it](https://github.com/meren/anvio/pull/191).

{:.warning}
**Our latest docker images are quite big** (~8 Gb), because images after `v6` come with lots of bells and whistles such as assemblers, mapping software, binning tools (that seems to be quite hard to install in general), on other useful tools. Installing them is pain, so we go through it once to create you working recipes, and storage is cheap. But if you need a minimal anvi'o image for your specific needs, please let us know, and we can help you with that.

## Running an anvi'o container

Here we assume that you have Docker installed, and the docker service is up and running. Just to make sure, type this in your terminal, and expect to see an output similar to mine below:

```
docker --version

Docker version 19.03.2, build 6a30dfc
```

If you get a command not found error, you may need to install docker, or start its service.

The following command display all the anvi'o docker images available to you:

``` bash
$ docker pull meren/anvio
```

If you are seeing something like this when you are done pulling everything, you are golden:

``` bash
$ docker images

REPOSITORY          TAG        IMAGE ID     CREATED    SIZE
meren/anvio         6          (...)        (...)      7.9GB
meren/anvio         5.5        (...)        (...)      1.37GB
meren/anvio         5          (...)        (...)      1.45GB
meren/anvio         4          (...)        (...)      900.0 MB
meren/anvio         3          (...)        (...)      932.0 MB
meren/anvio         2.4.0      (...)        (...)      822.0 MB
meren/anvio         2.3.2      (...)        (...)      806.0 MB
meren/anvio         2.3.0      (...)        (...)      804.7 MB
meren/anvio         2.1.0      (...)        (...)      667.4 MB
meren/anvio         2.0.2      (...)        (...)      672.0 MB
meren/anvio         2.0.1      (...)        (...)      671.9 MB
meren/anvio         1.2.3      (...)        (...)      595.2 MB
```

The best way to start an anvi'o instance is to first changing your working directory in your terminal to where the files you wish to work with are, and then run this command:

``` bash
docker run --rm -it -v `pwd`:`pwd` -w `pwd` -p 8080:8080 meren/anvio:latest
```

Following this command, you should see a new command line like this one:

``` bash
:: anvi'o v6 :: /WORK/DIR/PATH >>>
```

You're done! If you type `anvi-` and press `TAB` twice, you see all the anvi'o programs available to you.

{:.warning}
Please note that if you run an anvi'o interactive interface, it will not automatically fire up in your browser. You will have to manually visit the address http://localhost:8080 to interact with your data. Once you are done with your interactive interfaces you can press `CTRL`+`C` to kill the server.

{:.notice}
If you would like to leave docker environment and go back to your terminal, you can press `CTRL`+`D`.

## Installing new things into your image

The image comes with a minimal set of programs to make sure it is not too much of a burden to download, but since it is a running Ubuntu system, you can indeed install and run anything. Example:

``` bash
apt-get update
apt-get install vim
```

Done!

## Rebuilding the anvi'o docker image (for hackers)

We are keeping our docker recipes in the main anvi'o repository:

[https://github.com/merenlab/anvio/blob/master/Dockerfile](https://github.com/merenlab/anvio/blob/master/Dockerfile)

Start by getting a copy of the repository:

```
git clone https://github.com/merenlab/anvio.git
cd anvio/
```

In theory, at any given time you should be able to rebuild a docker image for anvi'o `master` this way:

```
docker build -t meren/anvio:master .
```

Alternatively, you can switch to a tag, such as anvi'o `v6`, and build a copy of that:

If you would like to rebuild the Docker image for anvi'o on your own server, you can simply get a copy of the file:

``` bash
git checkout v6
docker build -t meren/anvio:6 .
```

Optionally, you can push it to your account on the hub to allow other people run it easily (i.e., this is what we do to push new images to the anvi'o accound):

``` bash
docker push meren/anvio
```

Don't hesitate to get in touch if you have any questions!

{% include _join-anvio-slack.html %}


## Rebuilding the anvi'o docker image (for versions older than v6)

Before we started generating giant docker images with everything in them, we used to keep our docker files on gist and build our images only from anvi'o tarballs we submitted to pip. You can find a docker file for each version of anvi'o here:

[https://gist.github.com/meren/65b1f1bfea1b53e87e10f025d1e4c29a](https://gist.github.com/meren/65b1f1bfea1b53e87e10f025d1e4c29a)

For instance, if you would like to build a docker image for `v5.5`, you can get a copy of the docker recipe,

``` bash
wget https://gist.githubusercontent.com/meren/65b1f1bfea1b53e87e10f025d1e4c29a/raw/b0a841c4be8873ba73faba217e798d3d9e207823/Dockerfile_v5.5.sh -O Dockerfile
```

Add/remove things you want to the image, and build it:

``` bash
docker build -t meren/anvio:5.5 .
```


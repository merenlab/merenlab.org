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

If you are reading this, you probably already know what [Docker](https://www.docker.com/) is. If you don't know, here is a quote from their web site:

<blockquote>
An open platform for distributed applications for developers and sysadmins.
</blockquote>

If you do not want to deal with the [oldschool installation manual]({% post_url anvio/2016-06-26-installation-v2 %}), you can try anvi'o using its Docker image (special thanks go to to my good friend [Çağlar Onur](https://twitter.com/caglar10ur) for forcing us to [have it](https://github.com/meren/anvio/pull/191)). We now maintain [anvi'o docker images online](https://hub.docker.com/r/meren/anvio/tags/) as closely as possible to the latest stable release.

If you want to quickly run anvi'o on a server or a desktop system with zero pain, this could be a good solution for you.

## Running anvi'o docker image

Here we assume that you have Docker installed, and the docker service is up and running. Then you can simply type this command to get all the anvi'o docker images available to you:

``` bash
$ docker pull meren/anvio
```

If you are seeing something like this when you are done pulling everything, you are golden:

``` bash
$ docker images
REPOSITORY          TAG        IMAGE ID     CREATED    SIZE
meren/anvio         5.4        (...)        (...)      1.37GB
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

You are now ready to run anvi'o. The best way to do it is to first change your working directory in your terminal to where the files you wish to work with using anvi'o are, and then run this command:

``` bash
docker run --rm -it -v `pwd`:`pwd` -w `pwd` -p 8080:8080 meren/anvio:latest
```

Following this command, you should see a new command line like this one:

``` bash
:: anvi'o v5.4 :: /WORK/DIR/PATH >>>
```

You're done! If you type `anvi-` and press `TAB` twice, you see all the anvi'o programs available to you.

{:.warning}
Please note that if you run an anvi'o interactive interface, it will not automatically fire up in your browser. You will have to manually visit the address http://localhost:8080 to interact with your data. Once you are done with your interactive interfaces you can press `CTRL`+`C` to kill the server.

{:.notice}
If you would like to leave docker environment and go back to your terminal, you can press `CTRL`+`D`.

## Installing new things into your image

The image comes with a minimal set of programs to make sure it is not too much of a burden to download, but since it is a running Ubuntu system, you can indeed install and run anything. For instance, if you are a fan of `vi` and can't do without it, this is what you could do:

``` bash
apt-get update
apt-get install vim
```

Done!

## Rebuilding the anvi'o docker image (for hackers)

For your reference, We keep our Docker files for each anvi'o version here:

[https://gist.github.com/meren/65b1f1bfea1b53e87e10f025d1e4c29a](https://gist.github.com/meren/65b1f1bfea1b53e87e10f025d1e4c29a)

If you would like to rebuild the Docker image for anvi'o on your own server, you can simply get a copy of the file:

``` bash
wget https://gist.githubusercontent.com/meren/65b1f1bfea1b53e87e10f025d1e4c29a/raw/49e7a0b55474c37e271d8188865d2ba423923330/Dockerfile_v5.4.sh -O Dockerfile
```

Add/remove things you want to the image, and build the new docker image:

``` bash
docker build -t meren/anvio:MY_FANCY_ANVIO_IMAGE .
```

Optionally, you can push it to your account on the hub to allow other people run it easily (i.e., this is what we do to push new images to the anvi'o accound):

``` bash
docker push meren/anvio
```

Don't hesitate to get in touch if you have any questions!

{% include _join-anvio-slack.html %}
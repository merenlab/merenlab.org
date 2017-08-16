---
layout: post
authors: [meren]
title: "A Docker image for oligotyping"
excerpt: "For people who does not want to install the pipeline"
modified: 2017-08-15
tags: [virtual, installation]
categories: [oligotyping]
comments: true
---

{% include _toc.html %}

{:.notice}
We used to have a virtual box, but we no longer serve it.

You can start using oligotyping from within a [Docker](https://www.docker.com/https://www.docker.com/) image without installing anything on your computer.

## Running the oligotyping docker image (for Linux users)

I assume the server you are working on has Docker installed, and the docker service is up and running. I run these commands on my Ubuntu Linux box to install Docker, and running it:

``` bash
$ sudo apt-get install docker.io
$ sudo service docker.io start
```

Now you can simply type this to get all the oligotyping tags available to you:

``` bash
$ sudo docker pull meren/oligotyping
```

Here are the available tags as of today:

``` bash
$ sudo docker images
REPOSITORY          TAG        IMAGE ID     CREATED    SIZE
meren/oligotyping   latest     (...)        (...)      1.17GB
meren/oligotyping   2.1        (...)        (...)      1.17GB
```

Now you can simply type this to get the latest version running (or you can replace `latest` with the version number you see in the list of TAGs):

``` bash
$ sudo docker run --rm -it meren/oligotyping:latest
```

This will give you a new prompt in your terminal that will look like this:

``` bash
:: oligotyping ::  / >>>
```

And this is what you should see when you type 'oligotype' in the terminal:

``` bash
 :: oligotyping ::  / >>> oligotype
usage: oligotype [-h] [-o OUTPUT_DIRECTORY] [-c NUMBER_OF_AUTO_COMPONENTS]
                 [--qual-scores-file QUAL SCORES FILE]
                 [--qual-scores-dict QUAL SCORES DICT]
                 [--qual-stats-dict QUAL STATS DICT] [-q MIN_BASE_QUALITY]
                 [-C SELECTED_COMPONENTS] [-s MIN_NUMBER_OF_SAMPLES]
                 [-a MIN_PERCENT_ABUNDANCE] [-A MIN_ACTUAL_ABUNDANCE]
                 [-M MIN_SUBSTANTIVE_ABUNDANCE] [-t SAMPLE_NAME_SEPARATOR]
                 [-l LIMIT_REPRESENTATIVE_SEQUENCES]
                 [--limit-oligotypes-to LIMIT_OLIGOTYPES_TO]
                 [-e EXCLUDE_OLIGOTYPES] [--quick] [--no-figures]
                 [--blast-ref-db BLAST_REF_DB]
                 [--colors-list-file COLORS_LIST_FILE] [--do-blast-search]
                 [--no-display] [--skip-gen-html] [--generate-sets] [-K]
                 [-S COS_SIM_TR] [-E FILEPATH] [--project PROJECT]
                 [--skip-check-input-file] [--skip-basic-analyses]
                 [--skip-gexf-network-file] [-T] [-N INTEGER] [--version]
                 INPUT ALIGNMENT ENTROPY
oligotype: error: too few arguments
```

Good? Good. You can press `CTRL`+`D` to go back to your host terminal.


## Running the oligoytping docker image (for Mac OSX / Windows users)

First, you will need to install the [Docker Toolbox](https://www.docker.com/toolbox). Once you have it installed, run 'Docker Quickstart Terminal' application. At the end, this is what you want to see:

<div class="centerimg">
<a href="{{ site.url }}/images/anvio/2015-08-22-docker-image-for-anvio/docker-terminal.png"><img src="{{ site.url }}/images/anvio/2015-08-22-docker-image-for-anvio/docker-terminal.png" /></a>
</div>

Note the IP address (which is `192.168.99.100` in my case), becasue you will need it later. At this point, all you need to do is to run this command (note that you may need to change the version number):

``` bash
docker run -p 8080:8080 -it meren/oligityping:latest
```

Docker will download all the necessary images, and will finally give you a prompt which will look like this:

``` bash
:: oligotyping ::  / >>>
```

Done!

## How about my files?

You will soon realize that this new virtual environment in which you run oligotyping does not have access to the rest of your disk.

Let's assume the data you would like to analyze using oligotyping is in a directory at `/home/meren/my_data`. To have access to this directory from within the docker image, you can _bind_ it to a directory in the virtual environment. Let's assume, we want to have access to the content of `/home/meren/my_data`, from `/my_data` directory from within the docker image:

``` bash
$ sudo docker run --rm -v /home/merenbey/my_data:/my_data -it meren/oligotyping:latest
```

That's it. Here is a longer story:

``` bash
$ pwd
/home/meren
$ ls my_data/
test.fa
$ sudo docker run --rm -v /home/merenbey/my_data:/my_data -it meren/oligotyping:latest
:: anvi`o ::  / >>> ls /my_data/
test.fa
```


## Rebuilding the docker image (for hackers)

This is the current Docker image for v2.1:

[https://gist.github.com/meren/3b16c98b830e3796216b9d623ea2c2f4](https://gist.github.com/meren/3b16c98b830e3796216b9d623ea2c2f4)

If you would like to rebuild the Docker image for oligotyping on your own server, you first need to get the Docker file:

``` bash
wget https://gist.githubusercontent.com/meren/3b16c98b830e3796216b9d623ea2c2f4/raw/38e4188b0c21de89dd9eace5da361c9f9492f77d/Oligotyping_Dockerfile_v2.1.sh -O Dockerfile
```

Add/remove things you want, do your changes in the code, and build the new docker image (replace the username and tag with your preferences):

``` bash
docker build -t meren/oligotyping:2.1 .
```

And optionally, push it to your account on the hub to allow other people run it easily (i.e., this is what I do to push the new images to my account):

``` bash
docker push meren/oligotyping
```

Don't hesitate to get in touch if you have any questions.

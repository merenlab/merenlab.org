---
layout: post
authors: [meren]
title: "Oligotyping in a virtual box"
excerpt: "For people who does not want to install the pipeline"
modified: 2014-09-02
tags: [virtual, installation]
categories: [oligotyping]
comments: true
---

{% include _toc.html %}

[Çağrı Ulaş](http://about.me/cagriulas), a second year computer engineering undergrad from Canakkale Onsekiz Mart University, Turkey, created an [Arch Linux](https://www.archlinux.org/)-based [VirtualBox](https://www.virtualbox.org/) image for oligotyping. Now you can download this image, and run it on any operating system without having to install the oligotyping pipeline and its dependencies. The only thing you need to have installed is [VirtualBox](https://www.virtualbox.org/), which is free, and provides installers for any major operating system. I thank Çağrı for his useful contribution.

You can download the Arch Linux image of the oligotyping pipeline from here (1.4Gb compressed, ~6Gb uncompressed, MD5 sum 2aba169ddb811fd1e4a8ca0acd63c9b5):

[https://jbpc.mbl.edu/media/meren/oligotyping.vdi.tar.xz](https://jbpc.mbl.edu/media/meren/oligotyping.vdi.tar.xz)

## Running the virtual machine

Here is a list of instructions on how to get it running:

* Download and install VirtualBox. Download and uncompress the image file. The image is compressed with xz. Your system most probably already has the necessary tools to uncompress it (for instance I simply type “tar xvf oligotyping.vdi.tar.xz” and it works), but if not, it shouldn’t be hard to get xz working for you as it is pretty standard. What you want is to achieve the oligotyping.vdi file on your disk.

* Next, run VirtualBox, and click “New”, enter the following information, and click continue:

<figure>
	<a href="{{ site.url }}/images/oligotyping/vm-01.png"><img src="{{ site.url }}/images/oligotyping/vm-01.png"></a>
</figure>

* Virtual machines require a lot of system resources. How much you can spare is up to you. But I wouldn’t recommend allocating less than 2Gb of memory. Do your selection, and continue:

<figure>
	<a href="{{ site.url }}/images/oligotyping/vm-02.png"><img src="{{ site.url }}/images/oligotyping/vm-02.png"></a>
</figure>

* Next, it is time to select the image we just uncompressed. Once you are done, click create:

<figure>
	<a href="{{ site.url }}/images/oligotyping/vm-03.png"><img src="{{ site.url }}/images/oligotyping/vm-03.png"></a>
</figure>

* When you create the virtual machine, you will see an entry for the oligotyping on the left. Before we start, there are two more things to adjust: graphics memory and number of CPUs. If you allow your VM to use only one core, you need to remember to add `--no-threading` option at the end of every oligotyping command. If you increase the number here, you are safe. Click `settings` in the right click menu of the oligotyping entry, and navigate to System and Video tabs to set following configurations, and click OK to go back to the main screen:


<figure>
	<a href="{{ site.url }}/images/oligotyping/vm-04.png"><img src="{{ site.url }}/images/oligotyping/vm-04.png"></a>
</figure>

<figure>
	<a href="{{ site.url }}/images/oligotyping/vm-05.png"><img src="{{ site.url }}/images/oligotyping/vm-05.png"></a>
</figure>

* Now you can start the VM, which should welcome you with a clean looking desktop. If you click the terminal and type `oligotype`, you should get the help menu, in which case you can get your files into your virtual machine and start working on them:


<figure>
	<a href="{{ site.url }}/images/oligotyping/vm-06.png"><img src="{{ site.url }}/images/oligotyping/vm-06.png"></a>
</figure>

The final thing you must do is to upgrade the oligotyping pipeline on the virtual box. You can do it by typing this command in the terminal:

    sudo pip install --upgrade oligotyping

{: .notice}
This command will ask you to enter the root password. The default root password (as well as the user password) is "**oligo**".

 

I hope this would be helpful to people who would like to take quick look at the oligotyping pipeline without having to invest too much time for installation.

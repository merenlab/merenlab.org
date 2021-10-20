---
layout: page
title: configuration-ini [artifact]
categories: [anvio]
comments: false
redirect_from: /7.1/configuration-ini
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/TXT.png" alt="TXT" style="width:100px; border:none" />

A TXT-type anvi'o artifact. This artifact is typically provided **by the user** for anvi'o to import into its databases, process, and/or use.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


There are no anvi'o tools that generate this artifact, which means it is most likely provided to the anvi'o ecosystem by the user.


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-script-gen-short-reads](../../programs/anvi-script-gen-short-reads)</span></p>


## Description

This describes an INI formatted file used to configure a program. 

If you're unsure what the INI format is, you can check out its [Wikipedia page](https://en.wikipedia.org/wiki/INI_file). But it is essentially a text file with sections that defines the values of several keys.

As of now, it is only required in Anvi'o by the program <span class="artifact-n">[anvi-script-gen-short-reads](/software/anvio/help/7.1/programs/anvi-script-gen-short-reads)</span>, where the file describes various parameters for the short reads that you want to generate, such as the desired length and coverage. Take a look at the page for that program for an example. 


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/configuration-ini.md) to update this information.


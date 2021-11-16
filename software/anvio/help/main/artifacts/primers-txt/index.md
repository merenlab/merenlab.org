---
layout: page
title: primers-txt [artifact]
categories: [anvio]
comments: false
redirect_from: /m/primers-txt
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


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-script-get-primer-matches](../../programs/anvi-script-get-primer-matches)</span></p>


## Description

A **TAB-delimited** file to describe primer sequences. A primer sequence can be exact (such as `ATCG`), or fuzzy (such as `AT.G`, which would match any combination of `ATAG`, `ATTG`, `ATCG`, or `ATGG`). Fuzzy primers are defined by regular expressions, properties of which are explained best in [the Python documentation](https://docs.python.org/3/library/re.html), or cheatsheets like [this one](https://www.debuggex.com/cheatsheet/regex/python).

This file type includes two required and any number of user-defined optional columns.

The following three columns are **required** for this file type:

* `name`: a single-word name for a given primer,
* `primer_sequence`: the primer sequence.

### An example primers-txt file

Here is an example file with three primers:

|name|primer_sequence|
|:--|:--|
|PR01|AA.A..G..G..G.CCG.C.A.C|
|PR02|AACACCGCAGTCCATGAGA|
|PR03|A[TC]A[CG]T[ATC]TCGAGC|


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/primers-txt.md) to update this information.


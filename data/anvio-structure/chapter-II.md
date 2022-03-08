---
layout: page
title: Chapter II - Reproducing Kiefl et al, 2022
modified: 2021-10-21
excerpt: "A complete reproducible workflow of the manuscript 'Structure-informed microbial population genetics elucidate selective pressures that shape protein evolution' by Kiefl et al"
comments: true
authors: [evan]
redirect_from:
  - 
---

{% capture images %}{{site.url}}/data/anvio-structure/images{% endcapture %}
{% capture command_style %}background: #D7484822; border: 4px solid #D74848;{% endcapture %}
{% capture analysis_style %}background: #E6DBE4{% endcapture %}

{:.warning}
This document is **UNDER CONSTRUCTION**. It is not in a state where you can yet reproduce our work. We anticipate this workflow will be finalized by late March, and will remove this message when it is complete.

## Quick Navigation

- [Chapter I: The prologue]({{ site.url }}/data/anvio-structure/chapter-I)
- [Chapter II: Configure your system]({{ site.url }}/data/anvio-structure/chapter-II) ← _you are here_
- [Chapter III: Build the data]({{ site.url }}/data/anvio-structure/chapter-III)
- [Chapter IV: Analyze the data]({{ site.url }}/data/anvio-structure/chapter-IV)
- [Chapter V: Reproduce every number]({{ site.url }}/data/anvio-structure/chapter-V)


## Step X: Setting up the required environment

<div class="extra-info" markdown="1">
<span class="extra-info-header">Step X Info</span>
‣ **Prerequisite steps:** None FIXME add real example  
‣ **Checkpoint datapack:** None  
‣ **Internet:** Yes
</div>

This is the start ofthe workflow. It is the least fun and most important part of this workflow: installing all required programs and resolving all dependency issues and versioning.

Constructing an identical computing environment to mine is essential for replicating this analysis. Though reproducibility depends critically upon this step, it is contradictorily the part most commonly absent in "reproducible" workflows.

The computing environment for this study, which houses all of the programs and their correct versions, is sort of like the axioms of this study. Start with different axioms and you're going to get a different result, and the most probable result will be errors. Errors everywhere.

For example, if you go rogue and say, I will use version `1.0.1` of the Python package `pandas` instead of the required version `0.25.1`, you won't replicate this study--you'll instead get errors such as

```
TypeError: concat() got an unexpected keyword argument 'join_axes'" as join_axes is deprecated starting from Pandas v1.0.
```

or variants thereof. It is really that simple. Everything must be identical, otherwise you run the extreme risk of getting stuck.

Fortunately, creating an environment identical to mine is _way_ easier than it sounds, thanks to environment control systems like Docker and conda. So don't be scared.

### Overview of required programs and libraries

- list of required programs and a little description of why each is required (preface list with note that instructions are below)

(
Notes to self:
- pymol
- conda install -y -c salilab dssp
- conda install -c bioconda paml
- pip install xlrd==1.2.0
- pip install openpyxl==3.0.9
- pip install prody==2.0.1
- pip install tmscoring   # for TM score
- pip install 'iminuit<2' # for TM score
- run all required anvi-setup programs
- install instructions (include all R packages)
- install.packages("spatstat") # used for weighted KS test
- must have perl for pal2nal.pl (dnds stuff)
- run tests to make sure you got everything installed right
)

### Option #1: Docker

- much preffered
- if this solution works for you, use it

### Option #2: conda

- sometimes docker ain't an option
- that means you're forced to do the next best thing, conda
- this requires more setup

"_The privilege needed to start docker containers is basically super-user level, making it inappropriate for shared computing environments with “ordinary” users who need access to the software there.

Conda environments are basically able to do what the user can do, and importantly neither requires nor grants administrative privilege outside of the environment._"

- https://www.quora.com/What-is-the-difference-between-conda-and-Docker

### Option #3: rogue

- Good luck

### Perform all tests

It is much, much better to figure out that something is misconfigured now, rather than later. Please take the time ensure you can run all of these commands without error.

{:.notice}
Are you unable to complete the tests and need help? Please reach out to me (I'm Evan, you can find my contact info [here](https://merenlab.org/people/)), so I can (1) help you and (2) make these instructions clearer for all others who follow in your trailblazing steps.

{:.notice}
Did you complete all the tests but you've ran into downstream errors? Please reach out to me (I'm Evan, you can find my contact info [here](https://merenlab.org/people/)) so I can (1) help you and (2) make these instructions clearer for all others who follow in your trailblazing steps.


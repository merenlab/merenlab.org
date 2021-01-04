---
layout: page
title: anvi-script-transpose-matrix [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Transpose a TAB-delimited file.

See **[program help menu](../../../../vignette#anvi-script-transpose-matrix)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[view-data](../../artifacts/view-data) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[functions-txt](../../artifacts/functions-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[misc-data-items-txt](../../artifacts/misc-data-items-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[misc-data-layers-txt](../../artifacts/misc-data-layers-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[gene-calls-txt](../../artifacts/gene-calls-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[linkmers-txt](../../artifacts/linkmers-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[view-data](../../artifacts/view-data) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[functions-txt](../../artifacts/functions-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[misc-data-items-txt](../../artifacts/misc-data-items-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[misc-data-layers-txt](../../artifacts/misc-data-layers-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[gene-calls-txt](../../artifacts/gene-calls-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[linkmers-txt](../../artifacts/linkmers-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>

## Usage


This is a script that transposes tab-delimited files. That's it. 

It's helpful to get your inputs to line up with the types of inputs that anvi'o expects. Some programs have the `--transpose` flag, which will run this program for you, but some don't, and that's when you'll have to run it yourself. 

For example, anvi'o expects <span class="artifact-n">[view-data](/software/anvio/help/7/artifacts/view-data)</span> to have each column representing a sample. If the file that you want to integrate into your anvi'o project has the samples as rows and the data attribute as the columns, then you'll need to <span class="artifact-n">[anvi-script-transpose-matrix](/software/anvio/help/7/programs/anvi-script-transpose-matrix)</span> it. 

### An Example Run 

If you have an input ile `INPUT.txt` that looks like this: 

    1   2   3   
    4   5   6   
    7   8   9
    10  11  12
    
And you run this:

<div class="codeblock" markdown="1">
anvi&#45;script&#45;transpose&#45;matrix &#45;o INPUT_transposed.txt \
                             &#45;i INPUT.txt 
</div>

You'll get a file called `INPUT_transposed.txt` that looks like 

    1   4   7   10
    2   5   8   11
    3   6   9   12
    



{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-script-transpose-matrix.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-transpose-matrix) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

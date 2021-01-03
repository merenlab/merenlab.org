---
layout: page
title: anvi-script-snvs-to-interactive [program]
categories: [anvio]
comments: false
image:
  featurerelative: ../../../images/header.png
  display: true
---

Take the output of anvi-gen-variability-profile, prepare an output for interactive interface.

See **[program help menu](../../../../vignette#anvi-script-snvs-to-interactive)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[interactive](../../artifacts/interactive) <img src="../../images/icons/DISPLAY.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[variability-profile-txt](../../artifacts/variability-profile-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>

## Usage


This programs takes a <span class="artifact-n">[variability-profile-txt](/software/anvio/help/7/artifacts/variability-profile-txt)</span> and generates the information necessary to visualize its contents with <span class="artifact-n">[anvi-interactive](/software/anvio/help/7/programs/anvi-interactive)</span>. 

Specifically, this program outputs a directory that contains a <span class="artifact-n">[profile-db](/software/anvio/help/7/artifacts/profile-db)</span>, a <span class="artifact-n">[view-data](/software/anvio/help/7/artifacts/view-data)</span> artifact, and a <span class="artifact-n">[dendrogram](/software/anvio/help/7/artifacts/dendrogram)</span>. For example, if you ran this program like so: 

<div class="codeblock" markdown="1">
anvi&#45;script&#45;snvs&#45;to&#45;interactive &#45;o OUTPUT_DIR \
                                <span class="artifact&#45;n">[variability&#45;profile](/software/anvio/help/7/artifacts/variability&#45;profile)</span>
</div>

Then, you can open the interactive interface by running 

<div class="codeblock" markdown="1">
anvi&#45;interactive &#45;&#45;manual&#45;mode \
                 &#45;p OUTPUT_DIR/profile.db \
                 &#45;&#45;tree OUTPUT_DIR/tree.txt \
                 &#45;&#45;view&#45;data OUTPUT_DIR/view.txt
</div>

## Other parameters 

### Using Only a Subset of the Input

By default, all variability positions in your variability profile are considered. However, if the input is too large (i.e. more than 25,000 variability positions), the runtime on this program will be very long and the results won't display well. So, there are several ways to remove variability positions from  the input to get under this threshold: 

1. Ignore positions with with certain departures from the consensus sequence (with `--min-departure-from-consensus` and `--max-departure-from-consensus`)
2. Ignore positions with with certain departures from the reference sequence (with `--min-departure-from-reference` and `--max-departure-from-reference`)
3. Ignore positions in all non-coding regions with the flag `--only-in-genes`. 

If you still have more positions than you can tell the program to pick a random subset of the input with the parameter `--random` followed by a seed integer.

### Modifying the Output

By the default, the output data will use the departure from consensus values. If instead you want to look at the departure from the reference, just add the falg `--display-dep-from-reference`


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-script-snvs-to-interactive.md) to update this information.


## Additional Resources


* [Use in the Infant Gut Tutorial](http://merenlab.org/tutorials/infant-gut/#visualizing-snv-profiles-using-anvio)


{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-snvs-to-interactive) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

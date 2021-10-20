---
layout: page
title: anvi-report-linkmers [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-report-linkmers
image:
  featurerelative: ../../../images/header.png
  display: true
---

Reports sequences stored in one or more BAM files that cover one of more specific nucleotide positions in a reference.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[bam-file](../../artifacts/bam-file) <img src="../../images/icons/BAM.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[linkmers-txt](../../artifacts/linkmers-txt) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Usage


Reports sequences stored in a <span class="artifact-n">[bam-file](/software/anvio/help/7.1/artifacts/bam-file)</span> file that cover one of more specific nucleotide positions in a reference.

### Basic mode of operation

Assume you wish to recover reads stored in one or more BAM files, where the matching reads contain at least one nucleotide position that align to a nucleotide position `nucleotide_position_N` in a contig `contig_name_X`. In that case, the user would first generate a two column TAB-delmited file, for example `positions_for_linkmers.txt` with no header line,


<table>
  <tbody>
    <tr>
      <td> contig_name_X </td>
      <td> nucleotide_position_N </td>
    </tr>
  </tbody>
</table>

And run the program this way to recover the short reads from this way:

```
anvi-report-linkmers --contigs-and-positions positions_for_linkmers.txt \
                     -i SAMPLE_01.bam SAMPLE_02.bam SAMPLE_03.bam (...) \
                     -o linkmers.txt
```

The user can define multiple contigs in the input file, and one or more nucleotide positions for each one of them:

<table>
<tbody>
<tr>
      <td> contig_name_X </td>
      <td> nucleotide_position_01,nucleotide_position_02,nucleotide_position_03</td>
</tr>
<tr>
      <td> contig_name_Y </td>
      <td> nucleotide_position_04 </td>
</tr>
<tr>
      <td> contig_name_Z </td>
      <td> nucleotide_position_05,nucleotide_position_06 </td>
</tr>
</tbody>
</table>

The resulting <span class="artifact-n">[linkmers-txt](/software/anvio/help/7.1/artifacts/linkmers-txt)</span> would include all short reads that match any of these critera

### Complete or incomplete links?

Using the `--only-complete-links` flag, the user can enforce whether only complete links should be reported where each reported short read must cover each nucleotide position for a given contig.

Please note that if the nucleotide positions chosen for a given contig are too distant from each other given the short read length, zero reads may satisfy the complete links criterion.

Having complete links, however, will enable [oligotyping](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12114) analyses on **metagenomic reads** through the anvi'o program <span class="artifact-n">[anvi-oligotype-linkmers](/software/anvio/help/7.1/programs/anvi-oligotype-linkmers)</span>.

### See this program in action

[http://merenlab.org/2015/12/09/musings-over-commamox/](http://merenlab.org/2015/12/09/musings-over-commamox/)

{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-report-linkmers.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-report-linkmers) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

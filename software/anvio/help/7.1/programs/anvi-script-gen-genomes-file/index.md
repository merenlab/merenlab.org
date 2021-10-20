---
layout: page
title: anvi-script-gen-genomes-file [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-script-gen-genomes-file
image:
  featurerelative: ../../../images/header.png
  display: true
---

Generate an external genomes or internal genomes file.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/ivagljiva.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Iva Veseli</span><div class="page-author-social-box"><a href="mailto:iveseli@uchicago.edu" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://github.com/ivagljiva" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[collection](../../artifacts/collection) <img src="../../images/icons/COLLECTION.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[external-genomes](../../artifacts/external-genomes) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[internal-genomes](../../artifacts/internal-genomes) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Usage


The primary purpose of this script is to reduce the amount of labor required to generate <span class="artifact-n">[external-genomes](/software/anvio/help/7.1/artifacts/external-genomes)</span> or <span class="artifact-n">[internal-genomes](/software/anvio/help/7.1/artifacts/internal-genomes)</span> files anvi'o typically uses to learn about your bins and/or genomes.

## Generating an external genomes file

If you provide an input directory and a name for the output file, then every <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> in that directory will get a line in the resulting <span class="artifact-n">[external-genomes](/software/anvio/help/7.1/artifacts/external-genomes)</span> file:

```
anvi-script-gen-genomes-file --input-dir path/to/dir \
                             --output-file external_genomes.txt
```

Names for genomes in the the resulting external genomes file will be set based on the `project_name` variable, and the `contigs_db_path` column will contain absolute paths.

{:.notice}
You can learn the current `project_name` and/or change it for a given <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span> using the program <span class="artifact-n">[anvi-db-info](/software/anvio/help/7.1/programs/anvi-db-info)</span>. This variable is set by the program <span class="artifact-n">[anvi-gen-contigs-database](/software/anvio/help/7.1/programs/anvi-gen-contigs-database)</span>.

You can also instruct `anvi-script-gen-genomes-file` to include all subdirectories under a given directory path:

```
anvi-script-gen-genomes-file --input-dir path/to/dir \
                             --output-file external_genomes.txt \
                             --include-subdirs
```

## Generating an internal genomes file

To get an <span class="artifact-n">[internal-genomes](/software/anvio/help/7.1/artifacts/internal-genomes)</span> file containing all bins from a collection, provide a <span class="artifact-n">[profile-db](/software/anvio/help/7.1/artifacts/profile-db)</span>, its corresponding <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>, and the <span class="artifact-n">[collection](/software/anvio/help/7.1/artifacts/collection)</span> name:

<div class="codeblock" markdown="1">
anvi&#45;script&#45;gen&#45;genomes&#45;file &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                             &#45;p <span class="artifact&#45;n">[profile&#45;db](/software/anvio/help/7.1/artifacts/profile&#45;db)</span> \
                             &#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/7.1/artifacts/collection)</span> \
                             &#45;&#45;output&#45;file internal&#45;genomes.txt
</div>

The name of each internal genome will be the same as the bin name, and the path columns will contain absolute paths.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-script-gen-genomes-file.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-gen-genomes-file) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

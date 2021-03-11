---
layout: page
title: anvi-script-gen-genomes-file [program]
categories: [anvio]
comments: false
redirect_from: /m/anvi-script-gen-genomes-file
image:
  featurerelative: ../../../images/header.png
  display: true
---

Generate an external genomes or internal genomes file.

See **[program help menu](../../../../vignette#anvi-script-gen-genomes-file)** or go back to the **[main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Can provide

<p style="text-align: left" markdown="1"><span class="artifact-p">[external-genomes](../../artifacts/external-genomes) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span> <span class="artifact-p">[internal-genomes](../../artifacts/internal-genomes) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>

## Can consume

<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[collection](../../artifacts/collection) <img src="../../images/icons/COLLECTION.png" class="artifact-icon-mini" /></span></p>

## Usage


This script can automatically generate an external or internal genomes file.

## Generating an external genomes file
If you provide an input directory and a name for the output file, then every <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span> in that directory will get a line in the resulting <span class="artifact-n">[external-genomes](/software/anvio/help/main/artifacts/external-genomes)</span> file:

```
anvi-script-gen-genomes-file --input-dir path/to/dir -e external_genomes.txt
```

The name of each database will be whatever string is in front of the `*.db` extension, and the `contigs_db_path` column will contain absolute paths.

## Generating an internal genomes file
To get an <span class="artifact-n">[internal-genomes](/software/anvio/help/main/artifacts/internal-genomes)</span> file containing all bins from a collection, provide a <span class="artifact-n">[profile-db](/software/anvio/help/main/artifacts/profile-db)</span>, its corresponding <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span>, and the <span class="artifact-n">[collection](/software/anvio/help/main/artifacts/collection)</span> name:

```
anvi-script-gen-genomes-file -i internal_genomes.txt -c CONTIGS.db -p PROFILE.db -C default
```

The name of each internal genome will be the same as the bin name, and the path columns will contain absolute paths.


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-script-gen-genomes-file.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-script-gen-genomes-file) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

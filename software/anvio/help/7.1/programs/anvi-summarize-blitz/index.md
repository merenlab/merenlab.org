---
layout: page
title: anvi-summarize-blitz [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-summarize-blitz
image:
  featurerelative: ../../../images/header.png
  display: true
---

FAST summary of many anvi&#x27;o single profile databases (without having to use the program anvi-merge)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/ivagljiva.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Iva Veseli</span><div class="page-author-social-box"><a href="mailto:iveseli@uchicago.edu" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://github.com/ivagljiva" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[single-profile-db](../../artifacts/single-profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span></p>


## Can provide


<p style="text-align: left" markdown="1"><span class="artifact-p">[quick-summary](../../artifacts/quick-summary) <img src="../../images/icons/TXT.png" class="artifact-icon-mini" /></span></p>


## Usage


This program is a quicker, but less comprehensive, alternative to <span class="artifact-n">[anvi-summarize](/software/anvio/help/7.1/programs/anvi-summarize)</span>. It is used to summarize basic read recruitment statistics (like detection and coverage) from many single profiles that are all associated with the same <span class="artifact-n">[contigs-db](/software/anvio/help/7.1/artifacts/contigs-db)</span>.

Given a list of samples (single profiles) and a collection, `anvi-summarize-blitz` will compute the per-sample weighted average of each statistic for each bin in the collection. This is an average of the statistic value over each split in the bin, _weighted by the split length_.

The output will be a text file, and you can find details about its format by clicking on <span class="artifact-n">[quick-summary](/software/anvio/help/7.1/artifacts/quick-summary)</span>.

### Basic usage

In addition to your list of <span class="artifact-n">[single-profile-db](/software/anvio/help/7.1/artifacts/single-profile-db)</span>s, you must provide this program with their corresponding contigs database and a collection name.

<div class="codeblock" markdown="1">
anvi&#45;summarize&#45;blitz PROFILE_1.db PROFILE_2.db PROFILE_3.db [...] \
                     &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                     &#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/7.1/artifacts/collection)</span>
</div>

The program will summarize the same collection across all of your single profile databases. However, it will use only the first profile database in the argument list to learn about what is in the collection, so it is not exactly necessary to have this collection defined for all of the other profile databases (though one could argue that it is a good idea to do this regardless...). The collection name you provide to this program must be a collection that is present in at least the first profile database in the argument list. In the example above, only `PROFILE_1.db` is strictly required to include the collection you wish to summarize (though all other profiles must contain the same splits as this first profile, which should not be a problem if you generated them all in the same way).

### Choosing a different output prefix

If nothing is provided, the output file name will be the collection name, suffixed with `-SUMMARY-BLITZ.txt` (although the user can specify the output file name as they should using the parameter `--output-file`):

<div class="codeblock" markdown="1">
anvi&#45;summarize&#45;blitz PROFILE_1.db PROFILE_2.db PROFILE_3.db [...] \
                     &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                     &#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/7.1/artifacts/collection)</span> \
                     &#45;o OUTPUT.txt
</div>

### Choosing which statistics to summarize

The default statistics that will be summarized are detection and something called 'mean_coverage_Q2Q3' (which is [this](https://merenlab.org/2017/05/08/anvio-views/#mean-overage-q2q3)). You can choose which statistics to summarize by providing them as a comma-separated list (no spaces in the list) to the `--stats-to-summarize`, or `-S`, parameter:

<div class="codeblock" markdown="1">
anvi&#45;summarize&#45;blitz PROFILE_1.db PROFILE_2.db PROFILE_3.db [...] \
                     &#45;c <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/7.1/artifacts/contigs&#45;db)</span> \
                     &#45;C <span class="artifact&#45;n">[collection](/software/anvio/help/7.1/artifacts/collection)</span> \
                     &#45;S std_coverage,mean_coverage,detection
</div>

Each statistic will get its own column in the output file.

If you are not sure which statistics are available to choose from, just provide some ridiculous, arbitrary string (that cannot possibly be a name of a statistic) to this flag, and you will get an error message that includes a list of the available statistics. Or, you can just look at this example error message (but no guarantees that the list in this example will be the same as whatever you would get by doing it yourself. Just sayin'.)
```
Config Error: The statistic you requested, cattywampus, does not exist. Here are the options
              to choose from: std_coverage, mean_coverage, mean_coverage_Q2Q3, detection,
              abundance, variability
```

If you are curious about the statistics in the list, many of them have definitions in [this blog post](https://merenlab.org/2017/05/08/anvio-views).

## Common errors

### Existing file error

If the output file already exists, you will encounter the following error:
```
File/Path Error: AppendableFile class is refusing to open your file at test-quick_summary.txt
                 because it already exists. If you are a user, you should probably give Anvi'o a
                 different file name to work with. If you are a programmer and you don't want
                 this behavior, init this class with `fail_if_file_exists=False` instead.
```
You can either provide a different file prefix using the `-O` parameter, as the error message suggests, or you can simply delete the existing file and re-run your command.

### Missing table error

If you get an error that looks like this:
```
Config Error: The database at [PROFILE.db] does not seem to have a table named
              `detection_splits` :/ Here is a list of table names this database knows:
              [...]
```

That means your profile databases are not the correct version. The tables we are accessing in this program were introduced in profile database version 36. So the solution to this error is to update your databases to at least that version, using <span class="artifact-n">[anvi-migrate](/software/anvio/help/7.1/programs/anvi-migrate)</span>. :)


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-summarize-blitz.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-summarize-blitz) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

---
layout: page
title: anvi-migrate [program]
categories: [anvio]
comments: false
redirect_from: /7.1/anvi-migrate
image:
  featurerelative: ../../../images/header.png
  display: true
---

Migrate an anvi&#x27;o database or config file to a newer version.

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.


{% include _toc.html %}
<div id="svg" class="subnetwork"></div>
{% capture network_path %}{{ "network.json" }}{% endcapture %}
{% capture network_height %}{{ 300 }}{% endcapture %}
{% include _project-anvio-graph.html %}


## Authors

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/meren.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">A. Murat Eren (Meren)</span><div class="page-author-social-box"><a href="http://meren.org" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:a.murat.eren@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/merenbey" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/meren" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/ozcan.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Özcan C. Esen</span><div class="page-author-social-box"><a href="http://blog.ozcanesen.com/" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:ozcanesen@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/ozcanesen" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/ozcan" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/ekiefl.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Evan Kiefl</span><div class="page-author-social-box"><a href="mailto:kiefl.evan@gmail.com" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/evankiefl" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/ekiefl" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/ivagljiva.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Iva Veseli</span><div class="page-author-social-box"><a href="mailto:iveseli@uchicago.edu" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://github.com/ivagljiva" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>

<div class="page-author"><div class="page-author-info"><div class="page-person-photo"><img class="page-person-photo-img" src="../../images/authors/semiller10.jpg" /></div><div class="page-person-info-box"><span class="page-author-name">Samuel Miller</span><div class="page-author-social-box"><a href="https://semiller10.github.io" class="person-social" target="_blank"><i class="fa fa-fw fa-home"></i>Web</a><a href="mailto:samuelmiller@uchicago.edu" class="person-social" target="_blank"><i class="fa fa-fw fa-envelope-square"></i>Email</a><a href="http://twitter.com/smiller_science" class="person-social" target="_blank"><i class="fa fa-fw fa-twitter-square"></i>Twitter</a><a href="http://github.com/semiller10" class="person-social" target="_blank"><i class="fa fa-fw fa-github"></i>Github</a></div></div></div></div>



## Can consume


<p style="text-align: left" markdown="1"><span class="artifact-r">[contigs-db](../../artifacts/contigs-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[profile-db](../../artifacts/profile-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[pan-db](../../artifacts/pan-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[genes-db](../../artifacts/genes-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[genomes-storage-db](../../artifacts/genomes-storage-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[structure-db](../../artifacts/structure-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[modules-db](../../artifacts/modules-db) <img src="../../images/icons/DB.png" class="artifact-icon-mini" /></span> <span class="artifact-r">[workflow-config](../../artifacts/workflow-config) <img src="../../images/icons/JSON.png" class="artifact-icon-mini" /></span></p>


## Can provide


This program does not seem to provide any artifacts. Such programs usually print out some information for you to see or alter some anvi'o artifacts without producing any immediate outputs.


## Usage


This is a multi-talented program that seamlessly updates any anvi'o database to the latest version.

You can provide one or more anvi'o databases as command line parameters to this program, and it will migrate each one. However, you must choose whether to migrate the databases safely, or quickly.

If you choose to migrate safely, anvi'o will first make a copy of each database and save it as a backup. In case something goes wrong during the migration, it let you know what happened and will restore your original database from the copy it made. Then you can go on your merry way. (The copy is deleted after the migration script is finished running.)

This is how you migrate safely:
```
anvi-migrate --migrate-dbs-safely *.db
```

Of course, we will always suggest that migrating safely is better, because fewer people get angry at us when we do that. In practice though, making those backup copies takes up extra time and it is unlikely that the migration will fail anyway, so if you have a lot of databases to migrate and are okay with a bit of risk, you have the option to migrate quickly instead. In this case, anvi'o will _not_ copy your databases before starting the migration.
```
anvi-migrate --migrate-dbs-quickly *.db
```
Please remember that by living life in the fast lane, you forego your safety net. On the rare occasion that the migration does fail, this program will let you know what happened, leave you with a database that has a `.broken` file extension, and sassily remind you that this all could have been avoided if you chose the other option. In this unlikely event, you can always reach out to us. We will probably be sassy to you, too, but we will still see if we can help you unbreak things. :)

### Migrating to a specific version
If your database is a few versions behind the highest available version but for whatever reason you don't want to migrate it all the way, you can specify which version to update your database to. Just use the `-t` flag (note: migrating with this parameter only works on ONE database at a time):
```
anvi-migrate --migrate-dbs-safely -t 15 CONTIGS.db
```
Then anvi'o will update your database until it is whatever version you specified and stop. Of course, you cannot provide a version number that is higher than the highest available version. Nor can you provide a number that is lower than your database's current version (ie, backwards migration is not possible).

Not sure what your database's current version is? Try <span class="artifact-n">[anvi-db-info](/software/anvio/help/7.1/programs/anvi-db-info)</span>.
Not sure what the highest available version is? Run any anvi'o command with the `-v` option to see the version information for all database types (we recommend `anvi-interactive -v` for no particular reason).


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/programs/anvi-migrate.md) to update this information.


## Additional Resources



{:.notice}
Are you aware of resources that may help users better understand the utility of this program? Please feel free to edit [this file](https://github.com/merenlab/anvio/tree/master/bin/anvi-migrate) on GitHub. If you are not sure how to do that, find the `__resources__` tag in [this file](https://github.com/merenlab/anvio/blob/master/bin/anvi-interactive) to see an example.

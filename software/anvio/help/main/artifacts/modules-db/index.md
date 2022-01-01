---
layout: page
title: modules-db [artifact]
categories: [anvio]
comments: false
redirect_from: /m/modules-db
image:
  featurerelative: ../../../images/header.png
  display: true
---


{% include _toc.html %}


<img src="../../images/icons/DB.png" alt="DB" style="width:100px; border:none" />

A DB-type anvi'o artifact. This artifact is typically generated, used, and/or exported **by anvi'o** (and not provided by the user)..

🔙 **[To the main page](../../)** of anvi'o programs and artifacts.

## Provided by


<p style="text-align: left" markdown="1"><span class="artifact-p">[anvi-setup-kegg-kofams](../../programs/anvi-setup-kegg-kofams)</span> <span class="artifact-p">[anvi-setup-user-modules](../../programs/anvi-setup-user-modules)</span></p>


## Required or used by


<p style="text-align: left" markdown="1"><span class="artifact-r">[anvi-migrate](../../programs/anvi-migrate)</span></p>


## Description

A type of database containing information from either A) the [KEGG MODULE database](https://www.genome.jp/kegg/module.html), or B) user-defined metabolic modules, for use in metabolism estimation and/or functional annotation of KEGG Orthologs (KOs).

These databases are part of the <span class="artifact-n">[kegg-data](/software/anvio/help/main/artifacts/kegg-data)</span> and <span class="artifact-n">[user-modules-data](/software/anvio/help/main/artifacts/user-modules-data)</span> directories. You can get one on your computer by running <span class="artifact-n">[anvi-setup-kegg-kofams](/software/anvio/help/main/programs/anvi-setup-kegg-kofams)</span> or <span class="artifact-n">[anvi-setup-user-modules](/software/anvio/help/main/programs/anvi-setup-user-modules)</span>. Programs that rely on this type of database include <span class="artifact-n">[anvi-run-kegg-kofams](/software/anvio/help/main/programs/anvi-run-kegg-kofams)</span> and <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/main/programs/anvi-estimate-metabolism)</span>.

Most users will never have to interact directly with this kind of database. However, for the brave few who want to try this (or who are figuring out how anvi'o works under the hood), there is some relevant information below.

## Database Contents

### The modules table

In the current implementation, data about each metabolic pathway (from the KEGG MODULE database, or from user-defined modules) is present in the `modules` table, which looks like this:

| module | data_name | data_value | data_definition | line |
|:--|:--|:--|:--|:--|
| M00001 | ENTRY	| M00001 | Pathway | 1 |
| M00001 | NAME	| Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate | _NULL_ | 2 |
| M00001 | DEFINITION | (K00844,K12407,K00845,K00886,K08074,K00918) (K01810,K06859,K13810,K15916) (K00850,K16370,K21071,K00918) (K01623,K01624,K11645,K16305,K16306) K01803 ((K00134,K00150) K00927,K11389) (K01834,K15633,K15634,K15635) K01689 (K00873,K12406) | _NULL_ | 3 |
| M00001 | ORTHOLOGY | K00844	| hexokinase/glucokinase [EC:2.7.1.1 2.7.1.2] [RN:R01786] | 4 |
| M00001 | ORTHOLOGY | K12407	| hexokinase/glucokinase [EC:2.7.1.1 2.7.1.2] [RN:R01786] | 4 |
| (...) | (...) | (...) | (...) | (...) |

For the MODULES.db that comes out of <span class="artifact-n">[anvi-setup-kegg-kofams](/software/anvio/help/main/programs/anvi-setup-kegg-kofams)</span>, these data correspond to the information that can be found on the KEGG website for each metabolic module - for an example, you can see the page for [M00001](https://www.genome.jp/dbget-bin/www_bget?md:M00001) (or, alternatively, its [flat text file version](http://rest.kegg.jp/get/M00001) from the KEGG REST API).

The USER_MODULES.db that comes out of <span class="artifact-n">[anvi-setup-user-modules](/software/anvio/help/main/programs/anvi-setup-user-modules)</span> contains similar information, but defined by the user instead of downloaded from the KEGG website.

In either case, the `module` column indicates the module ID number while the `data_name` column indicates what type of data the row is describing about the module. These data names are usually fairly self-explanatory - for instance, the `DEFINITION` rows describe the module definition and the `ORTHOLOGY` rows describe the enzymes belonging to the module - however, for an official explanation, you can check [the KEGG help page](https://www.genome.jp/kegg/document/help_bget_module.html).

The `data_value` and `data_definition` columns hold the information corresponding to the row's `data_name`; for `ORTHOLOGY` fields these are the enzyme accession number and its functional annotation, respectively. Not all rows have a `data_definition` field.

Finally, some rows of data originate from the same line in the original KEGG MODULE text file; these rows will have the same number in the `line` column. Perhaps this is a useless field. But it is there.

### The database hash value

In the `self` table of this database, there is an entry called `hash`. This string is a hash of the contents of the database, and it allows us to identify the version of the data within the database. This value is important for ensuring that the same MODULES.db is used both for annotating a contigs database with <span class="artifact-n">[anvi-run-kegg-kofams](/software/anvio/help/main/programs/anvi-run-kegg-kofams)</span> and for estimating metabolism on that contigs database with <span class="artifact-n">[anvi-estimate-metabolism](/software/anvio/help/main/programs/anvi-estimate-metabolism)</span>.

You can easily check the hash value by running the following:

<div class="codeblock" markdown="1">
anvi&#45;db&#45;info <span class="artifact&#45;n">[modules&#45;db](/software/anvio/help/main/artifacts/modules&#45;db)</span>
</div>

It will appear in the `DB Info` section of the output, like so:
```
DB Info (no touch also)
===============================================
num_modules ..................................: 443
total_entries ................................: 13720
creation_date ................................: 1608740335.30248
hash .........................................: 45b7cc2e4fdc
```

If you have annotated a <span class="artifact-n">[contigs-db](/software/anvio/help/main/artifacts/contigs-db)</span> using <span class="artifact-n">[anvi-run-kegg-kofams](/software/anvio/help/main/programs/anvi-run-kegg-kofams)</span>, you would find that the corresponding hash in that contigs database matches to this one:

<div class="codeblock" markdown="1">
anvi&#45;db&#45;info <span class="artifact&#45;n">[contigs&#45;db](/software/anvio/help/main/artifacts/contigs&#45;db)</span>
</div>

```
DB Info (no touch also)
===============================================
[....]
modules_db_hash ..............................: 45b7cc2e4fdc
```

### Other important values in the self table

The `data_source` key will tell you if the current database was generated from KEGG data using <span class="artifact-n">[anvi-setup-kegg-kofams](/software/anvio/help/main/programs/anvi-setup-kegg-kofams)</span> or from user-defined metabolic modules using <span class="artifact-n">[anvi-setup-user-modules](/software/anvio/help/main/programs/anvi-setup-user-modules)</span>.

The `annotation_sources` key will list the functional annotation sources that are required to annotate all enzymes found in the module definitions.

Here is an example of what these fields look like for a KEGG MODULES.db:
```
DB Info (no touch also)
===============================================
data_source ..................................: KEGG
annotation_sources ...........................: KOfam
```

And here is an example of what they look like for a USER_MODULES.db:
```
DB Info (no touch also)
===============================================
data_source ..................................: USER
annotation_sources ...........................: KOfam,UpxZ,COG20_FUNCTION
```

## Querying the database

If you want to extract information directly from a modules database, you can do it with a bit of SQL :)

Here is one example, which obtains the name of every module in the default KEGG database:

```
# learn where the MODULES.db is:
export ANVIO_MODULES_DB=`python -c "import anvio; import os; print(os.path.join(os.path.dirname(anvio.__file__), 'data/misc/KEGG/MODULES.db'))"`
# get module names:
sqlite3 $ANVIO_MODULES_DB "select module,data_value from modules where data_name='NAME'" | \
    tr '|' '\t' > module_names.txt
```

## Loading the database in Python

The modules database class has plenty of helpful functions defined for it. You can easily load one in Python and use these functions to access the data within. Here is how you load the database:

```python
import anvio
import argparse
import os
from anvio import kegg

args = argparse.Namespace()
# CHANGE THIS PATH IF YOU WANT TO LOAD A MODULES DB AT A NON-DEFAULT LOCATION
path_to_db = os.path.join(os.path.dirname(anvio.__file__), 'data/misc/KEGG/MODULES.db')
db = kegg.KeggModulesDatabase(path_to_db, args)
```
Once you have done this, you can start to use the helper functions. For example, the following function will return a list of all paths through a module:
```python
db.unroll_module_definition('M00001')
```


{:.notice}
Edit [this file](https://github.com/merenlab/anvio/tree/master/anvio/docs/artifacts/modules-db.md) to update this information.


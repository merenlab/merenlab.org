---
layout: post
title: Working with item additonal data tables
modified: 2017-12-11
excerpt: "Getting additional data in an out of pan and profile databases like a pro."
comments: true
authors: [meren]
categories: [anvio]
---

{% include _toc.html %}

{% capture images %}{{site.url}}/images/anvio/2017-12-11-items-additional-data-tables{% endcapture %}

As you know, the anvi'o interactive interface is quite flexible, and allows you to add many different kinds of data into a single display.

Unless you are using the interactive interface in [`--manual` mode]({{ site.url }}/tutorials/interactive-interface/){:target="_blank"}) (by providing every bit of data yourself in an *ad hoc* manner), most of the data displayed in the interface comes from pre-computed information stored in anvi'o contigs, profile, or samples databases (such as GC-content of contigs, or coverage values from samples, or number of reads mapped to a given sample, etc).

Depending on the project, however, we often want to add additional stuff into these displays. For instance the display below comes from [one of our publications](https://peerj.com/articles/1839/), where we had to add two additional layers of information:

[![Tardigrade](https://i.imgur.com/3vhKxcI.jpg){:.center-img .width-70}](https://peerj.com/articles/1839/#fig-1){:target="_blank"}

In this display, every data point in this display comes from the anvi'o profile and contigs database for this project except the last two layers under 'Genomic selections' layer. These project-specific information comes in the form of an *additional data file* we generated separately, and provided to the `anvi-interactive` with the `--additional-layers` parameter to see them in the display. If you want, you can take a look at this [tutorial on anvi'o interactive interface]({{ site.url }}/tutorials/interactive-interface/#additional-data-for-the-items){:target="_blank"} to get more information on how to extend anvi'o displays in similar ways, but we are here because we think those are not enough.

We know it is helpful to be able to extend any anvi'o display in any direction by adding new information into it in an *ad hoc* manner without much pain using a TAB-delimited file through the `--additional-layers` parameter. But when the only option to do is through this flag, two bad things happen: one, it becomes necessary to carry around the additional data files with profile databases to make things reproducible, and two, the lack of a standard, programmable way to add or remove such additional information to databases prevents the implementation of cool ideas that can use this to generate insights to display smoothly.

That's why, while keeping the *ad hoc* workflow in place, we recently added a way to work with item additional data, which will be available with anvi'o `v4`.

## Dealing with item additional data tables as a user

As a user, you can add, remove, export, or update any data into your pan and profile databases through a couple of command-line programs: **anvi-import-item-additional-data**, **anvi-export-item-additional-data**, and **anvi-delete-item-additional-data**.

Here is an example. If you wish to follow the example on your own anvi'o `v4` installed computer, you can do this:

``` bash
 $ wget http://merenlab.org/files/items_additional_data_example.tar.gz
 $ tar -zxvf items_additional_data_example.tar.gz
 $ cd items_additional_data_example/
 $ ls
additional_data.txt tree.txt            view_data.txt
```

We could visualize the contents of the `view_data.txt` given the `tree.txt` file in manual mode the following way:

``` bash
anvi-interactive -d view_data.txt \
                 -t tree.txt \
                 -p profile.db \
                 --title "Test" \
                 --manual
```

And clicking 'Draw' would have given us this:

[![image]({{images}}/01.png)]({{images}}/01.png){:.center-img .width-70}

There is also an additional data file in the same directory that goes like this:

|contig|categorical_1|categorical_2|text_layer_01|numerical|bars_main!A;B;C|
|:--|:--:|:--:|:--:|:--:|:--:|
|backrest|b|y|nmwje|2.78|278;23;1|
|backward|b|x|bqmyujr psrd doefhi|2.49|249;52;2|
|backwind|b|y|hkfer lchpmzix|2.69|269;32;3|
|backyard|b|x|advoe bfkyhmg|2.05|205;96;4|
|bacteria|b|x|lqmcwn hywco|2.63|263;38;5|
|bacterin|b||vxqdmn|2.98|298;3;6|
|baetylus|b|x|fkgpydi owgyhfx xwlpj|2.19|219;82;7|
|bagpiped|b|y|ijmnur|2.12|212;89;8|
|balconet|b|y|ecizgs|2.89|289;12;9|
|(...)|(...)|(...)|(...)|(...)|(...)|

This file can be imported into the profile database the following way:

``` bash
$ anvi-import-item-additional-data additional_data.txt -p profile.db

New additional data...
===============================================
Key "categorical_1" ..........................: Predicted type: str
Key "categorical_2" ..........................: Predicted type: str
Key "text_layer_01" ..........................: Predicted type: str
Key "numerical" ..............................: Predicted type: float
Key "bars_main!A;B;C" ........................: Predicted type: stackedbar

WARNING
===============================================
You (or the programmer) asked anvi'o to NOT check the consistency of item names
between your additional data and the profile database you are attempting to
update. So be it.

New data added to the db .....................: categorical_1, categorical_2, text_layer_01, numerical, bars_main!A;B;C.
```

Now running the interactive interface again will give us something extra:

``` bash
anvi-interactive -d view_data.txt \
                 -t tree.txt \
                 -p profile.db \
                 --title "Test" \
                 --manual
```

[![image]({{images}}/02.png)]({{images}}/02.png){:.center-img .width-70}

All the layers in the additional data file appears in the same order in the interface. Fine.

You can export item additional data from any profile database:

``` bash
$ anvi-export-item-additional-data -p profile.db \
                                   -o exported_additional_file.txt
                                   
Output file ..................................: exported_additional_file.txt
```

Or you can delete the contents of this table as a whole, or by using a specific data key. If you don't know the available keys, you can simply ask anvi'o:

``` bash
 $ anvi-delete-item-additional-data -p profile.db \
                                    --list-available-keys

AVAILABLE DATA KEYS (5 FOUND)
===============================================
* bars_main!A;B;C (stackedbar, describes 300 items)
* categorical_1 (str, describes 300 items)
* categorical_2 (str, describes 300 items)
* numerical (float, describes 300 items)
* text_layer_01 (str, describes 300 items)

```

And then delete one of those:

``` bash
 $ anvi-delete-item-additional-data -p profile.db \
                                    --keys-to-remove categorical_1,categorical_2

WARNING
===============================================
Data for the following keys removed from the
database: 'categorical_1, categorical_2'.

 $ anvi-delete-item-additional-data -p profile.db \
                                    --list-available-keys

AVAILABLE DATA KEYS (3 FOUND)
===============================================
* bars_main!A;B;C (stackedbar, describes 300 items)
* numerical (float, describes 300 items)
* text_layer_01 (str, describes 300 items)

```

Or you can simply delete everything:

``` bash
$ anvi-delete-item-additional-data -p profile.db

WARNING
===============================================
All data from the items additional data table is removed.

 $ anvi-delete-item-additional-data -p profile.db --list-available-keys
 
* There are no item additional data in this database.

```

So that's that.


## Dealing with item additional data tables as a programmer

If you are writing a Python program, you can simply deal with additional data items the following way:

``` python
import argparse
import anvio.dbops as dbops

args = argparse.Namespace(pan_or_profile_db="/path/to/profile.db")

item_additional_data_table = dbops.TableForItemAdditionalData(args)

# add data:
item_additional_data_table.add(new_keys_list, new_data_dict)

# read data:
items_additional_data_keys, items_additional_data_dict = items_additional_data.get()

# remove data:
item_additional_data_table.remove(keys_list)
```

For instance, this code should work perfectly in the directory above:

``` python
import argparse

import anvio.dbops as dbops
import anvio.utils as utils

args = argparse.Namespace(profile_db="profile.db")

item_additional_data_table = dbops.TableForItemAdditionalData(args)

# make sense of the additional data file:
keys = utils.get_columns_of_TAB_delim_file('additional_data.txt')
data = utils.get_TAB_delimited_file_as_dictionary('additional_data.txt')

# add data (because this is a blank profile, we use `skip_check_names`
# flag. This flag is not necessary when working with regular profile
# and pan databases:
item_additional_data_table.add(keys, data, skip_check_names=True)

# get data from the database:
keys, data = item_additional_data_table.get()

# remove data:
item_additional_data_table.remove(keys)
```

That's it!

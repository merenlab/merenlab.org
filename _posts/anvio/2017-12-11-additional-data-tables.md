---
layout: post
title: Working with anvi'o additional data tables
modified: 2017-12-11
excerpt: "Getting additional data in an out of pan and profile databases like a pro."
comments: true
authors: [meren]
categories: [anvio]
redirect_from: /2015/11/10/samples-db/
---

{% include _toc.html %}

{% capture images %}{{site.url}}/images/anvio/2017-12-11-additional-data-tables{% endcapture %}

{:.notice}
This post will only work for anvi'o `v4` or later.

As you know, the anvi'o interactive interface is quite flexible, and allows you to add many different kinds of data into a single display.

Unless you are using the interactive interface in [`--manual` mode]({{ site.url }}/tutorials/interactive-interface/){:target="_blank"}) (by providing every bit of data yourself in an *ad hoc* manner), most of the data displayed in the interface comes from pre-computed information stored in anvi'o contigs, profile, or samples databases (such as GC-content of contigs, or coverage values from samples, or number of reads mapped to a given sample, etc).

Depending on the project, however, we often want to add additional stuff into these displays. The program `anvi-interactive` [allows]({{ site.url }}/tutorials/interactive-interface/#additional-data-for-the-items){:target="_blank"} its users to enrich their display with TAB-delimited data files rather quickly, but we are here because we think that flexibility is not enough.

While it is helpful to be able to extend any anvi'o display in any direction by adding new information through TAB-delimited files with `--additional-layers` parameter, this practice requires the users to carry around the additional data files with profile databases to make things reproducible. Plus, the lack of a user and programmer-friendly way to add or remove such additional information to anvi'o pan and profile databases complicates the implementation of cool ideas.

That's why, while keeping the *ad hoc* workflow in place, we have extended anvi'o in `v4` with a completely new design to work with additional data.

{:.notice}
As a part of this design, we ended up killing the [samples database]({% post_url anvio/2015-11-10-samples-db %}){:target="_blank"}. If you are too young to remember those days, that's fine, you will not miss the samples database. If you are stuck with a an anvio'o project that has a samples database, don't worry, `anvi-migrade-db` program will help you to import the data in the samples database into the new tables in the profile database with no effort.

## Views, items, layers, orders: some anvi'o terminology 

We are as confused as you are, so let's start with this example display to explain what is what:

![Example anvi'o display](https://anvi-server.org/static/img/example.png){:.center-img .width-40}

In this display, the deondrogram denoted by **(1)** shows the organization of **items**. So anything that appears there are what we call 'items'. Concentric circles identified by **(2)** represent the view data. View data is often computed by anvi'o itself, and stored in pan or profile databases. Additional data for items will decorate things around view data display (such as those green things in this particular example). The dendrogram identified by **(3)** shows how those concentric circles, or **layers** should be organized. Data we will import as 'layer orders' will appear there. What is shown by **(4)** is the additional data for layers, and you will also learn in this post how to import data to appear there.

Basically, the purpose of this post is to show you how to annotate a display with additional data for items, layers, and layer orders. For all these tasks, we will use the same three programs, **anvi-import-misc-data**, **anvi-export-misc-data**, and **anvi-delete-misc-data** with different target tables (such as `items` to decorate items, `layers` to make **(4)** appear, or `layer_orders` make data for **(3)** available for our pan or profile databases.

Throughout this post, I will use a simple dataset for demonstration purposes. If you would like to follow it on your anvi'o `v4` or later installed computer, first run these commands:


``` bash
wget http://merenlab.org/files/anvio_additional_data_tables_example.tar.gz
tar -zxvf anvio_additional_data_tables_example.tar.gz
cd anvio_additional_data_tables_example/
```


## Dealing with data tables as a user

This section will show step by step which table is good for what in the interface. Although these examples will use a blank anvi'o profile database, everything will work the same way for regular profile databases and pan databases.

### Items additional data table

**This is the table you want to work with if you would like to show things for each of your contigs, gene clusters, or any other item you have in the center. The target table name for items is `items`.**

Let's start by visualizing the contents of the `view_data.txt` given the `tree.txt` file in manual mode the following way:

``` bash
anvi-interactive -d view_data.txt \
                 -t tree.txt \
                 -p profile.db \
                 --title "Test" \
                 --manual
```

And clicking 'Draw' would have given us this:

[![image]({{images}}/01.png)]({{images}}/01.png){:.center-img .width-70}

There isn't much to look at. Fine. Let's assume we have the following information for each item displayed in here which goes like this:

|item_name|categorical_1|categorical_2|text_layer_01|numerical|bars_main!A;B;C|
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

If you take a look at the `view_data.txt` you would realize that the first columns of both files are pretty identical. 

This additional data shown for items can be imported into the profile database the following way:

``` bash
$ anvi-import-misc-data items_additional_data.txt \
                              -p profile.db \
                              --target-data-table items

New items additional data...
===============================================
Key "categorical_1" ..........................: Predicted type: str
Key "categorical_2" ..........................: Predicted type: str
Key "text_layer_01" ..........................: Predicted type: str
Key "numerical" ..............................: Predicted type: float
Key "bars_main!A;B;C" ........................: Predicted type: stackedbar


WARNING
===============================================
You (or the programmer) asked anvi'o to NOT check the consistency of the names
of your items between your additional data and the profile database you are
attempting to update. So be it.

New data added to the db for your items ......: categorical_1, categorical_2, text_layer_01, numerical, bars_main!A;B;C.
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

All the columns in the additional data file appears in the same order in the interface. Plus, if you click open the Settings panel, you can see the items order combo box in the main panel is already automatically populated with new orders of your items based on this new additional data:

[![image]({{images}}/02.1.png)]({{images}}/02.1..png){:.center-img .width-50}

{:.notice}
While we are here, let's save a default state by clicking Save State so next rounds we don't have to click Draw.

You can export data from a given table:

``` bash
$ anvi-export-misc-data -p profile.db \
                        --target-data-table items \
                        -o exported_additional_file.txt
                                   
Output file for items ........................: exported_additional_file.txt
```

Or you can delete the contents of a given table as a whole, or only a specific data items by specifying one or more data keys. Anvi'o can tell you about the available keys in a given table with `--list-available-keys` flag:

``` bash
 $ anvi-delete-misc-data -p profile.db \
                         --target-data-table items \
                         --list-available-keys

AVAILABLE DATA KEYS FOR ITEMS (5 FOUND)
===============================================
* bars_main!A;B;C (stackedbar, describes 300 items)
* categorical_1 (str, describes 300 items)
* categorical_2 (str, describes 300 items)
* numerical (float, describes 300 items)
* text_layer_01 (str, describes 300 items)

```

Say we delete one of those:

``` bash
 $ anvi-delete-misc-data -p profile.db \
                         --target-data-table items \
                         --keys-to-remove categorical_1,categorical_2

WARNING
===============================================
items data for the following keys removed from the database: 'categorical_1,
categorical_2'.
```

And now it's gone:

``` bash

 $ anvi-delete-misc-data -p profile.db \
                         --target-data-table items \
                         --list-available-keys

AVAILABLE DATA KEYS FOR ITEMS (3 FOUND)
===============================================
* bars_main!A;B;C (stackedbar, describes 300 items)
* numerical (float, describes 300 items)
* text_layer_01 (str, describes 300 items)

```

If you do not specify a data key, the entire content of the table would go bye bye:

``` bash
$ anvi-delete-misc-data -p profile.db \
                        --target-data-table items


WARNING
===============================================
All data from the items additional data table is removed.
```

And then you would see nothing:

``` bash
 $ anvi-delete-misc-data -p profile.db \
                         --target-data-table items \
                         --list-available-keys
 
* There are no item additional data for items in this database.

```

To continue with some data in our items additional data table, let's repopulate it and continue with the layer add:

``` bash
$ anvi-import-misc-data items_additional_data.txt \
                        -p profile.db \
                        --target-data-table items
```


### Layers additional data table

**This is the table you want to work with if you would like to show things for each of your layers. These layers could be your metagenomic, metatranscriptomic samples, or your genomes or any other layer identified as (2) in the figure above. The target table name for layers is `layers`.**


Access to layer additional data tables is conceptually identical to the way we work with additional data for items, and it requires a small change in the command line. For instance, take the following file:

|samples|numerical_01|numerical_02|categorical|stacked_bar!X;Y;Z|
|:--|:--:|:--:|:--:|:--:|
|c1|100|5|A|1;2;3|
|c2|200|4|B|2;3;1|
|c3|300|3|B|3;1;2|

Now the first column of this file is identical to our layer names, and every column describes a property of a given layer.

We could add this into the profile database this way :

```
 $ anvi-import-misc-data layers_additional_data.txt \
                         -p profile.db \
                         --target-data-table layers

New layers additional data...
===============================================
Key "numerical_01" ...........................: Predicted type: int
Key "numerical_02" ...........................: Predicted type: int
Key "categorical" ............................: Predicted type: str
Key "stacked_bar!X;Y;Z" ......................: Predicted type: stackedbar

New data added to the db for your layers .....: numerical_01, numerical_02, categorical, stacked_bar!X;Y;Z.
```

After this, more information for each layer should show up on the right hand side when you re-run the interactive interface:

``` bash
anvi-interactive -d view_data.txt \
                 -t tree.txt \
                 -p profile.db \
                 --title "Test" \
                 --manual
```

If you run the interactive interface again, you should see a new addition to the display:

[![image]({{images}}/03.png)]({{images}}/03.png){:.center-img .width-70}

The layer additional data in the input file is displayed with the same order they appeared in the file. In fact, if you click open the settings panel, and switch to the Samples tab, you can see that the combo box for sample orders is already populated with some automatic orders to organize your layers based on these data:

[![image]({{images}}/04.png)]({{images}}/04.png){:.center-img .width-50}


### Layer orders additional data table

**This is the table you want to work with if you would like to store specific orderings of your layers, such as phylogenetic trees, or orders in basic form. What you can do this table corresponds to the part identified as (3) in the example figure shown at the beginning of this post. The target table name for layers is `layer_orders`.**

The file format for layer orders data is this:

|item_name|data_type|data_value|
|:--|:--:|:--|
|test_tree|newick|(c2:0.0370199,(c1:0.0227268,c3:0.0227268)Int3:0.0370199);|
|test_list|basic|c3,c2,c1|
|(...)|(...)|(...)|

Each layer order could be either in basic or newick form, and you may have as many of those in a layer orders file as you like, of course. When you import a layer orders file the following way:

```bash
 $ anvi-import-misc-data layers_order.txt \
                         -p profile.db \
                         --target-data-table layer_orders

New layer_orders data...
===============================================
Data key "test_tree" .........................: Type: newick
Data key "test_list" .........................: Type: basic

New order data added to the db for layer_orders : test_tree, test_list.
```

Visualize it again,

``` bash
anvi-interactive -d view_data.txt \
                 -t tree.txt \
                 -p profile.db \
                 --title "Test" \
                 --manual
```

Now if you click open the settings panel again, and switch to the Samples tab, you can see your new orders in the combo box for sample orders:

[![image]({{images}}/05.png)]({{images}}/05.png){:.center-img .width-50}

Selecting a tree type order, and re-drawing the display will show your dendrogram on the side:

[![image]({{images}}/06.png)]({{images}}/06.png){:.center-img .width-70}

Done! Now you know how to extend anvi'o interactive interface displays!


### Dealing with item additional data tables as a programmer

If you are writing a Python program, you can simply deal with additional data items the following way:

``` python
import argparse
import anvio.dbops as dbops

args = argparse.Namespace(pan_or_profile_db="/path/to/profile.db", target_data_table="items")

# MiscDataTableFactory will give you the right object based on your `target_data_table`
# argument --alternatively you can directly access to the relevant class. The factory
# pattern makes it easier to seamlessly select the right inheritance route.
items_additional_data_table = dbops.MiscDataTableFactory(args)

# add data:
item_additional_data_table.add(new_keys_list, new_data_dict)

# read data:
items_additional_data_keys, items_additional_data_dict = item_additional_data_table.get()

# remove data:
item_additional_data_table.remove(keys_list)
```

For instance, this code should work perfectly in the directory above:

``` python
import argparse

import anvio.dbops as dbops
import anvio.utils as utils

args = argparse.Namespace(profile_db="profile.db", target_data_table="items")

# get some data:
keys = utils.get_columns_of_TAB_delim_file('items_additional_data.txt')
data = utils.get_TAB_delimited_file_as_dictionary('items_additional_data.txt')

# add it to the database (because this is a blank profile, we use `skip_check_names`
# flag. This flag is not necessary when working with regular profile
# and pan databases.
dbops.MiscDataTableFactory(args).add(data, keys, skip_check_names=True)

# get data from the database:
keys, data = dbops.MiscDataTableFactory(args).get()

# remove data:
dbops.MiscDataTableFactory(args).remove(keys)
```

That's it! By just changing the `target_data_table` variable between `items`, `layers`, or `layer_orders`, you can work with different tables. That said, the following would have given identical results to the one above, but the table selection would be explicitly done rather than through the `target_data_table` argument:

``` python
import argparse

import anvio.dbops as dbops
import anvio.utils as utils

args = argparse.Namespace(profile_db="profile.db")

# get some data:
keys = utils.get_columns_of_TAB_delim_file('items_additional_data.txt')
data = utils.get_TAB_delimited_file_as_dictionary('items_additional_data.txt')

# add data
dbops.TableForItemAdditionalData(args).add(data, keys, skip_check_names=True)

# get data from the database:
keys, data = dbops.TableForItemAdditionalData(args).get()

# remove data:
dbops.TableForItemAdditionalData(args).remove(keys)
```

Please feel free to ask any questions.
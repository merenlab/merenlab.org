Don't let this scare you! It is very simple. It is all copy-paste, and should not take less than 20 minutes if you have a decent internet connection.

# Installing and setting up centrifuge

Decide where you want to keep centrifuge and related data:

```
export CENTRIFUGE="/path/to/a/directory"
```

You can also put this in your `~/.bashrc` or `~/.bash_profile` to make sure your terminal always knows whre it is.


```
cd $CENTFIGURE
git clone https://github.com/infphilo/centrifuge
cd centrifuge && make
```

This compiles everything, but does not install anything. To make sure binary files are available directly, you can run this (and again, you can add this line into your `~/.bashrc` or `~/.bash_profile` to make sure every new terminal session remembers where they are):

```
export PATH=$PATH:$CENTRIFUGE/centrifuge
```

If everything is alright so far, this is what you should see if you run the following command:

```
$ centrifuge --version | head -n 1
centrifuge-class version v1.0.1-beta-27-g30e3f06ec3
```

Good? Good.

Next, you will need to download pre-computed indexes (unless you want to go Voldemort and compile your own indexes). The compressed indexes for Bacteria, Viruses, Human genome is 6.3 Gb, and it will take about 9 Gb on your disk uncompressed. You will download these only for once:

```
$ cd $CENTRIFUGE
$ wget ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/b+h+v.tar.gz
$ tar -zxvf b+h+v.tar.gz && rm -rf b+h+v.tar.gz
$ centrifuge-inspect --name-table $CENTRIFUGE/b+h+v/b+h+v > b+h+v-names-table.txt
```

The following command will make sure you have the necessary conversion table ready for anvi'o when the time comes. This is also something you will have to do only once:

```
centrifuge-inspect --name-table $CENTRIFUGE/b+h+v/b+h+v > $CENTRIFUGE/b+h+v-names-table.txt
```

If everything went alright, you should see something like this when you run the following command:

```
$ wc -l $CENTRIFUGE/b+h+v-names-table.txt
   13455 $CENTRIFUGE/b+h+v-names-table.txt
```

See? You are totally doing this!


# Importing centrifuge results into anvi'o

If you have a working centrifuge setup, the rest should be pretty easy.

Assuming you generated an anvi'o contigs database. To import taxonomy into this contigs database, first you will export all gene calls: 

```
anvi-get-sequences-for-gene-calls -c CONTIGS.db -o anvio-gene-calls.fa
```

Then you will run the following command:

```
$ centrifuge -f -x $CENTRIFUGE/b+h+v/b+h+v anvio-gene-calls.fa -S centrifuge-classification-report.txt
```

This takes about one minute on my laptop for 40,000 genes. If everything worked alright, at this point you have everything you need to update your contigs database! If you get a **command not found** error, you need to go back to the setup phase, and run those two `export` commands in your terminal. 


# Resources

https://github.com/infphilo/centrifuge

http://www.ccb.jhu.edu/software/centrifuge/

http://biorxiv.org/content/biorxiv/early/2016/05/25/054965.full.pdf
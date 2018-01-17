# How to run a copy of this website locally

You can run a local version of this website to see how your changes will look like when they are pushed. For this, you need to install `jekyll`, which can be doen the following way on a Mac OSX:

```
sudo gem install jekyll
```

If the command `gem` is not found, or your Ruby version is outdated, try to install it the following way first:

```
brew install ruby
```

Once the installation is complete, you can get a copy of the web site from GitHub:

```
git clone https://github.com/merenlab/web.git
cd web
```

To make sure you have the necesary gems, run:

```
bundle install
```

And run the following command in this directory:

```
bundle exec jekyll serve --incremental
```

You should see a similar output in your terminal:

```
Configuration file: /Users/evan/Software/web/_config.yml
Configuration file: /Users/evan/Software/web/_config.yml
            Source: /Users/evan/Software/web
       Destination: /Users/evan/Software/web/_site
 Incremental build: enabled
      Generating...
                    done in 6.585 seconds.
 Auto-regeneration: enabled for '/Users/evan/Software/web'
Configuration file: /Users/evan/Software/web/_config.yml
    Server address: http://127.0.0.1:4000/
  Server running... press ctrl-c to stop.
```

This basically runs a local server for you to see the changes you've made locally. You can access to this server by visiting the URL `http://127.0.0.1:4000/` with your browser.

With the `--incremental` flag every change you will make in any of the files will be reflected to your local website automatically.

# How to get yourself fired from the lab

It is simple.

You can make the following changes in the web site:

```
find ./_posts/ -type f | xargs sed -i '' -e  "s/An anvi'o/A microbial dark matter \(MDM\)/g"
find ./_posts/ -type f | xargs sed -i '' -e  "s/an anvi'o/a microbial dark matter \(MDM\)/g"
find ./_posts/ -type f | xargs sed -i '' -e  "s/Anvi'o/Microbial dark matter \(MDM\)/g"
find ./_posts/ -type f | xargs sed -i '' -e  "s/anvi'o/microbial dark matter \(MDM\)/g"
find ./_posts/ -type f | xargs sed -i '' -e  "s/anvi'o/microbial dark matter \(MDM\)/g"
find ./_posts/ -type f | xargs sed -i '' -e  "s/anvi'server/MDM'server/g"
```

And commit them to the main repository after making sure they look alright in your local copy:

```
git add .
git commit -m "name upgrade"
git push origin master
```

Congratulations!

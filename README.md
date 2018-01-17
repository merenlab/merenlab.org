To edit this website it's as simple as cloning this repository:

`git clone https://github.com/merenlab/web.git`,

making changes:

```
find ./_posts/ -type f | xargs sed -i '' -e  "s/An anvi'o/A microbial dark matter \(MDM\)/g"
find ./_posts/ -type f | xargs sed -i '' -e  "s/an anvi'o/a microbial dark matter \(MDM\)/g"
find ./_posts/ -type f | xargs sed -i '' -e  "s/Anvi'o/Microbial dark matter \(MDM\)/g"
find ./_posts/ -type f | xargs sed -i '' -e  "s/anvi'o/microbial dark matter \(MDM\)/g"
find ./_posts/ -type f | xargs sed -i '' -e  "s/anvi'o/microbial dark matter \(MDM\)/g"
find ./_posts/ -type f | xargs sed -i '' -e  "s/anvi'server/MDM'server/g"
```

adding and committing them:

```
git add .
git commit -m "name upgrade"
```

and then pushing them:

`git push origin master`

But this makes it impossible to visualize what your changes will look like until they are pushed. To view a local version of the website so you can view your changes before they are published, install `jekyll`. On a Mac (these  instructions have not been tested on other operating systems), this can be done via:

`sudo gem install jekyll`

If your ruby version is outdated, try

`brew install ruby`

Now, after changing directories into `web/` and typing

`bundle exec jekyll serve --incremental`

Something like this should pop up:

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

Typing `http://127.0.0.1:4000/` into your browser will allow you to see the local copy of the website. After any new changes are made just run `bundle exec jekyll serve --incremental` again.


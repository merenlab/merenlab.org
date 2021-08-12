# How to run a copy of this website locally

Installing ruby on Mac is always a pain. But you can start with this (newer versions of ruby may cause you trouble with some of the gems):

```
brew install ruby@2.7
```

This step may ask you to add things to your `~/.bash_profile` to update your `$PATH`, which you should do, and then either open a new terminal, or run `source ~/.bash_profile`, so the output for `which ruby` looks something like this:

```
/usr/local/opt/ruby@2.7/bin/ruby
```

Once you're done with these initial steps, you can get a copy of the web site from GitHub:

```
mkdir -p ~/github
cd ~/github/

git clone https://github.com/merenlab/web.git
```

To make sure you have the necesary gems, run:

```
cd web
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

If you are not seeing some of the changes you expect to see, press `ctrl-c` to stop the server on your termianl, clean out the static web directory by running the command `rm -rf _site/`, and re-run the server using the command above.

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

# Notes

If you want to show/hide content, you can use this notation in your markdown files:

```
<details markdown="1"><summary>Show/Hide SOME CONTENT</summary>

SOME CONTENT GOES HERE

</details>
```

If you want to show summary sections with a different background color, you can use this notation:

```
<div class="extra-info" markdown="1">

<span class="extra-info-header">Smart title for extra info</span>

EXTRA INFO GOES HERE

</div>
```

You should feel free to use warning and notice statements:

```
{:.warning}
A warning messages goes here.

{:.notice}
A notice statement goes here.
```

When naming new posts, make sure the are not **POST-DATED**!

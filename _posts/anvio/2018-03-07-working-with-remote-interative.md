---
layout: post
title: "Running remote anvi'o interactive interfaces on your local computer"
excerpt: "A hacker's guide to anvi'o networking. Are you scared? You should be!"
modified: 2018-03-07
tags: [tutorial, interactive]
categories: [anvio]
comments: true
authors: [ozcan]
---


{% include _toc.html %}


We had decided to use web browsers for our visualization needs in anvi'o from the get go. Which served us very well. Thanks to that decision, and our follow up decision to implement our visualization library from scratch, we were able to interactively work with tens of thousands of SVG objects even on a laptop computer. Kudos to all web browser developers out there who put a lot of effort into making browsers perform well, even for extreme cases of stress like the one anvi'o creates.

That said, not all steps in anvi'o workflows can be run on our laptops, so we often use server systems elsewhere to run computationally intensive or time consuming taks. The main problem this creates is that we end up having to transfer tremendous amount of data back and forth to able to work with them interactively. In addition to that, the changes we make on our local copy of the data (for example the bins or annotations we have modified) do not become immediately available for others who have access to the data on the server.

There are multiple ways to solve this. Such as setting up a [VNC server](https://en.wikipedia.org/wiki/Virtual_Network_Computing), or, as we wrote about before, using an [SSH tunnel with port forwarding]({% post_url anvio/2015-11-28-visualizing-from-a-server %}).

I have another suggestion that will likely work much seamlessly for you, and it is easy enough for anyone to setup with their server system without bothering their sys admins if they are using anvi'o in a virtual environment as they should. This is how we work with our server systems at MerenLab.

## A proposal for a better workflow ##

This part will walk you through how to set your local and remote machines for a seamless workflow. Most of the following commands are simply copy-pasteable, but please read everything carefully just to make sure you are on top of everything. This recipe will only work if you have setup your anvi'o installation within a virtual Python environment. If you haven't done that, you may need extra steps.

### The remote: the fake webbrowser module

This is the part where you will feel like a hacker, but don't worry and bear with me.

Anvi'o uses Python's `webbrowser` module when it needs a new browser window for interactive interfaces. But on a server system, this module is not very useful since there is no graphical user interface. For this setup to work, we will need to *trick* this module.

When you are using a virtual Python environment, Python first looks for modules inside the active virtual environment before looking at other locations. So first, we will create a new webbrowser module in our virtual environment **on the server side**:

``` python
cat << EOF > $(python -c "import os; print(os.path.join(os.path.dirname(os.path.abspath(os.__file__)), 'webbrowser.py'))")
def open_new(url):
    print("OPEN_ON_LOCAL[" + url + "]")

def register(*args, **kwargs):
    pass

def BackgroundBrowser(*args, **kwargs):
    pass
EOF
```

So now, when I run this command,

``` bash
python -c 'import webbrowser; print(webbrowser.__file__.rstrip("c"))'
```

I see this output,

``` bash
{path-to-anvio-virtual-env}/lib/python3.5/site-packages/webbrowser.py
```

This means every time anvi'o will request a new browser window, our fake module will just print `OPEN_ON_LOCAL[url]` to the terminal. This is the information we will catch from our local machine for magic to happen.

In order to test this, you can run the command below:

```python
python -c "import webbrowser; webbrowser.open_new('https://www.google.com/')"
```

It should print

```bash
OPEN_ON_LOCAL[https://www.google.com/]
```

If you made it this far, you're golden.


### The remote: allocating port numbers

When you run `anvi-interactive`, `anvi-display-contig-stats`, `anvi-display-pan`, or any other anvi'o interactive program, they first pick a port number to serve the interface. Unless you specified which port to use via the `--port` parameter, anvi'o will look for available ports starting from the port number `8080`.

In the [latest release](https://github.com/merenlab/anvio/releases) of anvi'o we extended this behavior. Now it first checks the `ANVIO_PORT` environment variable, and if it is set, uses the value stored in that variable instead of the default.

Using this feature, we can allocate port numbers for each user if there are multiple users, so everyone can request a different port number on the server without stepping each others' toes.

To set this up for our lab, I put the code below at the end of our `{path-to-anvio-virtual-env}/bin/activate` script on the server: 

```
[[ "$(whoami)" = "meren" ]] && export ANVIO_PORT=8080
[[ "$(whoami)" = "oesen" ]] && export ANVIO_PORT=8090
[[ "$(whoami)" = "awatson" ]] && export ANVIO_PORT=8100
[[ "$(whoami)" = "ashaiber" ]] && export ANVIO_PORT=32120
[[ "$(whoami)" = "ekiefl" ]] && export ANVIO_PORT=8120
[[ "$(whoami)" = "slee" ]] && export ANVIO_PORT=8130
```

If this is setup right, everytime you login to your server, you should be able to run this command,

``` bash
echo $ANVIO_PORT
```

and see the port number you allocated for yourself. For me, it goes like this:

``` bash
echo $ANVIO_PORT
8090
```

With this, we are done in the remote. Now we will setup the local.


### The local: creating an alias

For a proper setup in our local computers we need two things: creating an alias for the connection mode with correct port forwarding setup, and to pipe the output of the SSH connection to a script so it can *catch* the "open on local" message from the remote and actually do something with it.

First, the alias:

```
alias barhali="ssh -L 8090:localhost:8090 ozcan@barhal.mbl.edu | tee /dev/tty | python ~/.ssh/run_webbrowser.py"
```

{:.notice}
Our server's name is `barhal`, so `barhali` is an alias to open an *interactive* connection. Don't forget to change the name to your own liking, and update the connection protocol to your server, for which I used `ozcan@barhal.mbl.edu`. If you follow [this recipe](http://www.linuxproblem.org/art_9.html), you can setup your SSH for a paswordless login, in which case everything would be even more smoother.

An important note here is to setup the port numbers right for each user's local computer. I used `-L 8090:localhost:8090` in the example above becasue it is what I put in for my user in the server.

You can also repeat this parameter multiple times if you may need to run multiple interactive interfaces simultaneously:

```
alias barhali="ssh -L 8090:localhost:8090 8091:localhost:8091 8092:localhost:8092 ozcan@barhal.mbl.edu | tee /dev/tty | python3 ~/.ssh/run_webbrowser.py"
```

The last step of this command uses `python3`, please make sure that you have `python3` on your system. If you see Version 3.x when you run `python` you can also use just `python` instead `python3`.

One thing you should remember in this case is that when you open the second connection by typing `barhali` in another terminal, the SSH client will print warning messages on the screen mentioning that the port number you're trying to use is already in use. You can ignore those messages. Everything will work.

### The local: creating the final piece

The final thing is to create the script that handles the rest. You can copy the following content:

``` python
import sys
import webbrowser

for line in sys.stdin:
    if "OPEN_ON_LOCAL[" in line:
        line = line.replace("0.0.0.0", "127.0.0.1")
        webbrowser.open(line.split("OPEN_ON_LOCAL[")[1].split("]")[0])
```

and paste it into the following file path on your local machine:

```
~/.ssh/run_webbrowser.py
```

---

That's it!

With this setup we can type `barhali` on our laptops to connect to our server machine, and when we run interactive interfaces there, things would show up on our local machines automatically.

Try it, and let us know how it went for you.

If you have any questions find us on Slack, or send a comment down below.


---
layout: page
noleftpanel: true
title: "Programming for Life Scientists"
author: Course Plan
authors: [meren]
date: February, 2026
geometry: "left=3cm,right=3cm,top=2cm,bottom=2cm"
output: syllabus.pdf
colorlinks: true
urlcolor: blue
monofont: DejaVuSansMono.ttf
---

## Introduction to PFLS

The purpose of Programming for Life Scientists (PFLS) is to provide you with a general understanding of,

* The UNIX terminal environment,
* Common tools in the UNIX shell and how they are used,
* Shell scripting and its application to real-world problems,
* General concepts and algorithmic thinking in programming through exercises,
* Basics of Large Language Models and AI-assisted program solving,
* Orchestrating multiple of these tools to solve complex tasks.

This document includes **everything** you will need for this course.

## Faculty and Communication

The course is directed by [Prof. Dr. A. Murat Eren](https://merenlab.org), who goes by Meren, and [Prof. Dr. Sarahi Garcia](https://miint.org/), who goes by Sarahi. The lectures and exercises will be primarily delivered by Meren, but the following table lists everyone, including the course TAs for the current iteration of the course:

|Name|Role|Expertise|Contact information|
|:---|:---|:---|:---|
|**Meren**|Professor|Microbial Ecology, Computer Science|meren@hifmb.de|
|**Sarahi**|Professor|Microbiology, Microbial Ecology|sarahi.garcia@uni-oldenburg.de|
|**Ahmed Belfaqih**|TA|Microbiology, Bioinformatics|ahmed.belfaqih@uni-oldenburg.de|
|**Ghazaleh Sheikhi**|TA|Microbiology, Bioinformatics|ghazaleh.sheikhi.ghahi@uni-oldenburg.de|

Throughout the course (and beyond) you can reach out to Meren with questions, who should be your first contact for anything related to the activities that will follow.


## Course Responsibilities and Grading

{:.notice}
This section is relevant to you *only* if this is as a part of your coursework at the University of Oldenburg in a degree program and you wish to be graded for it.

Course responsibilities, such as the *class citizenship* emails and the details of grading is in [this document](grading).

## Course Facts and Logistics

The design of the course assumes that you know *next to nothing* about any of the topics it covers. High expectations poison everything: so please assume that you will unlikely become an expert of the tools and approaches that we will cover throughout the course, but hopefully, our discussions will provide you with enough understanding of the fundamentals of programming, helping you think about how to make computers work for you, which will give you enough *foundation* to become an expert of everything the course covers. If you are taking this course in person, we will spend about 60 hours together. That is not nothing, but it is hardly enough time to become truly skilled at anything. In his book "*[Outliers: The Story of Success](https://en.wikipedia.org/wiki/Outliers_(book))*", Malcolm Gladwell, a Canadian author and thinker, argues that a person could become an expert in nearly any field as long as they were willing to devote 10,000 hours to studying and practicing the subject or skill. Reminiscing his own journey, he says, and I am paraphrasing here from his book, "*I was completely overwhelmed at the beginning, but by the end, I felt like an expert. It took me 10 years — exactly that long*".

The previous cohort of students who took this course shared in their feedback that the pace was occasionally fast, that two weeks was not enough time, and that they sometimes felt overwhelmed. I suspect you may feel the same way. But here is what I want you to recognize: every single one of them also made it through, and most of them thought that got something important out of this. In my opinion what made their experience worthwhile at the end despite the challenges they experienced throughout the course was not that they understood everything perfectly, but that they started to ask questions when they were confused, started to feel more comfortable sitting with discomfort when they felt that they were not getting something, and they trusted the process. I will do my very best to answer every question you have, no matter how basic it seems to you. I will revisit concepts multiple times if needed, and try to explain them differently until they click. What matters most is not that you memorize every command or understand every line of code by the end of our time together. What matters is that you develop a broad appreciation of what is possible. And that you learn to recognize the shape of a problem that can be solved with a general strategy. That recognition will be your foundation, and it is all I can offer here. Absolute mastery of concepts we will discuss here is also possible, and it will come later with practice if you want that. But you cannot practice what you do not know exists, and that is why we are here. To recognize what exists. So please, ask questions. Participate in discussions. Tell me when something does not make sense. The only way to fail this course is to stay silent when you are lost.

While the course is designed for life scientists and especially microbiologists in mind, it should be useful for anyone to develop a sufficient understanding of the topics the course aims to communicate. But we will be using mock and real-world datasets and problems that convey typical characteristics of what researchers often encounter in the data-enabled era of microbiology.

The course is designed to be delivered within about two weeks (cross your fingers), and it will feel like a sprint rather than a marathon. But I anticipate that each one of you will learn some new things and enjoy your experience (cross remaining fingers?).

The plan is that the **first week** will offer insights into the terminal environment, common UNIX tools, effective use of shell, and shell scripting in general. Throughout the **second week** we will discuss version control, cloud-based solutions for collaborative coding, and AI-assisted problem solving through hands-on problem solving sessions. Throughout the course Meren will have this document open, as well as a terminal, and while following the content he will often go back to the terminal to demonstrate how things work and help you gain hands-on experience by helping you try and understand everything happening in his terminal window.

In-person attendance is **extremely important** since those active discussions are what will matter most here for your learning experience, and not this document without a narrator attached to it.

The following list offers a more detailed list of goals of the course:

* UNIX basics & file navigation
* File manipulation (`cat`, `less`, `cp`, `mv`, etc.)
* Searching files and working with data streams (`grep`, `find`, `awk`, etc.)
* Pipes and redirection (`|`, `>`, `>>`)
* Shell scripting basics (`for`, `while`, `if`, etc.)
* Basics of Git (`init`, `clone`, `commit`, `push`, etc.)
* Creating repositories & pushing scripts
* Application of shell scripting to real world, non-trivial challenges
* Using Large Language Models such as ChatGPT or DeepSeek for programming
* Helping you get hands-on experience through interactive exercises
* Helping you put your imagination and learnings in use through assignments
* Gaining insights into reproducible reporting
* And getting you to participate in lots and lots of discussions
* If we have time, we will also covering Python basics and discuss variables, data types, operators, flows, and controls

It may look like this course does not have a conventional structure and it is all over the place. You are correct, though this is intentional. I hope that this *modern structure* will work for most of you since it has multiple qualities:

* **Balanced Approach**. Rather than focusing on a single topic in great depth, the course structure goes in and out of fundamentals of programming, UNIX tools, BASH scripting, AWK programming language, Git, GitHub, and Large Language Models (LLMs) in a structured way.
* **AI Tools Introduction** – Placed at an optimal time to help you learn efficient debugging and script generation with ChatGPT, Claude, or DeepSeek before diving into difficult tasks.
* **Real-World Applications** – Exercises resemble real-world problems researchers often run into during their day-to-day workflows rather than unrelateable hypothetical programming tasks.
* **Version Control & Reproducibility** – Introduction to version control systems and cloud solutions for collaborative work ensures that you will develop good coding habits early in your journey.
* **Quick and Comprehensive Recall**. You will get to apply everything you have learned to solve an actual problem (that literally happened and someone had to solve it in their daily work).

The last time the course was delivered to a group of students with no terminal exposure, we dynamically changed the direction of the course based on what is working and what is not working. For instance, revisiting some of the earlier concepts once the participants had a working vocabulary and familiarity with scripting proved to be extremely effective for deeper understanding of the topic, so active input from the participants regarding what they are struggling with is quite important.

## Technical Setup and Recommendations

Here are a few recommendations that will help you throughout this course and beyond if you plan to apply what you have learned. You may not be familiar with some of the terms below, but we will get to them one by one.

### Shell Access

We will make quite a heavy use of the 'terminal environment' (any *terminal* that gives access to a *UNIX* [shell](https://en.wikipedia.org/wiki/Unix_shell), like [BASH](https://en.wikipedia.org/wiki/Bash_(Unix_shell))).

If you are using Linux or Mac OSX, you have native access to a reasonable shell. Please take a moment to find out how to open your terminal now.

If you are using Mac OSX, I would strongly recommend you to install [iTerm2](https://iterm2.com/) and use it instead of the default terminal application on Mac OSX. If you are using Windows, you need to install [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install) to have access to a UNIX shell.

We will spend a lot of time learning about the UNIX shell, but here are a few resources if you would like to take a brief look ahead of the course:

* [Beginner's Guide to the BASH](https://www.youtube.com/watch?v=oxuRxtrO2Ag) (a video introduction to the command line environment -- although Joe Collins is talking about Linux, the topics are relevant to anyone who uses a command line environment and Meren strongly recommends everyone to take a look).
* [Learning the Shell](http://linuxcommand.org/lc3_learning_the_shell.php) (a chapter from the open book "*The Linux Command Line*" by William Shotts -- Meren highly recommends).

### Code Editors

Throughout this course you will be writing and editing lots of code. Writing good code comfortably requires a good and comfortable editor that is *designed* to write code, which excludes text editors such as Microsoft Word or NotePad which are more appropriate for daily writing tasks. While Meren exclusively uses [vim](https://www.youtube.com/watch?v=-txKSRn0qeA) for his everyday coding tasks, students are encouraged to consider using a **graphical code editor with syntax highlighting** that supports multiple programming languages.

You can consider [Sublime Text](https://www.sublimetext.com/), which works on Mac, [Geany](https://www.geany.org/), which works on Linux, or [Notepad++](https://notepad-plus-plus.org/) which works on Windows. These are all very nice lightweight editors. Alternatively you can consider [VS Code](https://code.visualstudio.com/), which is a behemoth that works on Windows, Mac, and Linux.

{:.notice}
During our previous experience with this course, most Windows users found **VS Code** to be the most convenient way of following the course since it gives access to an editor as well as a terminal in the same window, and connects to the WSL environment which makes life very easy. Some of the remaining participants quickly learned basic commands in **vim**, and started using it with a combination of `CTRL+Z` (to go back to the terminal to test their script) and `fg` (to quickly come back to the editor to continue editing their script).

### Version Management System

The course will require you to use [Git](https://en.wikipedia.org/wiki/Git) version management system, and have an account on [GitHub](https://github.com), which is an online service built to store Git repositories in the cloud. We will discuss both git and GitHub in detail later, but you should open an account on GitHub unless you already have one as soon as you read these sentences.

## Course Data Package

Hands-on exercises throughout the course will make use of a previously prepared data package. You need to download this data package on your computer, and uncompress it. You will do it only once (unless the content changes mid-course and you are asked to do it again). Every time you open your terminal you will go into the data pack directory throughout the course since everything we will do in the terminal environment will assume that you are in the data pack directory.

Since at this point everyone has their terminals ready, all you need to do is to open a terminal, and paste the following commands in it to get a copy of the course data pack:

```sh
## change your work directory to your 'home' directory
cd ~

## download the package on your computer using 'curl'
curl -L 'https://www.dropbox.com/scl/fi/kjmkzv35fx4vv4pqo6qhl/PFLS-DATA-PACKAGE.tar.gz?rlkey=brll5t615ubsjx2bbdvjg9gme' -o PFLS-DATA-PACKAGE.tar.gz

## unpack the archive
tar -zxvf PFLS-DATA-PACKAGE.tar.gz

## go into the data pack directory (the part that you will do over
## and over again when you open a new terminal with the intention
## to follow the course content
cd PFLS-DATA-PACKAGE
```

If you got an error in any of the lines above, please do not continue before addressing the issue by working with Meren.

## Final Checks

This is time for Meren to make sure that every participant,

* Has access to a properly set up computer (with WSL for Windows users, etc),
* Has a working terminal environment that runs BASH,
* Has the data package downloaded and ready to go,
* Is able to run `git`,
* Has a GitHub account,
* Is enthusiastic about this course and puts their war paint on to deal with whatever it will bring into their life.

If we are all good to go, let's start now with the course content.

## Course Content

### The Shell

The purpose of this section is to familiarize you with your terminal and some of the features of the shell environment.

At the end of this section, you will have an understanding of the power of the command line environment and how to use multiple programs in tandem to get things done interactively.

#### What is shell, exactly?

This is a good place to start because the shell is one of the most powerful yet under-appreciated programs on your computer. In simplest terms, the shell is a command-line interface that sits between you and the operating system kernel —the hidden core of your system responsible for running program processes, executing their instructions on the CPU, managing memory, and handling storage and peripherals. At any given time, your computer runs a gazillion programs; from those that respond to your inputs through a keyboard and mouse to those that literally show where your cursor is on your screen or keeps track of all the windows you have opened, and all of them go through the kernel which allows your computer to run more programs than the number of available CPU cores, gives all programs the illusion of unlimited memory despite its physical limits, and takes care of priorities among them so all programs big and small can use the hardware resources on your computer in harmony to serve you. The shell is another program running on your computer. Like many others, it operates in what we call 'user space' and communicates with core processes running in 'kernel space' through system calls, ensuring that our commands and requests integrate smoothly into the operating system's workings without disrupting its stability. In many ways, the shell is your command center that gives you direct control over your system's most powerful functions; in a way, it is the computer equivalent of an airplane cockpit or a nuclear reactor control room. On top of its abilities to run programs, modern shells come with scripting languages and control structures for us to run things even more efficiently.

What happens when you type something in your terminal is a fascinating story, and has many many layers. Once you press Enter, the terminal sends your input to the shell, which parses the command, expands special characters and variables, and determines whether it is a built-in command or an external program. If it's an external program, the shell searches the relevant directories to locate it, and if the program is found and has the right permissions, the shell creates a new process by invoking the kernel, and puts the target program into the driver seat, upon which the kernel loads the program's code into memory, registers it with the scheduler, and begins execution. Throughout its runtime the program interacts with the kernel via system calls until it either terminates normally (successful or not) or the kernel forcefully kills it due to a myriad of reasons such as unhandled fatal signals, resource exhaustion, or illegal memory access. Regardless of how it ends, the kernel cleans up the kitchen and notifies the shell, which retrieves the exit status and resumes its eternal wait for your next command, ready to start the cycle all over again.

Shells are complex, and you don't need to learn any of these things, obviously. But as an undergraduate student of computer science I *did* implement a shell from scratch. My shell was neither as good nor as talented as any of the modern shells. But doing that, and forcing myself to go through that suffering created many new synapses in my brain that eventually afforded me the framework I needed to understand and solve more complex problems elsewhere. Which, in a big plot twist, brings me to the use of AI.

One may argue that almost none of the things we will discuss throughout this course are necessary to learn. Indeed, you can solve the vast majority of programming challenges, including those that we will cover during our exercises, using popular LLM clients such as ChatGPT, Claude, Gemini, Grok, DeepSeek, or many others almost instantaneously. But there are no cheat codes to life -- without going through the pain of truly understanding something, there is no way to achieve mastery, lead with confidence, or create something original.

<blockquote>
There is a tide in the affairs of people,<br/>
Which, taken at the flood, leads on to fortune;<br/>
Omitted, all the voyage of their life<br/>
Is bound in shallows and in miseries.

<div class="blockquote-author">William Shakespeare (Julius Caesar, Act 4, Scene 3)</div>
</blockquote>

I am asking you to not use any AI tools until I ask you to do so. Do not prevent yourself from taking advantage of the tide that is rising in your academic life due to this course. Once you have some fundamental understanding of the topics we cover here, you will be much more efficient using AI-assistance for your dealings with shell and any other coding tasks.

Now, please open your terminals and let's start with the simplest shell functions and most commonly used UNIX programs in it that help us navigate through files and directories, study their contents, and special characters that make life easier.

{:.notice}
As we go through these examples Meren will interactively demonstrate their utility in various ways, and he encourages you to also try them by literally typing them in your terminal since this activity will already help you develop a level of familiarity with the command line interface.

#### Navigating files and directories

`pwd` – Print working directory

Displays the absolute path of the current directory.

```sh
pwd  # Show full path of the current directory
```

---

`ls` – List directory contents

Lists files and directories in the current directory.

```sh
ls        # List files in the current directory
ls -l     # Long listing format (permissions, owner, size, date)
ls -a     # Show hidden files
ls -lh    # Human-readable file sizes
ls -R     # Show files recursively
```

---

`cd` – Change directory

Moves between directories.

```sh
cd /home/meren/Documents  # Change to Documents folder
cd ..                     # Go one level up
cd /                      # Go to root directory

cd -                      # Go back to the directory you were in
                          # right before this current one

cd                        # And here is a little surprise -- if
                          # you don't provide any parameters,
                          # cd will take you to your home
                          # directory
```

---

`mkdir` – Create a new directory

Creates a new empty directory.

```sh
mkdir my_folder       # Create a directory named my_folder
mkdir -p parent/child # Create nested directories
```

---

`touch` – Create an empty file

Creates a new empty file.

```sh
touch newfile.txt  # Create an empty file called newfile.txt
```

---

`rm` – Remove files or directories

Deletes files or directories.

```sh
rm file.txt         # Delete a file
rm -r my_folder     # Delete a directory and its contents
rm -rf my_folder    # Force delete (dangerous, but fun)
```

---

`cp` – Copy files and directories

Copies files and folders.

```sh
cp file1.txt file2.txt         # Copy file1.txt to file2.txt
cp -r dir1 dir2                # Copy directory dir1 to dir2
cp file1.txt file2.txt backup/ # Copy both files into the backup folder
```

---

`mv` – Move or rename files

Moves or renames files and directories.

```sh
mv oldname.txt newname.txt   # Rename a file
mv file.txt my_folder/       # Move a file into my_folder
mv folder1 folder2           # Rename folder1 to folder2
```

---

`basename` - Strip directory and suffix from filenames

Extremely useful when you have to deal with full paths for files.

```sh
## this will return 'file.txt'
basename path/to/a/file.txt

## this will return 'file'
basename path/to/a/file.txt .txt
```

---

`echo` - Display a line of text, or print the value of a variable

This is not quite related to files and directories, and more of a very important special function in shell, but it is important that we appreciate what it is before the next section.

```sh
echo "hello everyone!"
```

`echo` simply prints things out for us. I know it looks dumb, but you will soon see its utility.

---

`find` – Search for files and directories

Finds files in a directory hierarchy.

```sh
find /home -name "document.txt"  # Search for document.txt in /home
find /var -size +100M            # Find files larger than 100MB in /var
find . -type f -name "*.log"     # Find all .log files in current directory
```

The last command uses `*`, a special character (so-called 'wildcard') that is extremely useful to target multiple files that match to a particular pattern. This is not the only special character, and most shells will process user input commands and interpret a series of special characters when they are found. Before we continue with more fun programs, let's take a look at a list of commonly used special characters first.


#### Special Characters

`$` – Variable substitution and command substitution

Used to reference variables and execute commands. Extremely important character that we will use many many times.

```sh
echo $HOME       # Prints the home directory
echo $(date)     # Runs the date command and prints the result
```

---

`#` – Comment

Everything after `#` on a line is ignored by BASH. It is mostly useful when writing shell scripts as you will see soon.

```sh
## This is a comment
echo "Hello World"  # Prints Hello World
```

---

`*` – Wildcard (matches multiple characters)

Matches all files and directories in a given location. A life saver and your biggest enemy.

```sh
ls *.txt    # Lists all files ending in .txt
```

---

`?` – Single character wildcard

Matches any single character.

```sh
ls file??.txt  # Matches file01.txt, file02.txt, etc, but not fileX.txt
```

---

`.` – Current directory

Refers to the current directory.

```sh
ls .          # List files in the current directory
./script.sh   # Execute script.sh in the current directory
```

---

`..` – Parent directory

Refers to the directory one level above the current directory.

```sh
cd ..   # Move to the parent directory
```

---

`~` – Home directory

Represents the current user's home directory.

```sh
cd ~           # Go to the home directory
ls ~/Downloads # List files in the "Downloads" directory inside home
```

---

`>` – Redirect output (overwrite)

Sends output to a file (overwrites if file exists).

```sh
echo "Hello" > file.txt  # Writes "Hello" into file.txt (overwrites)
```

---

`>>` – Append output to a file

Appends output to the end of a file instead of overwriting.

```sh
echo "Hello" >> file.txt  # Writes "Hello" into file.txt (appends)
```

---

`<` – Input redirection

Reads input from a file.

```sh
sort < unsorted.txt   # Reads unsorted.txt as input for the sort command
```
---

`&&` – Logical AND (run second command only if first succeeds)

Runs the second command **only if the first one succeeds**.

```sh
mkdir new_dir && cd new_dir  # Creates and moves into `new_dir` *if* successful
```

---

`||` – Logical OR (run second command if first fails)

Runs the second command **only if the first one fails**.

```sh
mkdir my_dir || echo "Directory creation failed"  # Prints message *if* mkdir fails
```

---

`;` – Command separator

Allows multiple commands on the same line.

```sh
echo "Hello"; echo "World"  # Prints Hello then World
```

---

`"` (Double Quotes) – Allow variable expansion

Expands variables inside.

```sh
echo "Home: $HOME"  # Prints: Home: /home/user
```
---

`'` (Single Quotes) – Preserve literal text

Prevents expansion of variables and special characters.

```sh
echo '$HOME'  # Prints: $HOME (without expanding)
```

---

`\` – Escape special characters

Prevents special interpretation of characters.

```sh
echo "This is a quote: \"Hello\" :)"  # Prints: This is a quote: "Hello"
```

---

`!$` – The last argument from the last command

Repeats arguments from the previous command.

```sh
ls -l /etc
echo !$  # Expands to: echo /etc
```

---

`>` and `2>` – Redirect standard and error output

Redirects normal output (`>`) and error output (`2>`).

```sh
ls > output.txt    # Saves normal output
ls nonexistent 2> error.txt  # Saves error output
ls nonexist &> all.txt  # Redirects both normal and error output
```

---

`|` – Pipe (send output of one command to another)

This one implements one of the most important and powerful concept in the shell environment: it passes the output of one command as input to another. By doing so, it enables the most effective use of the shell environment through a dynamic orchestration of multiple tools as you will see more and more throughout the rest of the course.

```sh
ls | grep "file"   # List files and filter for those containing "file"
ps aux | grep ssh  # Show running processes related to SSH
```

#### Working with files

`cat` – Display file contents

Concatenates and displays file content.

```sh
cat file.txt         # Show full file content
cat file1.txt file2.txt > merged.txt  # Merge two files into one
```

---

`tac` – Display file in reverse order

Prints a file **from bottom to top**.

```sh
tac file.txt  # Show file.txt in reverse order
```

---

`less` – View file contents page by page

Allows **scrolling** through large files.

```sh
less largefile.txt   # Open a large file for viewing
```

You can press `q` anytime to quit and go back to your terminal, or type `/search_term` to find a term.

---

`head` – Show the first few lines of a file

Displays the first 10 lines by default.

```sh
head file.txt       # Show the first 10 lines
head -n 5 file.txt  # Show the first 5 lines
```

---

`tail` – Show the last few lines of a file

Displays the last 10 lines by default.

```sh
tail file.txt        # Show the last 10 lines
tail -n 5 file.txt   # Show the last 5 lines
tail -f log.txt      # Continuously show new lines (useful for logs)
```

---

`grep` – Search for text in a file

Finds lines containing a specific pattern.

```sh
grep "error" logfile.txt       # Find lines containing "error"
grep -i "warning" logfile.txt  # Case-insensitive search
grep -r "function" /home/code  # Search in all files under /home/code
```

---

`sed` – Stream editor (modify file content)

Used for **text replacement** and processing.

```sh
sed 's/old/new/g' file.txt   # Replace "old" with "new" in file.txt
sed -i 's/foo/bar/g' file.txt # Edit file in place
```

---

`sort` – Sort lines in a file

Sorts text **alphabetically** or **numerically**.

```sh
sort names.txt         # Sort lines alphabetically
sort -r names.txt      # Reverse order sort
sort -n numbers.txt    # Sort numerically
sort -u names.txt      # Remove duplicate lines
```

---

`uniq` – Remove duplicate lines from sorted text

Filters out repeated lines in a file.

```sh
sort names.txt | uniq    # Remove duplicates after sorting
uniq -c names.txt        # Show duplicate counts
```

---

`wc` – Count words, lines, and characters in a file

Displays **line count, word count, and byte size**.

```sh
wc -l file.txt  # Show the number of lines
wc -w file.txt  # Show the number of words
wc -c file.txt  # Show the number of characters
```

---

`diff` – Compare two files line by line

Finds differences between two files.

```sh
diff file1.txt file2.txt    # Show line-by-line differences
diff -y file1.txt file2.txt # Show side-by-side comparison
```

---

`awk` – Pattern scanning & processing

AWK is a data-driven language that excels in text processing, data extraction, and reporting, and it is one of the most amazing little tools you will find in the UNIX shell environment (fun fact, AWK is not a meaningful acronym since it simply comes from the names of the developers: Alfred **A**ho, Peter **W**einberger, and Brian **K**ernighan :)).

In a nutshell, AWK scans a file line by line, splits input into fields based on a separator, enables pattern-based filtering, and allows users to perform actions on matching lines, and helps produce highly formatted reports.

It is really difficult to demonstrate the utility of AWK without a few examples, so I put together the following file of all chancellors of Germany where the columns indicate (1) the name of the chancellor, (2) their education, (3) the age at which they became a chancellor of Germany, (4) the number of years they served at this position, and (5) the year they assumed this position. You can save the contents of this file as `german_chancellors.txt` in your working directory, and follow the examples below:

```
Konrad_Adenauer   Law          73   14   1876
Ludwig_Erhard     Economics    65   3    1897
Kurt_Georg_Kiesinger  Law      62   3    1904
Willy_Brandt      History      55   5    1913
Walter_Scheel     Law          54   1    1919
Helmut_Schmidt    Economics    56   8    1918
Helmut_Kohl       History      52   16   1930
Gerhard_Schröder  Law          54   7    1944
Angela_Merkel     Physics      51   16   1954
Olaf_Scholz       Law          63   Ongoing  1958
```

Print the name, age when becoming chancellor, and birth year:

```sh
awk '{print $1, "became chancellor at age", $3, "and was born in", $5}' german_chancellors.txt
```

Find chancellors born before 1920:

```sh
awk '$5 < 1920 {print $1, "was born in", $5}' german_chancellors.txt
```

Calculate each chancellor's age today (assuming we are still in 2026 by the time you're seeing this):

```sh
awk '{print $1, "is (or would have been)", 2026 - $5, "years old today"}' german_chancellors.txt
```

Find chancellors born after women gained voting rights in Germany:

```sh
awk '$5 > 1918 {print $1, "was born in", $5}' german_chancellors.txt
```

---

Let's say we want to print the first letter of the first name and last name of each chancellor who studied law given this file. How can we do that?

First we can filter this file for chancellors who studied law. There are multiple ways to do it. For instance, this is one way to do it:

```sh
awk '{if ($2=="Law") print}' german_chancellors.txt
```

but a more stylish way to do it in AWK would make use of patterns:

```sh
awk '/Law/{print}' german_chancellors.txt
```

where `/Law/` is a pattern and means "run the following statement on every line that contains the word Law".

To get the first letter of the first name and last name of each chancellor, we first need to get their names. So we can do it by printing the exact column where names appear:

```sh
awk '/Law/{print $1}' german_chancellors.txt
```

We see where the first and the last names are, but how to access them separately? Well, one way to do it is to 'split' this text into fields. A special keyword in AWK is `FS`, which means field separator. By default, the field separator is whitespace, but we can replace it with anything BEFORE running AWK commands on anything by simply mentioning it to AWK at the very beginning of the process.

Here is an example:

```sh
awk '/Law/{print $1}' german_chancellors.txt | awk 'BEGIN{FS="_"} {print $1 ", " $2}'
```

Cool. Now we can access name and lastname as separate columns. But how to get the first letter of each of these? Well, if you were to search on Google "how to get the first letter of a variable in AWK", you would learn that there is an AWK 'function' called `substr`, and we can use that to get what we want:

```sh
awk '/Law/{print $1}' german_chancellors.txt | awk 'BEGIN{FS="_"} {print(substr($1, 1, 1) substr($2,1,1))}'
```

---

OK. Here is a relatively difficult one since it requires some level of algorithmic thinking: find the average age at which chancellors took office. This will take care of it, and let's break it down to its individual components to discuss what we are looking at here:

```sh
awk '{sum_age += $3; count++} END { print "Average age when becoming chancellor:", sum_age/count}' german_chancellors.txt
```

### EXC-001

Let's do an exercise using a few of the things we have learned so far. For this exercise, please go into the relevant exercise directory in the course data pack:

```sh
cd EXC-001
```

Where you will find a FASTA file. The FASTA file was generated by an anvi'o user, and it describes all the genes that are found in a single bacterial genome that resolves to *Vibrio jasicida*. Please try to answer the following questions:

* Are you able to take a look at the file and see its general structure?
* How many genes are there in this FASTA file?
* How many genes are there from the contig `c_000000000023`?
* What is the name of the longest gene in the file?
* What is the total number of genes in the top ten contigs with most genes?
* Each gene sequence in this FASTA file can spread across multiple lines -- can you create a new FASTA file called `genes-one-line.fa` in which every gene sequence occupies a single line?

Try your best, and it is OK if you can't answer each one of them. If you try your hardest, the solutions we will go through and explain together will make much more sense even if you fail.

{:.warning}
**Please turn in your solutions the following way**: Copy-paste the questions above into your email client, under each question write your final answers along with the command line that led to that answer, and send the email to _meren@hifmb.de_ **and** _ahmed.belfaqih@uni-oldenburg.de_. The subject line of your email must be `PFLS EXC-001` :)

Once we are done, we will review the [solutions](solutions/EXC-001) together.

## Shell Scripting

The last section focused on how the command line environment and the common tools that are accessible to us in the UNIX shell can empower its users to perform tasks that would have taken much longer to do manually. For instance, one could find how many genes are in a given FASTA file by literally going through it line by line in their text editor, but the ability to perform this task with a single command is a life saver. Using individual lines of instructions interactively is very powerful, but not suitable to complete repetitive tasks or implement complex ideas that require multiple instructions to be run one after another.

Let's talk more about what shell scripting is and where it comes in handy through a realistic example. Consider the following.

You are interested in understanding the functional landscape of the *Wolbachia* genus. To do this, you downloaded all the representative genomes from the Genome Taxonomy Database (GTDB). Since there are no cultures of *Wolbachia*, most of the genomes will be reconstructed from metagenomes or will be single amplified, both of which can yield highly fragmented genomes. But for an appropriate analysis of functions, you need genomes that are relatively well put together. Ideally a single contig, but if not, let's say no more than 30 contigs. So after downloading these genomes, you will want to count the number of contigs in each one of them, and put aside the ones that you actually would like to use for your downstream analyses.

The genomes you have downloaded looks like this:

```
GCA_018224395.1_genomic.fna
GCA_019061405.1_genomic.fna
GCA_022836975.1_genomic.fna
GCA_023052945.1_genomic.fna
GCA_902636535.1_genomic.fna
GCF_000008385.1_genomic.fna
GCF_000306885.1_genomic.fna
GCF_000376585.1_genomic.fna
(...)
```

{:.notice}
These genomes are also available to you in the `EXC-002` directory of your data pack. Please go into that directory now, and confirm that you can see the files in it. Being there will become handy soon when you try to test some of the steps below.

You could indeed indeed use `grep` and `wc` to put together a command that gives you the number of contigs in a given genome. Here is an example that shows how that looks like on my terminal for one of those genomes:

```sh
meren $ grep '>' GCA_022836975.1_genomic.fna | wc -l
192
```

Let's stop here and try to actually generate a report from these FASTA files together. I want you to create a table that shows each FASTA file name and the number of sequences in it, so we can determine which genomes contain fewer than 30 contigs. Go ahead and start working on creating a table for me since you have everything you need.

---

This strategy requires you to write each of these command one by one, then read the number in the output, and perhaps put it in a table in an environment like EXCEL just so you can sort it based on the number of contigs to then determine the names of genomes that contain less than a certain number of contigs, then use those names to put the genomes of interest in a different directory.

This is exactly where shell scripts come handy: repetitive tasks where the output of the task (such as the number of contigs) is not the final answer but actually the input of another independent step (such as determining whether to use that genome further).

Let's first talk about some of the fundamentals.

### What is a shell script?

Generally speaking, a script is a text file that contains instructions to be executed by an 'interpreter'. Scripts are often used to automate tasks that could have been done manually, and in that sense they are often much less complex and much more readable than programs that implement full-fledged applications. So a shell script is such a text file that is interpreted by a shell (such as BASH).

Here we can create the simplest BASH script by instead putting all the `wc` and `grep` commands we run in the last section into a single file rather than running them one by one. For that, we can create a text file called 'get-num-contigs.sh' in the same directory with our genomes and save the following lines in it,

```sh
grep '>' GCA_018224395.1_genomic.fna | wc -l
grep '>' GCA_019061405.1_genomic.fna | wc -l
grep '>' GCA_022836975.1_genomic.fna | wc -l
grep '>' GCA_023052945.1_genomic.fna | wc -l
grep '>' GCA_902636535.1_genomic.fna | wc -l
grep '>' GCF_000008385.1_genomic.fna | wc -l
grep '>' GCF_000306885.1_genomic.fna | wc -l
grep '>' GCF_000376585.1_genomic.fna | wc -l
```
{% include CODEBLOCKFILENAME filename="get-num-contigs.sh" %}

Then we can go back to the terminal window, and run it:

```sh
bash get-num-contigs.sh
```

Which would give us the following output:

```
1
30
192
1
108
1
1
1
```

While this is indeed a BASH script, it is not really a good one. Programming and scripting languages include some very common ideas to realize complex and dynamic operations through variables, loops, and means to branch into different tasks through conditionals. By doing so, they help us avoid redundancy (such as having to type `grep` as many times as there are genomes), and scale up our ideas to very large datasets. Let's talk about those concepts a bit, and then come back to this.


### Variables

Variables are placeholders that store data such as text, numbers, or even the outputs of commands. They enable scripts to be dynamic, flexible, and reusable by allowing values to be assigned, modified, and referenced throughout execution.

Defining a variable is simple. You can define one the following way:

```sh
my_variable="Hello!"
```

To access the content of a variable, we use the `$` character in front of it:

```sh
echo $my_variable
```

We can talk about BASH variables in three main classes. One of those classes, is **user-defined variables**, where you define the variable, and it has no value before you do it.

Here are some examples:

```sh
## define a variable literally with a text you just wrote
my_name="Meren"
echo "My name is $my_name"

echo "##############################"

## capture the output of a command, and put it into a variable:
my_bd="2002-05-14 12:00:00"
seconds_to_my_bd=$(date -d "$my_bd" +%s)
seconds_to_now=$(date +%s)
seconds_since_my_bd=$((seconds_to_now - seconds_to_my_bd))
echo "It has been about $seconds_since_my_bd seconds since I was born!"

echo "##############################"

## you can also define a variable that stores lists of things
fasta_files=*.fna
echo $fasta_files

echo "##############################"

## the placeholder nature of variables can lead to creative
## applications! what do you think will happen when we run
## the next two lines here?
x=ls
$x
```
{% include CODEBLOCKFILENAME filename="variables-user-defined.sh" %}

The second class of variables are **environmental variables**: variables that were previously defined by various processes, including those that are set every time you open a terminal. You can see all of these variables by simply typing the command `env` in your terminal, and you can access any of these variables from within your BASH scripts:

```sh
echo "My username is '$USER'. There are $(ls $HOME | wc -l) files in my home folder, which is at '$HOME'."
```
{% include CODEBLOCKFILENAME filename="variables-environmental.sh" %}

The third class of variables are the **built-in variables** that are set every time a command or script is run. If you have the following shell script,

```sh
echo "Script Name ..................: $0"
echo "Number of Arguments ..........: $#"
echo "First Argument ...............: $1"
echo "Second Argument ..............: $2"
echo "Fifth Argument ...............: $5"
echo "All Arguments ................: $@"

xxx &> /dev/null

echo "Exit Status of command 'xx' ..: $?"

ls &> /dev/null

echo "Exit Status of command 'ls' ..: $?"
```
{% include CODEBLOCKFILENAME filename="variables-built-in.sh" %}

Running it in your terminal with the following arguments,

```sh
bash variables-built-in.sh apple orange banana strawberry
```

Will result in the following output:

```
Script Name ..................: variables-built-in.sh
Number of Arguments ..........: 4
First Argument ...............: apple
Second Argument ..............: orange
Fifth Argument ...............:
All Arguments ................: apple orange banana strawberry
Exit Status of command 'xx' ..: 127
Exit Status of command 'ls' ..: 0
```

### Loops

Loops make everything much more fun and dynamic. One of the most famous loop forming strategies in programming is to use the `for` loop, where pretty much every language will have a version of it. A `for` loop iterates over a list of values, and allows its user to define various commands to be run on each item separately. For instance, in shell, a `for` loop has the following general structure:

```sh
for variable in list
do
    # in this block, do something
    # with the variable
done
```

`for` loops shine whenever there is a task that must be run the same way on many distinct items, such as processing multiple files the same way *wink* *wink*. The most important part of building a good `for` loop requires one to understand the `list` in the structure above. Here are a few examples to think about:

```sh
for number in 1 2 3 4 5
do
    echo $number
done
```
{% include CODEBLOCKFILENAME filename="for-loop-examples.sh" %}

In this instance, the `list` consists of letters we manually typed. We could of course also define them as a variable, and pass that variable to the for loop as the `list` to iterate over:

```sh
numbers="1 2 3 4 5"
for number in $numbers
do
    echo $number
done
```
{% include CODEBLOCKFILENAME filename="for-loop-examples.sh" %}

If I were to show you how this logic would have been implemented in different languages, you probably would recognize the structural similarities. For instance, here is the same for loop in Python:

```python
numbers = [1, 2, 3, 4, 5]
for number in numbers:
    print(number)
```

or Swift:

```swift
let numbers = [1, 2, 3, 4, 5]
for number in numbers {
    print(number)
}
```

In R:

```r
numbers <- c(1, 2, 3, 4, 5)
for (number in numbers) {
    print(number)
}
```

or in Julia:

```julia
numbers = [1, 2, 3, 4, 5]
for number in numbers
    println(number)
end
```

OK. Going back to our example in shell, the items we wish to iterate over may come from another command. For instance, the following notation in BASH will give you numbers from 1 to 5 as a sequence (you can copy paste it in your terminal to see):

```sh
echo {1..5}
```

And our numbers could come from such a function that produce the same numbers:

```sh
for number in $(echo {1..5})
do
    echo $number
done
```
{% include CODEBLOCKFILENAME filename="for-loop-examples.sh" %}

Assuming you are in the `EXC-002` directory, What do you think this will do?

```sh
for f in $(ls *.fna)
do
    echo $f
done
```

Great. Do you remember our glorious first shell script we wrote to count the number of sequences in all FASTA files in that directory?

```sh
grep '>' GCA_018224395.1_genomic.fna | wc -l
grep '>' GCA_019061405.1_genomic.fna | wc -l
grep '>' GCA_022836975.1_genomic.fna | wc -l
grep '>' GCA_023052945.1_genomic.fna | wc -l
grep '>' GCA_902636535.1_genomic.fna | wc -l
grep '>' GCF_000008385.1_genomic.fna | wc -l
grep '>' GCF_000306885.1_genomic.fna | wc -l
grep '>' GCF_000376585.1_genomic.fna | wc -l
```
{% include CODEBLOCKFILENAME filename="get-num-contigs.sh" %}

I want you to help me write a version of this that produces this exact output:

```
GCF_003344345.1_genomic.fna 237
GCA_022836975.1_genomic.fna 192
GCF_023661085.1_genomic.fna 186
GCF_014534705.1_genomic.fna 182
GCA_902636535.1_genomic.fna 108
GCF_020278625.1_genomic.fna 93
GCF_013366805.1_genomic.fna 41
GCA_019061405.1_genomic.fna 30
GCF_020405475.1_genomic.fna 12
GCF_001752665.1_genomic.fna 12
GCF_936270435.1_genomic.fna 1
GCF_936270145.1_genomic.fna 1
GCF_918342435.1_genomic.fna 1
GCF_025021925.1_genomic.fna 1
(...)
```

OK. There is one more topic to cover, and after that we will come back to these genomes, which will serve as our next exercise!

### Conditionals

So far we discussed how to define and make use of variables, and how to build loops using `for` in our shell scripts. Our `for` loops run on all items, without having to make any decisions. But real-world tasks often require some sort of decision making and performing an operation only if a certain condition is, or a few of them are, met (or do other things if they don't!). Those conditions could include a variety of considerations such as doing something only if a file exists or is absent, comparing variables and taking action depending on whether they are equal or not, and thus controlling the flow of our script based on our expectations from it.

Shell scripts can use `if` statements to check for conditions. The most general structure of an `if` statement is the following:

```sh
if [ condition ]
then
    # in this block, do something
    # with whatever
fi
```

For instance, let's reconsider our simple example with the numbers:

```sh
for number in 1 2 3 4 5
do
    echo $number
done
```
{% include CODEBLOCKFILENAME filename="for-loop-examples.sh" %}

If we wanted to print out only the numbers larger than 2, our 'condition' statement in the general structure above would require a test that determines whether a given number is greater than 2. This particular example will satisfy that:

```sh
for number in 1 2 3 4 5
do
    if [ $number -gt 2 ]
    then
        echo $number
    fi
done
```
{% include CODEBLOCKFILENAME filename="if-else-examples.sh" %}

Here the *condition* is `$number -gt 2`, and the `-gt` operator, which means '*greater than*', does the heavy lifting of the entire operation. When a number is truly greater than 2, that *condition* becomes a true statement, and it meets the criterion for the `if` statement to continue.

If we wished to consider only the numbers that are *not* greater than 2, we could simply change this condition to become a true statement only when the `$number -gt 2` is false by UNO reversing the entire thing with an exclamation mark:

```sh
for number in 1 2 3 4 5
do
    if [ ! $number -gt 2 ]
    then
        echo $number
    fi
done
```
{% include CODEBLOCKFILENAME filename="if-else-examples.sh" %}

What if we didn't want to use `!`? What is it we would need to do to get the numbers that are *not* greater than 2?

This is a good time to have a look at all the operators that one can use when they are forming their conditions:

| Operator          | Description                                                                 | Example                          |
|-------------------|-----------------------------------------------------------------------------|----------------------------------|
| `-eq`             | Equal to (numeric comparison)                                               | `if [ "$a" -eq "$b" ]`           |
| `-ne`             | Not equal to (numeric comparison)                                           | `if [ "$a" -ne "$b" ]`           |
| `-gt`             | Greater than (numeric comparison)                                           | `if [ "$a" -gt "$b" ]`           |
| `-ge`             | Greater than or equal to (numeric comparison)                               | `if [ "$a" -ge "$b" ]`           |
| `-lt`             | Less than (numeric comparison)                                              | `if [ "$a" -lt "$b" ]`           |
| `-le`             | Less than or equal to (numeric comparison)                                  | `if [ "$a" -le "$b" ]`           |
| `=`               | Equal to (string comparison)                                                | `if [ "$a" = "$b" ]`             |
| `==`              | Equal to (string comparison, synonymous with `=`)                           | `if [ "$a" == "$b" ]`            |
| `!=`              | Not equal to (string comparison)                                            | `if [ "$a" != "$b" ]`            |
| `-z`              | String is null (has zero length)                                            | `if [ -z "$a" ]`                 |
| `-n`              | String is not null (has non-zero length)                                    | `if [ -n "$a" ]`                 |
| `-e`              | File exists                                                                 | `if [ -e "file.txt" ]`           |
| `-f`              | File exists and is a regular file                                           | `if [ -f "file.txt" ]`           |
| `-d`              | File exists and is a directory                                              | `if [ -d "dir" ]`                |
| `-r`              | File exists and is readable                                                 | `if [ -r "file.txt" ]`           |
| `-w`              | File exists and is writable                                                 | `if [ -w "file.txt" ]`           |
| `-x`              | File exists and is executable                                               | `if [ -x "file.txt" ]`           |
| `-s`              | File exists and has a size greater than zero                                | `if [ -s "file.txt" ]`           |
| `-h` or `-L`      | File exists and is a symbolic link                                          | `if [ -h "link" ]`               |
| `-O`              | File exists and is owned by the current user                                | `if [ -O "file.txt" ]`           |
| `-G`              | File exists and is owned by the current user's group                        | `if [ -G "file.txt" ]`           |
| `-N`              | File exists and has been modified since it was last read                    | `if [ -N "file.txt" ]`           |
| `-a`              | Logical AND (deprecated, use `&&` instead)                                  | `if [ "$a" -eq 1 -a "$b" -eq 2 ]`|
| `-o`              | Logical OR (deprecated, use `\|\|` instead)                                 | `if [ "$a" -eq 1 -o "$b" -eq 2 ]`|
| `!`               | Logical NOT                                                                 | `if [ ! "$a" -eq 1 ]`            |
| `&&`              | Logical AND (used between commands or conditions)                           | `if [ "$a" -eq 1 ] && [ "$b" -eq 2 ]` |
| `\|\|`            | Logical OR (used between commands or conditions)                            | `if [ "$a" -eq 1 ] \|\| [ "$b" -eq 2 ]` |


One can also extend `if` statements with additional conditions and `else` statements:

```sh
if [ condition ]; then
    # do something when 'condition_x' is true
elif [ condition_y ]; then
    # do something when 'condition_x' is false,
    # but condition_y is true
else
    # do something for cases that renders neither
    # 'condition_x' or 'condition_y' true
fi
```

So here is our example with numbers that makes the best use of everything `if/elif/else` statements have to offer:

```sh
for number in 1 2 3 4 5
do
    if [ $number -gt 3 ]; then
        echo "$number is greater than 3"
    elif [ $number -lt 3 ]; then
        echo "$number is less than 3"
    else
        echo "$number is neither greater nor less than 3 (program is confused)"
    fi
done
```
{% include CODEBLOCKFILENAME filename="if-else-examples.sh" %}

### EXC-002

You now have the theoretical knowledge required to write a shell script to solve the following problem and try to put your learnings in good use.

For this exercise, please go into the relevant exercise directory in the course data pack:

```sh
cd EXC-002
```

You have 34 *Wolbachia* genomes downloaded from the GTDB, and a text file that shows which genome matches to which host organism *Wolbachia* infects. Please write a shell script that does the following tasks when it is run in the `EXC-002` directory:

* Creates a directory called `WOLBACHIA-GENOMES`,
* Identifies *Wolbachia* genomes with less than `x` number of contigs, where `x` is sent to the shell script as a parameter,
* Creates a copy of each genome that matches to the above condition and puts it inside the `WOLBACHIA-GENOMES` directory,
* But uses the `wolbachia-hosts.txt` to rename each genome file to match the host name from which the *Wolbachia* was recovered.

{:.warning}
**Please turn in your solutions the following way**: Save your script as `process-wolbachia-genomes.sh`, add it as an attachment to an email with the subject line `PFLS EXC-002` and send your email to _meren@hifmb.de_ **and** _ahmed.belfaqih@uni-oldenburg.de_.

We will go through the [solution](solutions/EXC-002) together once you have given this exercise your best shot.


## Working with Git

### Introduction to Version Control

Probably everyone who is going through this document is familiar with the fact that scientific writing is not a linear process with a clear endpoint. If you have ever had multiple versions of the same document on your disk named like `final.doc`, `final_final.doc`, `really_final.doc`, `final_Jan292025,doc`, `final_Jan292025_final.doc`, you already know what I am talking about with firsthand experience. The reason we often have multiple versions of the same document like that (instead of changing the same file) is simple: we want to make sure that if we ever need to go back in time to recall a version of the text we are working on, those earlier versions of the text will not be overwritten with the newer versions of it and be forever lost. But of course this leads to a lot of confusion and redundancy.

The same happens when we write code, where the number of changes that occur in a single file can be very high. For instance, here are the top ten files in the anvi'o code based on how many times they were changed during the past 10 years.

|**Anvi'o file**|**Number of times it was changed**|
|:--|:--|
|kegg.py|1334|
|dbops.py|906|
|init.py|665|
|utils.py|575|
|interactive.py|486|
|reactionnetwork.py|460|
|bottleroutes.py|452|
|profiler.py|407|
|variabilityops.py|318|
|summarizer.py|307|
|contigops.py|304|
|trnaseq.py|299|
|sequence.py|240|
|merger.py|192|
|structureops.py|185|

It is crucial for the developers of anvi'o to be able to keep track of all these changes: if a particular change introduces a bug in the code that is discovered much later, there must be a way for anvi'o developers to go back to the version that was working well. But anvi'o developers cannot resort to creating multiple copies of every file for each major change.

This is precisely what version control systems do for you. They offer you means to record your changes to files over time, allowing you to keep track of everything, go back to a specific version of your file, or recall what was actually changed between two versions of it. This also enhances transparency, allowing others to track your work's evolution and collaborate while keeping a clear record of each change.

As tracking changes through version control systems is not limited to code, it is particularly useful for scientists whose work requires them to write things, whether those things are experimental protocols, code, or text, and whether their priority is to track the evolution of their work, collaborate with others, and ensure reproducibility.

In the past few sections of this course, we focused heavily on text-based communication with your computer. You wrote commands in your terminal, put a few of them in a file to run tasks in batch, and implemented comprehensive shell scripts. Now that we have explored the benefits of version control in managing changes, let's look at how it integrates with text-based workflows and cloud services.


### Introduction to Git

The history of version control systems is long and painful, and there are many options that you may have never heard about:

|**VCS**|**Type**|**Key Features**|**Best For**|**Free to Use**|**Personal Projects/Writing**|**Suitable for Academicians**|**Cloud Support**|**Active Since**|
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
|Git|Distributed|Most popular, fast branching, great for collaboration|General-purpose, open-source projects|Yes|Yes|Yes|Yes (GitHub, GitLab, Bitbucket)|2005|
|SVN (Subversion)|Centralized|Strict access control, linear workflow|Corporate environments|Yes|No|Limited|Limited|2000|
|Mercurial (Hg)|Distributed|Simpler than Git, better for large repos|Large-scale projects like Mozilla|Yes|Yes|Yes|Yes (Bitbucket)|2005|
|Perforce (Helix Core)|Hybrid|Handles large binary files, fine-grained permissions|Game development, large enterprises|No|No|No|Yes (Helix Core Cloud)|1995|
|Fossil|Distributed|Built-in issue tracker, single binary|Small projects needing integrated tools|Yes|Yes|Yes|No|2006|
|CVS|Centralized|Very old, outdated branching|Legacy projects|Yes|No|No|No|1986|

Luckily, Git comes as the clear choice for most of us as a version control system as a free and general purpose version control system supported by cloud solutions and suitable for academic use. We all should thank [Linus Torvalds](https://en.wikipedia.org/wiki/Linus_Torvalds) for it.

Knowing about Git and being able to use it is already enough to manage your files, track your modifications, and maintain an organized history of your work -- not only to manage workflows or shell scripting, but also for your papers, reproducible data analyses, and so on. If you combine Git with Markdown, there is nothing you can't do! For instance, this entire course was prepared using Git and [Markdown](https://www.markdownguide.org/basic-syntax/), a simple mark-up language to create well structured documents that can be ported into anything -- from web pages to PDF documents to Microsoft Word files using a tool like [pandoc](https://pandoc.org/) in the comfort of your terminal environment :)

### Frequently Used Git Commands

Git has a very large number of commands and options, however, to begin using it effectively requires only a few. The following list describes some of the most essential and/or commonly used Git commands.

Many of these commands will help us understand how to generate and maintain a local repository from scratch, but we will also get a copy of the anvi'o repository to demonstrate some functions that can be best demonstrated using a real repository with a colorful history.

Let's talk about them a bit while demonstrating their purpose in the command line.

---

`clone` – Clones an existing repository from a remote server

This command brings a copy of an already existing Git repository on your computer.

```sh
## you can clone a repository that is on the same computer
## you are working on
git clone file:///Users/meren/courses/PFLS

## or you can clone a repository that is hosted at a remote
## server (in this case GitHub, which we will talk about soon)
git clone https://github.com/merenlab/anvio.git
```
 
---

`init` – Initializes a new Git repository in the current directory

Using this command you can turn any directory into a Git repository. Running this command will not do anything apart from than saying "I want this directory to be a Git repository".

```sh
git init
```
 
---

`status` – Shows the current status of the working directory

Running this command in any directory that is a Git repository will show you what is going on in it. This is one of the most frequently used commands in Git. Just like the way when you enter a new directory you typically enter `ls` to see what is in it, when you are working with Git you often check its status to see *modified*, *staged*, *unstaged*, or *untracked* files. It is extremely important for us to go through what you see in different scenarios since the output of this command will help us determine what is it we need to do next.

```sh
git status
```

---

`add` – Stages a specific file for commit

When you create a new file in a Git repository or make changes to an existing file, these actions are not recognized by Git immediately, for which you need to explicitly tell Git that the new file or the changes you've made to an existing file is intentional and you would like them to be a part of the history of your repository. Yes. Every. Single. Time.

```sh
## this will 'stage' changes in this file
git add my_script.sh

## this will 'stage' all changes in the work
## directory and below -- which is not a good
## idea unless you are CERTAIN that everything
## you see in the git status output need to be
## added to the repository.
git add .
```
 
---

`commit` – Commits staged changes with a descriptive message

When you use `git add` to stage one or more files does not mean that these changes are committed to the repository. This command will finalize the staging process. Think of `git add` like selecting items for checkout at an online store. Think of this command as finalizing your purchase.

```sh
git commit -m "a meaningful message that describes the changes"

## alternatively you can run it as the following, in which
## case Git will show you a text editor to explain these
## changes -- once you save and exit the editor Git will
## complete the task, and store your staged changes
## permanently to the repository and associate them with
## your description to be kept in the logs
git commit
```

Just like the way you can group different sets of items for purchase, you can use `git add` to stage some files prior to committing them to the repository with different messages. This is a good practice when working with large number of files and create logical groupings of changes.

Every commit will be assigned a unique 'commit hash' which will enable going back to them if/when necessary. 
You can see 'commit hash' values for each change in the output of `git log`.

---

`diff` – Shows the changes made since last commit

This is another very important command that helps investigating what has changed in a given file since the last commit. Without looking at this output running `git add` or `git commit` on a given file often leads to problems, such as committing unintended changes to a repository.

```sh
## show all unstaged changes
git diff

## show unstaged changes in a single file
git diff my_script.sh

## show differences between two commits
git diff <a commit hash> <another commit hash>
```

---

`checkout` – Discards local changes in a file

If you have changed the contents of a file and you are not happy with those changes, running this command will always bring you back to the last committed version of an unstaged file.

```sh
## this goes back to the last committed version of
## one file
git checkout my_script.sh

## this one will do the same thing for all files in
## an entire directory and subdirectories in it
git checkout .
```

There is no going back from `checkout` and your uncommitted changes will be lost forever when you run it on a file. But there is a nice intermediary, which is `git stash` -- Meren shall demonstrate its utility for you.
 
---

`reset` – Unstages a file from the staging area

This command will remove a file that you staged with `git add`. Your changes will not be lost, it will just be unstaged.

```sh
git reset my_script.sh
```

If you run this one, however, it will reset the working directory and staging area to the last commit:

```sh
git reset --hard
```

---

`log` – Shows commit history

This is a very useful command to see what has been happening in a given repository (and looking at its output demonstrates how important it is to put the effort to include meaningful messages to individual commits).

```sh
git log

## this one will condense the output (great for `grep`):
git log --oneline

## and this one will show a fancy output with the inclusion
## of 'branching' events that took place
git log --graph --oneline --all
```

---

`revert` – Creates a new commit that undoes the specified commit

If you have committed changes to a repository and you are no longer happy with those changes, you can always go back to a specific time-point in the history of the repository.

```sh
git revert <commit hash>
```

---

`branch` – Enables working with branches in the repository

One of the most powerful Git commands that makes many many people work on the same repository without stepping on each others' toes. You can create branches in Git, and make changes in them without affecting other branches. Every time you change your branch, the latest committed versions of the files would be stored in your work directory. Better than magic. 

```sh
## this will list all branches in a repository
git branch

## this will create a new branch
git branch a_new_feature

## this will bring all the changes from another branch
## into the active branch, essentially synchronizing
## the active branch to the changes that happened in
## another branch
git merge main

## this will go from one branch to another
git switch main
```
 
---

`blame` - Track who made changes to each line of a file

Running this command on a file will show line-by-line changes to it along with who made that change and which commit brought that change into the repository. Blame is a harsh word, and it is actually often used to find out 'who broke things' in collaborative repositories, but in a more collegial and happy ways :)

```sh
git blame my_script.sh
```

---

`pull` - Get the latest updates (if you are working with a remote repository)

If you are working on a copy of a remote repository, running this command will bring into your local copy the latest changes committed to the remote.

```sh
git pull
```

If you also made changes to files that are changed in the remote repository, you will have to reconcile them. It is not fun (and we will do it now interactively).

---

`push` - Send your changes (if you are working with a remote repository)

Running this command will send all the local changes you've made and committed to your local repository to its remote origin, essentially making your local changes available to everyone who is tracking the remote repository.

```
git push
```

If there has been changes in the remote repository since your last `git pull` you will have to re-run it to make sure you are not overwriting existing changes up there.


### Introduction to GitHub

Git is an open-source, stand-alone program that enables you to create local repositories for version control -- local repositories you store on your own computer or the University servers, and that you can connect with your terminal.

GitHub, on the other hand, is a platform that enables you to store your Git repositories on the 'cloud.' It is a proprietary and for-profit developer platform that has over 1 billion dollars in revenue while itself is not open-source. Unlike Git, which is one of us, GitHub is one of *them*, if you will. But even though Linus Torvalds [hates](https://news.ycombinator.com/item?id=36123124) GitHub for various technical reasons, and I personally hate it for the gargantuan enterprise it represents, GitHub has successfully captured the attention of millions of developers worldwide as it made software development more accessible, manageable, interactive, and fun by adding key features to their cloud-based hosting service such as (1) graphical user interfaces for repository management and code reviews, (2) a service for 'pull requests' that help team members to [review and discuss code changes](https://github.com/merenlab/anvio/pull/2155) before merging them, (3) providing a platform to [report issues and discuss them](https://github.com/merenlab/anvio/issues/1248), (4) [action and workflow](https://github.com/merenlab/anvio/actions) support for continuous integration or testing, and many, many more. GitHub is free unless you wish to have private repositories, in which case you are asked to pay a fee (which is about 4 Euros as of 2026).

In addition to GitHub, there are other online services built on Git, such as [GitLab](https://gitlab.com/gitlab-org/gitlab), which is very similar to GitHub, or [Bitbucket](https://bitbucket.org/lbi-usp/anvio_refine_bin/src/master/), which is primarily used for private repositories. All of them enable you to store your Git repositories in the cloud, allow multiple people to work on the same project with access control, and manage issues and bugs.

We will stick with GitHub due to its convenience.

### Working with Git and GitHub

Your work environment should be ready to work with Git, and you should have a working username for GitHub. Let's go step by step some of the basic Git commands and discuss how to work with GitHub.

#### Setting Up Git

The first thing is to configure your identity, which is something you will only do once for each computer you are using.

```sh
## Set your username
git config --global user.name "A. Murat Eren"

## Set your email
git config --global user.email "a.murat.eren@gmail.com"
```

This information will be associated with your commits to your Git repositories.

#### Creating a New Repository

A **repository** (repo) is a directory that Git tracks.

```sh
## First create a directory on your computer where you
## will keep all your Git repositories. I usually like to
## keep them in a directory called `github` under my home
## directory, so let's create one for you:
mkdir -p ~/github

## and enter into that directory:
cd ~/github

## now create a new project directory, which will keep
## all our files that will be a part of this project.
## I will call it PFLS so you can store your exercises
## in it:
mkdir PFLS && cd PFLS

## Initialize this directory as a a Git repository:
git init
```

The last command tells Git that this directory is meant to be a Git repository, and so it creates a `.git` directory to store its metadata. There is nothing to see in the directory at this point, but if you were to type `ls -a` you can see the hidden Git directory.

#### Checking Repository Status

At any given time you can see which files are being tracked, and what has changed:

```sh
git status
```

You will use this command over and over again. Just so you know, Git keeps system-wide settings in another hidden file in your home directory at `~/.gitconfig`. You can edit this file to add aliases for some commands. Here is how mine looks like just for your reference:

```ini
[user]
	name = A. Murat Eren
	email = a.murat.eren@gmail.com
[alias]
	a = add
	s = status
	st = status
	ci = commit
	b = branch
	co = checkout
	re = remote
	d = diff
	dc = diff --cached
	lol = log --graph --decorate --pretty=oneline --abbrev-commit
	lola = log --graph --decorate --pretty=oneline --abbrev-commit --all
	ls = ls-files
	lg = log --color --graph --pretty=format:'%Cred%h%Creset -%C(yellow)%d%Creset %s %Cgreen(%cr) %C(bold blue)<%an>%Creset' --abbrev-commit
	lgi = log --color --graph --pretty=format:'%Cred%h%Creset -%C(yellow)%d%Creset %s %Cgreen(%ci) %C(bold blue)<%an>%Creset' --abbrev-commit
[color]
	ui = auto
	branch = auto
	diff = auto
	grep = auto
	interactive = auto
	status = auto
[core]
	editor = vim
```
{% include CODEBLOCKFILENAME filename="~/.gitconfig" %}


#### Adding Files to Git

Let's add our first file to this new repository:

```sh
## creating a new, empty file:
touch test.md

## let Git know that we wish this file to be a part of this repository:
git add test.md

## alternatively you could have asked Git to add every file in a
## given directory:
git add .
```

#### Committing Changes

A **commit** saves the changes in a given repository. You need to include a message that describes your changes. In an ideal world, these should be brief and clear:

```sh
git commit -m "my first commit!"
```

After this commit, please edit this file, and check its status with `git status`.

#### Viewing Commit History

At any given time you can review the commit history of a given repository:

```sh
git log
```

#### Associating Local Repo with a Remote One

Essentially you already have a working Git repository. You can leave things as is, and use this repository on your computer forever without having to ever send it anywhere. But you can also 'add' a remote repository to your local Git repository to link it to a remote repository, such as one on GitHub. This is something you do once for each newly generated repository:

```sh
git remote add origin git@github.com:meren/PFLS.git
```

This will allow you to push to and pull changes between your local and remote repositories. Assuming you have added your authentication keys properly, now you can run this command:

```sh
git push -u origin main
```

... and get an error from GitHub that says "*ERROR: Repository not found*". This is happening since as far as GitHub is concerned there is no such remote repository. To mitigate this you can [create a new, empty repository on GitHub](https://github.com/new) called PFLS first, and then re-run the last command, and you should get an output that looks more or less like this:

```
Enumerating objects: 3, done.
Counting objects: 100% (3/3), done.
Writing objects: 100% (3/3), 219 bytes | 219.00 KiB/s, done.
Total 3 (delta 0), reused 0 (delta 0), pack-reused 0
To github.com:meren/PFLS.git
 * [new branch]      main -> main
branch 'main' set up to track 'origin/main'.
```

which means we now have a local repository linked to GitHub! We can add files to our local repository, and push them to the remote repository, and pull changes from our online repository (made by us from other computers, or made by others we gave permission to access that repository) to our local one.

#### Pushing Changes to GitHub

The 'push' command sends all committed changes to your remote repository:

```sh
git push
```

#### Pulling Changes from GitHub

The 'pull' command fetches all committed changes from a remote repository:

```sh
git pull
```


#### Cloning an Existing Repository

You can also 'clone' remote repositories on your computer using Git, which we will use in a second for an actual purpose.

But here is an example:

```sh
## got to the good directory
cd ~/github

## clone a repository
git clone https://github.com/merenlab/reads-for-assembly.git
```

There is much more to Git than what we covered here, but I hope the explanations above gave you the first step into this dimension of computation you can integrate into your daily work and enhance the reproducibility and robustness of your data-enabled research activities. Please take a moment to quickly skim through the following resources to learn more

- Official Git Documentation: [https://git-scm.com/doc](https://git-scm.com/doc)
- GitHub Guides: [https://docs.github.com/en/get-started](https://docs.github.com/en/get-started)
- Git Interactive Tutorial: [https://learngitbranching.js.org/](https://learngitbranching.js.org/)

### EXC-003

You now have the theoretical knowledge required to write a shell script to solve a given problem, AND put it on GitHub to share it with the world.

I would like you to write a shell script that accepts any FASTA file as a parameter, and produces a short report about it that looks **exactly** like this,

```
FASTA File Statistics:
----------------------
Number of sequences: ______
Total length of sequences: _____
Length of the longest sequence: ______
Length of the shortest sequence: ______
Average sequence length: ______
GC Content (%): ______
```

Except where the blanks are, of course, which should be replaced with actual numbers your script will calculate for a given FASTA file. Most of the metrics should be self-explanatory, but for anyone who doesn't know, the GC-content is the proportion of 'G' and 'C' nucleotides in a set of sequences (represented as a percentage).

Once you are done with your script, please commit it to a GitHub repository called `PFLS`, and put it in a directory called `EXC-003`, with the file name `fasta-file-processor.sh`.

This is extremely important, because I will test your solutions by cloning your repository on my computer the following way,

```sh
git clone https://github.com/_______/PFLS.git
```

and testing it on a FASTA file of my choosing using a command like this:

```sh
bash PFLS/EXC-003/fasta-file-processor.sh test.fa
```

You can use this information to make sure things will work on my end (if `cd ~/github && ls PFLS/EXC-003/fasta-file-processor.sh` does not produce an error, you're good). You can test your program using any FASTA file, including those that we have in the data package, but you do not know the *FASTA files I will be using* to test your script on my end.

If the `git clone` step and the step of running your script on a FASTA file both work with the expected output format without a problem, you shall get a full points for this exercise *even if the numbers are not correct*, but if your script produces an output that does not match to the template shown above, you shall get 0 points.

---

A small tip! You can use the AWK function `gsub` to calculate the G and C bases in a given sequence:

```sh
echo "ATCGATCGCG" | awk '{gc_count += gsub(/[GgCc]/, "", $1)} END {print gc_count}'
```

But of course what you want is GC-content, which essentially is the ratio of GC bases to _all_ bases.

The solution I wrote almost entirely in AWK is [here](solutions/EXC-003), and usual, we will go through it together once you are done. Good luck!

{:.warning}
Please note that if your script captures the output of various commands and assign them to variables using `$()`, printing those variables in different operating systems may yield minor differences (such as the inclusion of `\t` or `\n` only in some systems). If you wish to use `$()`, please consider being more defensive to prevent them.

#### Evaluating EXC-003

Since you are such experts of BASH programming at this point, we will have some fun and interactively write a BASH script altogether to test whether *your* script produces an output that is exactly matching the output requested above. Our script will,

* Get a copy of your repository from `https://github.com/$username/PFLS.git` where `$username` will be the GitHub username for each one of you,
* Expect to find the directory `EXC-003/fasta-file-processor.sh` in the cloned repository,
* Capture the output of your program when it is run on a test FASTA file, and finally,
* Test whether the output produces matches to the output expected.

For the last step, we will make use of the following logic (the code is incomplete, but you *could* complete it to make sure it will give a green mark for your submission 😉):

```sh
echo "$output" | awk -v username="$username" '
    BEGIN {
        FS="\n";
        output_is_correct = 1;
        expected[1] = "FASTA File Statistics:";
        expected[2] = "----------------------";
        expected[3] = "Number of sequences: [0-9]+";
        expected[4] = "Total length of sequences: [0-9]+";
        expected[5] = "Length of the longest sequence: [0-9]+";
        expected[6] = "Length of the shortest sequence: [0-9]+";
        expected[7] = "Average sequence length: [0-9]+(\\.[0-9]+)?";
        expected[8] = "GC Content \\(%\\): .*$";
    }

    {
        if($0 !~ expected[NR]) {
            output_is_correct = 0;
            exit
        }
    }

    END {
        if (output_is_correct && NR == 8)
            print("✅ \033[4m" username "\033[0m: The output is correct! 😊")
        else if (output_is_correct && NR != 8)
            print("❌ \033[4m" username "\033[0m: The number of lines in the output is wrong 😞")
        else
            print("❌ \033[4m" username "\033[0m: Line " NR " is wrong 😞")
    }'
```
{% include CODEBLOCKFILENAME filename="check-exc-003-output.sh" %}

The code above includes some ideas we haven't discussed, such as the [tilde operator](https://www.gnu.org/software/gawk/manual/html_node/Regexp-Usage.html#index-_007e-_0028tilde_0029_002c-_007e-operator) (`~`) in AWK, which is used for *pattern matching* with regular expressions, or the use of 'arrays', which I will describe while we are going to interactively turn this snippet into a complete program to will go through submissions from each one of you automatically.

## Introduction to Large Language Models (LLMs)

LLMs like ChatGPT, Claude, Gemini, Grok, DeepSeek, and [others](https://en.wikipedia.org/wiki/List_of_large_language_models) are AI systems trained on vast amounts of text data to understand and generate human-like text. They are based on the 'transformer architecture', a breakthrough in deep learning that enables them to process and generate coherent and relevant text.

LLMs are generally trained on very diverse datasets, including books, news articles, Wikipedia, open-access scientific literature, code repositories, and more. Large and diverse training data allows LLMs to generalize concepts across domains and make connections that may be elusive. For instance, you can ask a general LLM a question like "can you describe a critical phenomenon, best practice, or risk that is shared between software development and modern agricultural practices?", and you *will* get an answer that will mostly likely be more satisfactory than an answer you may get from a person who is only a computer scientist or only a farmer. Once a foundation model is trained on general data, it can be fine-tuned to perform specific tasks with high performance, such as generating images, videos, code, or creative writing. Then the users of these models (whether people, or other programs or models) can interact with LLMs through 'prompts', i.e., the input text you put in into the text box that guides the model to generate the relevant output for you. Even though an average person hears a lot of about OpenAI/GPT (due to its success and convenient online 'chat' interface) or DeepSeek (thanks to its surprising benchmarks and push towards open source), there are already [many many models](https://huggingface.co/models) for general or specific applications.

### Transformers: How do they even?

Transformers are the foundation of modern LLMs. As you can imagine, it is almost impossible to fully describe how they work here since truly understanding transformers require a substantial understanding of linear algebra, probability and statistics, and concepts in machine learning for starters. But luckily we don't need to understand all the details here to have a general sense of what they do for us. You can drive a car without knowing anything about combustion engines, but knowing just a bit about the engines is the only way for one to understand why there is a nonlinear relationship between the speed of a given car and the fuel consumption, or why the fuel economy often drops significantly at very high speeds even if one ignores aerodynamic drag. The same applies to LLMs. You don't need to grasp every technical detail, but without even a basic understanding, the behavior of LLMs will come across as a mystery rather than engineered solutions with predictable shortcomings and limitations.

Transformers were first introduced in 2017 by a team of researchers affiliated with [Google Research](https://research.google/) through a seminal paper that was very aptly titled as "[Attention is All You Need](https://arxiv.org/pdf/1706.03762)". Prior to this architecture, models that were capable of 'natural language processing' were struggling with key limitations, such as difficulty capturing long-range dependencies in lengthy texts. As they processed each word one after another, pre-transformer models struggled with tasks that required understanding relationships between words that were far from one another in a sentence. The transformer architecture changed all that by introducing self-attention, which allows the model to process all tokens in a sequence in parallel, instead of processing tokens one by one. This made it possible to fully utilize GPUs, which excel at parallel matrix operations. Even though the paper was published less than a decade ago, it is already cited over 200,000 times, which is a meaningful metric in this case to understand how it revolutionized natural language processing. Here's a heavily watered down explanation of the core concepts of transformers with the hope that they provide you with some abstract ideas about the inner-workings of the transformer architecture that makes modern LLMs a reality.

#### Input Representation

Upon receiving an input prompt, transformers do three things:

* Tokenization: The input is first broken into smaller units called tokens. These tokens could be words, subwords, waveforms, or pixels depending on whether the model works with text, sound, images, or video. Let's stick with words for the sake of simplicity here. Think of the tokenization process like cutting a sentence into individual puzzle pieces. For instance, if the input is "Write a Python function to calculate the GC content of a DNA sequence", the tokenization step will split it into tokens ["Write", "a", "Python", "function", "to", "calculate", "the", "GC", "content", "of", "a", "DNA", "sequence", "."] for further processing, where each of these elements will have a unique index that connects them to the word recognized by the model (just to note for my own sanity, in reality LLMs do not often tokenize words as whole words but as subwords; so Python in this example may be tokenized as "Py" and "thon" as they always try to properly handle rare words efficiently).

{:.notice}
For a given model, the term token is also used to describe the training data size used to generate the model (i.e., models process many many tokens to learn patterns) or the capacity of it (since models generate responses one token at a time and they are limited how many tokens they can handle at once). For instance, based on numbers I found online at the beginning of 2025, GPT-3 model was trained in 300 billion tokens and can work with 128 thousand tokens at a time. In contrast, DeepSeek R1 model was trained on 671 billion tokens. The tokens in tokenization and tokens in token counts of models are coming from the same underlying concept: a unit of processing for the model. But they are not identical.

* Embeddings: Once the tokenization is done, each token is then converted into a numerical vector (a list of numbers) called an *embedding*. For instance, the word "Python" in the previous example might be represented as a vector that goes like [0.22, -0.3, 0.7, (so on)]. Each dimension here will be describing a linguistic feature of a given word, and the number of dimensions in the embedding, which is also called the embedding size, will depend on the model. For instance, in Google's BERT model the embedding size is 768 (i.e., there will be 768 numbers in that vector to describe the word "Python"), and in GPT-3 model the embedding size is 12,288. In a way, these vectors represent the meaning of the token in a way the computer can understand and work with it. You can think of this as translating a word from human language to model language by placing it at a unique coordinate in a hyper-dimensional universe of known tokens. These vectors make it possible to perform mathematical operations on words, or resolve the relationships between a given word and other words. Embeddings are learned during the model's training phase, where the model processes vast amounts of text to identify patterns and relationships between words. These initial embeddings remain fixed until a specific input sentence is processed, at which point self-attention dynamically refines them. So the initial embeddings associate each token with a value in that space, where the Python in our example does not know whether it is a programming language or a snake. This is where the self-attention mechanism comes in. As the model processes the full sentence, it updates each word's embedding dynamically using surrounding words. For example, since the phrase "write a function" strongly suggests programming, the embedding for "Python" shifts toward the region of the embedding space where programming-related terms like "Java" or "C++" are located. At this stage Python as a programming language and Python as a kind of snake will have entirely different vectors. This contextual adaptation is exactly how transformers revolutionized natural language processing, where meaning of words are dynamically refined based on context rather than a fixed dictionary-like mapping, and made it very similar to how you understand *Jrxaal must be a programming language* when you read the sentence "I need a Jrxaal function that could calculate Celsius from Fahrenheit". Even though you have never heard about Jrxaal before, and it doesn't exist. Refined embeddings following the self-attention step are not just used for meaning, but the new meaning they have gained influence the final outputs when they are passed to further layers of the Transformer to perform a task.

* Positional Encodings: A step that takes word embeddings, and turn them into position-aware word embeddings. This is necessary since unlike its predecessors that worked with each embedding one by one, transformers take all embeddings all at once to process them, which improves their performance dramatically, but at the expense of losing the original order of words. The step of updating the embeddings with positional encodings ensure that the sequence of tokens, which may dramatically influence meaning, are preserved. Think of this like numbering the puzzle pieces so you know which piece comes first, second, etc. This is a step that is done so very elegantly in the papers that introduce transformers, and I had very hard time truly appreciating the nuances there until I saw [this excellent YouTube video](https://www.youtube.com/watch?v=dichIcUZfOw) (if you are really interested in better understanding this step). This is one of the two very big innovations in the transformer architecture since it doesn't increase the input data size while maintaining the order information for parallel processing that enable significant gains in speed.

#### Self-Attention Mechanism

Self-attention mechanism is one of the biggest innovations in the transformer architecture, which allows the model to focus on the most relevant words in a sentence when processing a specific word, and has become a central aspect of every modern deep learning model. In a sense, this is the stage where the model takes each puzzle piece one by one and compares it all the other pieces in each step to resolve which combination makes the most sense through mathematical operations. For instance, in the example prompt "Write a Python function to calculate the GC content of a DNA sequence", when the model is processing the word "GC", the self-attention procedure will help the model to focus on "calculate" to understand the task at hand, "DNA" to understand the context, and "Python" to know that it needs to generate code, and assign very high attention scores to "calculate" and "Python" since  they specify the task and output format, "DNA" might get a medium score since it provides context. At the end of this step, every token will have sets of attention scores, and thus have a contextualized representation with information incorporated from all other tokens before it is passed on to the next step. Here is [an excellent article](https://sebastianraschka.com/blog/2023/self-attention-from-scratch.html) to understand this mechanism better if you are interested.

#### Multi-Head Attention

If self-attention mechanism is a dynamical spotlight shining on a stage where the model focuses on different actors (i.e., tokens) based on their relevance, multi-head attention is an extension of it where multiple spotlights shining at the same time, each focusing on different aspects of the scene, helping the viewer (i.e., the model) to develop a more comprehensive understanding of what is happening. If it helps, you can also think of this as having multiple teams working on different parts of the puzzle simultaneously, each focusing on a different aspect (e.g., one team looks for edges, another looks for colors, et). In our example of "Write a Python function to calculate the GC content of a DNA sequence", one attention head might focus on the relationship between "GC" and "DNA" to resolve the biological context, another attention head might focus on the relationship between "Python", "function" and "calculate" to resolve the coding task. Combining outputs from these heads with multi-head attention, i.e., a coding task for a biological question, might represent a more nuanced understanding of the input task or data or more diverse relationships in them (Meren found [this blog post](https://sebastianraschka.com/blog/2023/self-attention-from-scratch.html) extremely helpful to better understand some aspects of multi-head attention).

####  Feed-Forward Neural Networks, and other steps

After the attention mechanism, the output is passed through a feed-forward neural network (FFN). This is like a second layer of processing that refines the information further. Imagine taking the assembled puzzle pieces and smoothing out the edges to make them fit perfectly. Going back to our example prompt of "Write a Python function to calculate the GC content of a DNA sequence", the FFN processes the combined output from the attention heads, refining the understanding of the task. For example in this stage it might be reinforced that "GC content" refers to the proportion of "G" and "C" nucleotides in a DNA sequence, and "function" refers to a Python function.

Other steps that are even more machine learning heavy include *Layer Normalization* (ensuring the data through the model stable, such as ensuring that the puzzle pieces don't fall of the table), *Residual Connections* (ensuring model doesn't forget important things, such as keeping a reference picture of the puzzle to guide your efforts relevant), and *Stacking Layers* step to create a deep network to find complex patterns and relationships in the data, such as having multiple puzzle solving teams operate together to refine the solution further and further. In the first layer, the model might focus on understanding the task ("calculate GC content"), in the second layer it might focus on the context ("DNA sequence", in the third layer it might focus on the output format ("Python function"). By the final layer, the model may have a clear understanding of the prompt and be ready to generate the output.


#### Output Generation

LLM's ability to generate meaningful output is as complex as their ability to make sense of the input, and relies on similar principles such as keeping track of the context and making use of the learned patterns from existing text with which they were trained. Transformer uses a decoder architecture for text generation, precisely generating tokens one at a time using previously generated tokens as context. This is almost like you writing down something new: carefully considering each word in a sentence you are writing while keeping in mind what you have already written, and relying on your understanding of the language in which you are writing based on all the previous readings you have done before to benefit from patterns that look similar to what you are writing. So both the context and what is previously written play a role in your text generation, as do in LLMs.

For our input **"Write a Python function to calculate the GC content of a DNA sequence."**, the model does not generate the entire output in one go. Instead, it uses linguistic properties of the Python programming language to establish that the first token must be `def`. Then, it considers `def` as a context to predict the next most probable token, e.g., a function name, such as `calculate_gc_content` given the task, and then the `(`, and so on, building the output token by token. As this program in its final form that is emerging as it is being generated likely does not exist anywhere, each step assigns probability scores to possible next tokens, during which the model selects the most likely option, and sometimes introduces controlled randomness for creativity, which is controlled by a parameter generally called 'temperature'. As a user of general-purpose services such as ChatGPT or DeepSeek you don't have access to parameters such as temperature, 'nucleus sampling', 'frequency penalty', or others, which are set to certain default values in these services. But if you were to use LLMs from within other programs through what we call 'application programmer interfaces' or APIs, then you would get to do a lot of prompt engineering to see the impact of these parameters on output generation. Here is [an amazing resource](https://www.promptingguide.ai/) on prompt engineering that can help you understand a lot of these concepts in greater detail if you are particularly interested in these things.

At the end of the continuous generation of text that finds its way by choosing the most probable continuation after each token, the model reaches its stop condition with a complete result, which looks like this for our example prompt:

```python
def calculate_gc_content(sequence):
    gc_count = sequence.upper().count("G") + sequence.upper().count("C")
    total = len(sequence)
    return (gc_count / total) * 100
```

Depending on the key settings of the model, running the same question multiple times may not result in the same output, which is a desired feature of chatbots to give the impression of a level of dynamism that make them sound like jazz musicians rather than robots that read the same canned responses. 

The same principles in output generation apply to code generation, translating from one language to another, creative writing, or generating arguments for or against certain ideas, scientific or not. Thus, as far as LLMs and the transformers that bring them to life are concerned, there is no difference between writing Python code or translating from English to French except a few minor nuances.


## Using LLMs for Programming

The emergence of LLMs has been one of the most fascinating developments in computation given their broad range of applications. One of the by-products of these powerful tools was to make computers write code for us. An ability now democratizes programming at scales hard to comprehend since even your uncle can use ChatGPT to generate computer code in any language.

Making computers write code for humans has long been a dream of computer science, and it has a fascinating story with roots go way further back than you probably anticipate. For computers to generate code, first there had to be a way to formalize the rules of languages in which code can be written. In 1956 Noam Chomsky, a prominent linguist of that time and one of the most brilliant philosophers of our time, introduced formal grammars and the idea that grammars of languages can be formalized to understand the languages they generate. This concept, which is now known as [Chomsky Hierarchy](https://en.wikipedia.org/wiki/Chomsky_hierarchy), rapidly became a cornerstone of linguistics. Probably Chomsky had no such intentions at the time, but his mathematically rigorous effort to try to understand the underlying principles of human languages had enormous implications in computer science. The core idea captured by the Chomsky Hierarchy provided a classification of languages that would later influence how computers can parse and generate code, and led to the design of programming languages and compilers that could work with them. So in retrospect, we can say the human's desire to make machines write code goes as far back as 1950s. A lot has happened between then and now, and I wish I could talk more about UNCOL, PROLOG, biology-inspired genetic programming, and all the other exciting and inspiring history of our effort to bridge the gap between humans and computers. Perhaps another time.

It is safe to say that compilers, which take instructions written by humans in a programming language and turn them into machine instructions that are almost impossible to understand by humans, is already a form of code generation capability. Before compilers, most programming was done by writing machine code, i.e., *literal* instructions for the processor, by directly invoking opcodes and memory addresses. Since you will never see it anywhere throughout your journey in programming, I wanted to show you an example. For instance, the following is what a machine code that would calculate the sum of a bunch of numbers in an array on Intel 8086 CPU that was introduced in 1978:

```
10111000 00000000 00000000  ; MOV  SI, 0       (Index SI = 0)
10111011 00000000 00000000  ; MOV  BX, 0       (Sum BX = 0)
10110000 00000101           ; MOV  AL, [SI]    (Load array value)
00000011 11000000           ; ADD  BX, AX      (Add to sum)
01000110                    ; INC  SI          (Next element)
10000011 11110110           ; CMP  SI, 6       (Check end condition)
01110100 11111000           ; JNE  loop        (Repeat if not done)
```

The first programming language for which we had a high-level compiler was FORTRAN. FORTRAN took away much of the pain you see above by generating that kind of machine-level code from a language that is much easier to read and understand (well, relatively speaking):

```fortran
PROGRAM SUM_ARRAY
INTEGER :: SUM, NUMBERS(6)
DATA NUMBERS /1, 2, 3, 4, 5, 6/

SUM = 0
DO I = 1, 6
   SUM = SUM + NUMBERS(I)
END DO

PRINT *, "Sum is:", SUM
END
```

In fact, the CPU technology has advanced so much from the days of relatively simpler architectures such as Intel's x86 processors (which had less than 100 opcodes in contrast to over 1,500 we have in modern CPUs), it is practically impossible to write machine code anymore as the difference between doing it back then and now is akin to the difference between flying a kite and flying an F16 fighter jet. We depend on compilers that generate highly optimized, architecture specific code that can take advantage of a vast number of features new processing units entail.

Compilers were our first successful attempt to write human-understandable instructions and generate code that computers could execute, which started this ball rolling. That said, given your fresh exposure to the syntax that is necessary to write shell scripts, you can tell that even the human-understandable instructions are not equally or easily understandable by all humans. To truly unlock the power of computers for everyone, our journey of generating code needed to go from generating code for machines from instructions written in programming languages to generating instructions written in programming languages from ideas explained in human languages.

Before the advances in natural language processing that made it finally possible, there has been many many attempts to bridge that gap. These efforts included simplest imaginable solutions such as syntax highlighting or auto-completion and documentation lookups to slightly more advanced approaches such as boilerplate code generation in integrated development environments (which never really useful), tools that translated code from one programming language to another (which never really worked well), tools to create templates of common programming patterns (🤮), and even more sophisticated efforts in recent times with modern machine learning approaches that did not use LLMs. Perhaps one should look back fondly to those attempts since they were certainly helpful, but they did lack context awareness, had no abilities to understand human languages, could not adapt or learn, could not explain what they were doing, and they had very limited code generation capabilities. Needless to say our initial attempts to generate code never reached the transformative impact of coming up with compilers, and pale in comparison to what we have today with LLMs that deliver everything we wished for and more.

Most life scientists are familiar with the power of molecular tools to generate data to address complex biological questions, and analyze, interpret, and visualize biological data. Programming is nothing more than that even though it can often feel intimidating to those without a background in computation. I hope throughout that course we addressed that to a degree, and you can recognize LLMs as powerful 'assistants' that can make programming more accessible and intuitive for you. But one must first understand advantages and disadvantages of using LLMs for programming, and aspects of it one should pay attention to when incorporating them into their daily work.

### Advantages of using LLMs for programming

Even though advantages are likely extremely clear, let's just mention a few of them for the sake of being explicit.

* LLMs lower the barrier to entry to programming. As they can generate code snippets and explain core concepts of programming in plain language, they make it easier to start writing code and troubleshoot code without extensive training.

* They allow you to describe what you need in your own language, and produce code that accomplishes that. This is great for rapidly generating code to test ideas, and lower the burden of non-creative tasks (such as calculating the GC content of a piece of DNA or finding the longest sequences in a FASTA file, etc). In fact, everything we have done so far in this course can be done with an LLM effortlessly. This way they enable you to focus on science rather than programming best practices or learning about programming paradigms that require so much experience.

* LLMs can also act as great trainers. You can copy-paste a piece of code, and ask for line-by-line explanations. It is not a useful practice for someone who has no idea about programming, but with a relatively low level of familiarity with programming concepts, interactively working with LLMs can be very helpful to learn more. The best aspect of all is that thanks to their extensive training on existing code, LLMs can help you write code that adheres to the highest standards of coding. 

### Disadvantages of using LLMs for programming

Nothing comes without their trade offs and LLMs are no exception. 

* They will make you stupider if you are not careful. Relying on LLMs for everything will prevent you to develop programming skills that will rewire your brain to help it hold bigger and more complex ideas. You will likely fail to learn fundamentals of programming, which in return will prevent you from being able to troubleshoot issues and handle new challenges.

* LLMs make a lot of mistakes. Generate code with them without fully understanding the underlying principles will force you to blindly trust their outputs.

These disadvantages can be mitigated with a mindful approach.

* Learning the basics of programming like we attempted here and understanding the utility of variables, loops, conditionals, functions, etc will help you run, debug, update, and integrate code suggestions from LLMs.

* Combining LLMs with old-school learning practices, and continue forcing yourself to come up with solutions or seek for answers in traditional ways (such as looking things up in help menus and tutorials) will help you maintain your control over the tasks you use LLMs for.

LLMs are revolutionary, and by managing to born at the right time, you have the rare opportunity to witness this transformation firsthand. Imagine all the people who have lived and died without seeing what humans managed to create. But you have not been born in an era where it may be acceptable to let these advances make you redundant.

### Broader implications of the existence of LLMs

I also want to say a few words about this new world order in a different way, just to create space for those of you who may wish to engage in deeper thought experiments.

The transformer architecture fits on a napkin. You don't need to be a mathematician to figure that out since you can tell that it is the case by literally looking at the length of the paper that describes the core principles behind it that everyone relies upon. So you have the algorithm to create a new generalized model. Next, you need the data for training. It is right there, too: while massive, all the data you need can be scraped from the internet. That's what everyone else did. But if it is that simple, why are all the famous models produced and put in the market by such a small number of well-known profit seeking entities such as OpenAI, Google, Meta, and a handful of others? How do they keep this transformative power so consolidated? The answer is simple: the barrier to entry is pure capital.

To actually train a model like GPT-4, you need, and I am not exaggerating, a dedicated power plant, a warehouse full of the most advanced chips on Earth, a cooling infrastructure that could climate control a small city, and about a billion dollars for starters. Public data show that GPT-4 takes about 230,000 petaFLOP-days to train (which is about 1 billion laptop-days, if you will). What does that really mean? Well, it means that to be able to train such a model, say in one month of 24/7 work so there aren't five other models by the time you are done, you would need a computational infrastructure that can maintain about 7.7 exaFLOPS of compute power around the clock. For comparison, the largest supercomputer in Germany, [JUPITER](https://www.fz-juelich.de/en/jsc/jupiter), can reach about 1 exaFLOP at its peak. Which means the best Germany can do, and one of the most advanced compute resources in Europe that serves a few thousand researchers, has only a fraction of the compute power needed to train a single GPT-4-class model in a reasonable timeframe even if everyone else stopped using it for other purposes. The compute resource need is no joke. To build a system that can achieve and maintain 7 exaFLOPS, you would need around 20,000 NVIDIA H100 GPUs, which would cost you about a billion Euros at today's prices, if you can actually buy them. As you know, these chips are now so sought-after that they have become a source of geopolitical tension that result in export restrictions and waiting lists for large orders spanning years.. But even if you could buy them, you would realize that it is not free to operate them. First, you would need a physical space that is equivalent to a large IKEA store to fit all your server racks filled with GPUs, and your cooling systems, cables, networking equipment, and storage, which would collectively set you back another few hundred million Euros. The network bandwidth requirements of your new IKEA warehouse to handle the constant data exchange between all GPUs during training would be so high that it would surpass the entire internet infrastructure of a small country, and your GPUs would consume about 30 to 40 megawatts of electricity constantly, an amount of energy that is enough to power a city the size of Oldenburg. You would have to have all this in place before you start training your model, and that is before counting the monthly expenses for cooling, electricity, and specialized engineering staff.

This is why only a handful of organizations can build frontier models with deep reasoning capabilities and keep them up-to-date. The barrier to building these models is not intellectual. It is industrial.

And here you must ask yourselves, not only what this all means for a more equitable future for all humans, but also to what extent you wish to outsource your reasoning capabilities and ability to experience and make sense of life and its challenges to these costly 'products'. Products that must become indispensable to justify their cost (to you and the society), to sustain their unquenchable hunger for resources, and to feed the next generation of billionaires.

### Some questions to think about in the post-LLM era

The following is a list of questions I have received through the class citizenship emails during 2026, we discussed it together during the class, and I think it is a good idea to continue wondering about these, even in cases where we have no clear answer, and learn to sit with our tension:

* Could reliance on LLMs erode critical thinking skills?
* Should AI assistance with editing be disclosed in academic papers?
* Can the energy and resource demands of AI data centers be justified given climate change?
* How does training data bias toward widely-spoken languages affect LLM performance in others?
* Can LLMs create novelty, or only recombine patterns from training data?

Most of these questions and more are deeply explored in many disciplines of research, including economics, social science, and neuroscience. Some of them (such as the carbon footprint one) are heavily contested, and the jury is still out (or lost, if you will, given the wide range of interest groups involved for some of these topics). But other questions in this list have almost certain answers. One of such questions, and I think the one that is most immediately relevant to us here, is the first one: do LLMs erode critical thinking skills? I will let you guess the answer, and will share with you a few resources to follow up on if you are interested:

* [Your Brain on ChatGPT](https://arxiv.org/abs/2506.08872): This is a longitudinal study that uses EEG to track brain activity while performing writing tasks. After four months, participants who relied on ChatGPT for writing shows significantly reduced neural connectivity in regions crucial for attention, working memory, and language processing. Tthey term this trade off of short-term convenience at the cost of long-term cognitive development 'cognitive debt'. More alarmingly, participants who became accustomed to AI assistance during this phase struggled to re-engage the necessary neural networks in their future writing tasks without AI assistance.

* [Cognitive Offloading Is Real—And New Learners Are Most at Risk](https://compare.rm.com/blog/2025/09/cognitive-offloading-is-realand-new-learners-are-most-at-risk-how-rm-compare-keeps-the-thinking-human/): This is science of learning research synthesis that tells us that the risks of 'cognitive offloading' are not shared equally across all individuals. According the these authors, less experienced learners are particularly susceptible to the negative effects of AI-driven offloading compared to experts. While the experts can benefit from the shortcuts with little risk, for those that are just learning new topics may suffer from limited critical growth, self-monitoring, and resilience. Basically the most critical essential lifelong skills in learning.

* [From Offloading to Engagement](https://www.mdpi.com/2306-5729/10/11/172): This is a slighly positive (or more nuanced) note on the implications of AI on learning. In this experimental study that includes partipants from Germany, Switzerland, adn UK, researches benchmark the implications of AI usage given the *way* participants use it. The take home message from the stdy is that 'unguided' AI use fosters cognitive offloading without improving reasoning quality. But structured prompting (i.e., when users critically engage with model outputs) enhances both critical reasoning and reflective engagement.

While these are worth thinking about, my advice at this juncture is simple and pragmatic: use LLMs to support and advance your critical thinking, be thoughtful and responsible as you make them a part of your work, and turn them into your assistant that you *could* live without, and not your guide, without whom you could not search for answers on your own.

## An exercise using an LLM of your choice

Now that you have an general understanding of the terminal environment, shell scripting, version control, and the fundamentals of how LLMs work, we are going to get our hands dirty with some LLM work using one of the exercises we have already gone through during this course: [EXC-003](#exc-003), where you had to write a script to process a FASTA file and print out some features of it.

Now we will re-develop that solution using AI assistance with the same initial instructions for EXC-003, but extending it further for you to encounter various lessons about robust scripting and effective AI collaboration.

While solving EXC-003 with your favorite LLM, I suggest you explore the following additional considerations to make your script better, more robust, and more generally applicable and useful. For instance, thinking about its robustness and resilience against edge cases, can you figure out the best practices regarding

- What happens with an empty file, or a file with a single sequence?
- How should the script handle sequences containing ambiguous bases (N, R, Y, etc)? Should they be counted? Excluded from GC calculation? What would be the best practice there?
- What if someone passes a file that is not a FASTA at all? How should your script behave, and how should it handle the situation?

You can also consider extending your solution with additional statistics to report to the user of your porgram. Some ideas for that may include the following, but you can come up with more:

- Calculating N50 (a metric you will encounter constantly in genomics)
- Reporting the number of sequences that are above a certain length threshold
- Showing a simple histogram of sequence length distribution in the terminal for convenience

You can also consider usability improvements. For instance,

- You can add a `--help` flag that explains how to use the script
- You can make sure the script can accept multiple FASTA files and produce a combined or comparative report that the user can save as a TAB-delimited text file
- You can add a flag to report an HTML page rather than a TAB-delimited file

Once you are ready, I will have a few volunteers (or randomly selected individuals) come here and present to us their,

* Prompting strategy, so we can see their first prompt and the progression and their ways of refining the code
* Final solution, as they walk us through the code and explain what it does with example runs
* Lessons learned, where they can tell us if they observed anything interesting

So, please go back to EXC-003 now, think about it again from scratch with some help from your favorite tool, and then we will regroup to discuss things together.

## EXC-004

This is our final exercise together. And it is here to immerse you into the joy of writing code for real-world data tasks, and you are now free to use LLMs to do that.

This is the most realistic problem we are going to be working on so far, and comes from a real task in a real project I and Sarahi Garcia worked on together in the past.

Please go into the `EXC-004` directory in your data pack, and take a look at the contents of it. You will realize that

* There are multiple directories in `RAW-DATA` directory that goes like `DNA57`, `DNA58`, `DNA64`, etc. Each of these directories contains the results of a genome-resolved analysis of the metagenomic sequencing of a single culture generated from Lake Erken samples. Since each culture was started with just a few cells from the environment, the researchers who conducted this study were able to recover one or more metagenomic bins. The directory names match to library names rather than sample names, and the actual names of these cultures are stored in the file `sample-translation.txt`. Each of these directories has identical structure, and contain the following two files and another directory:
  - `bins/` -- Contains the FASTA files that emerged from the automatic binning of each sample.
  - `checkm.txt` -- Completion and redundancy estimates for the tentative genomes represented by individual FASTA files in the `bins/` directory based on bacterial and archaeal single-copy core genes.
  - `gtdb.gtdbtk.tax` -- Taxonomy for each of them based on the Genome Taxonomy Database.

One of the most critical next steps is to estimate the actual abundances of individual genomes. But this particular organization of these data is not very useful for such downstream analyses, and we need a shell script that can help us combine these data into a more meaningful representation. The script needs to do the following things:

* Generate a new directory at the same level of `RAW-DATA` called `COMBINED-DATA`.
* Process each FASTA file in individual directories to (1) copy each `bin-unbinned.fasta` into the to `COMBINED-DATA` directory as `XXX_UNBINNED.fa` and (2) copy every other FASTA file into the to `COMBINED-DATA` directory as `XXX_YYY_ZZZ.fa` where,
  * `XXX` is the culture name recovered from `sample-translation.txt`,
  * `YYY` is `MAG` if the completion is 50% or more and contamination is 5% less according to the information in stored in the relevant `checkm.txt` file, otherwise `YYY` is `BIN`,
  * `ZZZ` is `001`, `002`, and so on, where each `MAG` and each `BIN` for individual cultures have sequential numbering.
* Ensure that each sequence in each FASTA file has a unique *defline* associated with the culture (`XXX`) -- pro tip: you can use {% include PROGRAM name="anvi-script-reformat-fasta" %} for this.
* Copy `checkm.txt` and `gtdb.gtdbtk.tax` files in individual `bins/` directories into `COMBINED-DATA` as `XXX-CHECKM.txt` and `XXX-GTDB-TAX.txt`.

Based on these instructions, when you run your script in the directory `EXC-004`,

```sh
bash generate-combined-data.sh
```

or

```sh
python generate-combined-data.py
```

the contents of the newly generated `COMBINED-DATA` directory should look like this:

```sh
ls COMBINED-DATA/
CO64-CHECKM.txt
CO64-GTDB-TAX.txt
CO64_BIN_001.fa
CO64_BIN_002.fa
CO64_BIN_003.fa
CO64_BIN_004.fa
CO64_BIN_005.fa
CO64_MAG_001.fa
CO64_UNBINNED.fa
CO83-CHECKM.txt
CO83-GTDB-TAX.txt
CO83_BIN_001.fa
CO83_BIN_002.fa
CO83_BIN_003.fa
CO83_BIN_004.fa
CO83_BIN_005.fa
CO83_BIN_006.fa
CO83_BIN_007.fa
CO83_BIN_008.fa
CO83_BIN_009.fa
CO83_BIN_010.fa
CO83_MAG_001.fa
CO83_MAG_002.fa
CO83_UNBINNED.fa
CO86-CHECKM.txt
CO86-GTDB-TAX.txt
CO86_BIN_001.fa
CO86_BIN_002.fa
CO86_BIN_003.fa
(...)
```

Once you are done, please commit your script to your GitHub repository for `PFLS`, and put it in a directory called `EXC-004`, with the file name `generate-combined-data.sh` or `generate-combined-data.py`.

The solution I wrote in BASH is [here](solutions/EXC-004), and as always, we will go through it together once you are done.

## Closing Remarks and Final Discussions

This brings us to the end of our 2-weeks journey. Thank you very much for your attention and participation, and I hope this was a useful experience for you.

---
layout: page
noleftpanel: true
title: "Programming for Life Scientists"
author: Course Plan
date: February, 2025
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
* The utility of Python programming language,
* Large Language Models and AI-assisted program solving,
* Orchestrating multiple of these tools to solve complex, daily tasks.

This document includes **everything** you need for this course.

{:.warning}
If you are taking this course as a part of your educational program with the expectation to be graded for it, please pay particular attention to the [Course Responsibilities](#course-responsibilities) and [Grading](#grading) sections.

## Faculty and Communication

The course is directed by [Prof. Dr. A. Murat Eren](https://merenlab.org), who goes by Meren, and [Prof. Dr. Sarahi Garcia](https://miint.org/). The lectures and exercises will be primarily delivered by Meren, but the following table lists individuals who were directly or indirectly involved in the contents of the course:

|Name|Role|Expertise|Contact information|
|:---|:---|:---|:---|
|**Meren**|Professor|Microbial Ecology, Computer Science|meren@hifmb.de|
|**Sarahi**|Professor|Microbiology, Microbial Ecology|sarahi.garcia@uni-oldenburg.de|
|**Florian Trigodet**|Senior Scientist|Microbiology, Bioinformatics|florian.trigodet@hifmb.de|
|**Iva Veseli**|Postdoc| Microbial Ecology, Computer Science|iva.veseli@hifmb.de|

Throughout the course (and beyond) you can reach out to Meren with questions, who should be your first contact for anything related to the activities that will follow.

## Course Facts and Logistics

{:.notice}
The course code at the University of Oldenburg is `5.13.622`.

The design of the course assumes that you know *next to nothing* about any of these topics. High expectations poison everything: so please assume that you will unlikely become an expert of the tools and approaches that are covered here, but hopefully, our discussions will provide you with enough understanding of the fundamentals of programming, helping you think about how to make computers work for you, which will give you enough foundation to become an expert of everything here.If you’re taking this course in person, we’ll spend about 60 hours together. That’s not nothing, but it’s hardly enough time to become truly skilled at anything. In his book "*[Outliers: The Story of Success](https://en.wikipedia.org/wiki/Outliers_(book))*", Malcolm Gladwell, a Canadian author and thinker, argues that a person could become an expert in nearly any field as long as they were willing to devote the requisite 10,000 hours to studying and practicing the subject or skill. Reminiscing his own journey, he says, paraphrasing here from his book, "*I was completely overwhelmed at the beginning, but by the end, I felt like an expert. It took me 10 years—exactly that long*". You may feel overwhelmed at the start, throughout, and even at the end of this course, and that's OK. My goal is to help you recognize that bending computers to your will with freedom is something you can see yourself doing well if you were to invest time. That’s why I encourage you to ask questions, participate in discussions, and see this as more than just another course to get through.

While the course is designed for life scientists in mind, it should be useful for anyone to develop a sufficient understanding of the topics the course aims to communicate. But we will be using mock and real-world datasets and problems that convey typical characteristics of what researchers often encounter in the data-enabled era of life sciences.

The course is designed to be delivered within about two weeks (cross your fingers), and it will feel like a sprint rather than a marathon -- but I anticipate that each one of you will learn some new things and enjoy your experience (cross remaining fingers?).

The plan is that the **first week** will offer insights into the terminal environment, common UNIX tools, effective use of shell, and shell scripting in general. Depending on how we perform during the first week, we will use the **second week** to discuss AI-assisted problem solving, and learn about Python through hands-on problem solving sessions. there will be numerous exercises and a few assignments, and in-person attendance is extremely important.

The following list offers a more detailed goals of the course:

* UNIX basics & file navigation
* File manipulation (`cat`, `less`, `cp`, `mv`, etc.)
* Searching files and working with data streams (`grep`, `find`, `awk`, etc.)
* Pipes and redirection (`|`, `>`, `>>`)
* Shell scripting basics (`for`, `while`, `if`, etc.)
* Basics of Git (`init`, `clone`, `commit`, `push`, etc.)
* Creating repositories & pushing scripts
* Application of shell scripting to real world, non-trivial challenges
* Using Large Language Models such as ChatGPT or DeepSeek for programming
* Covering Python basics: variables, data types, operators, flows, and controls
* Helping you get hands-on experience through interactive exercises
* Helping you put your imagination and learnings in use through assignments
* Gaining insights into reproducible reporting
* And getting you to participate in lots and lots of discussions

It may look like this course does not have a conventional structure and it is all over the place. You are correct, though this is intentional. I hope that this, if I may, 'modern structure' will work for most of you since it has multiple qualities:

* **Balanced Approach**. Rather than focusing on a single topic in great depth, the course structure goes in and out of fundamentals of programming, UNIX tools, Git, BASH scripting, and Python in a structured way.
* **AI Tools Introduction** – Placed at an optimal time to help you learn efficient debugging and script generation with ChatGPT or DeepSeek before diving into Python.
* **Real-World Applications** – Assignments resemble real-world problems researchers often run into during their day-to-day workflows rather than unrelatable hypothetical programming tasks.
* **Version Control & Reproducibility** – Introduction to Git and GitHub ensures that you will develop good coding habits early in your journey.
* **Quick and Comprehensive Recall**. You will get to apply everything you have learned to solve an actual problem (that literally happened and someone had to solve it in their daily work).

We will see how everything goes, and you will tell me at the end what worked and what did not :)

## Technical Recommendations

Here are a few recommendations that will help you throughout this course and beyond if you plan to apply what you have learned.

### Shell Environment

We will make quite a heavy use of the terminal environment (any terminal that gives access to a UNIX [shell](https://en.wikipedia.org/wiki/Unix_shell), like [BASH](https://en.wikipedia.org/wiki/Bash_(Unix_shell))).

If you are using Linux or Mac OSX, you have native access to a reasonable shell. Please take a moment to find out how to open your terminal now.

If you are using Mac OSX, I would strongly recommend you to install [iTerm2](https://iterm2.com/) and use it instead of the default terminal application on Mac OSX. If you are using Windows, you need to install [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install) to access to a UNIX shell.

We will spend a lot of time learning about the UNIX shell, but here are a few resources if you would like to take a brief look ahead of the course:

* [Beginner's Guide to the BASH](https://www.youtube.com/watch?v=oxuRxtrO2Ag) (a video introduction to the command line environment -- although Joe Collins is talking about Linux, the topics are relevant to anyone who uses a command line environment and Meren strongly recommends everyone to take a look).
* [Learning the Shell](http://linuxcommand.org/lc3_learning_the_shell.php) (a chapter from the open book "*The Linux Command Line*" by William Shotts -- Meren highly recommends).

### Code Editors

Throughout this course you will be writing and editing lots of code. Writing good code comfortably requires a good and comfortable text editor that is *designed* to write code, which excludes text editors such as Microsoft Word or NotePad which are more appropriate for daily writing tasks. While Meren exclusively uses [vim](https://www.youtube.com/watch?v=-txKSRn0qeA) for his everyday coding tasks, students are encouraged to consider using a **graphical code editor with syntax highlighting** that supports multiple programming languages.

You can consider [Sublime Text](https://www.sublimetext.com/), which works on Mac, [Geany](https://www.geany.org/), which works on Linux, or [Notepad++](https://notepad-plus-plus.org/) which works on Windows. These are all very nice lighthweight editors. Alternatively you can consider [VS Code](https://code.visualstudio.com/), which is a behemoth that works on Windows, Mac, and Linux (Meren doesn't like it, but promises to not judge anyone).

### Version Management System

The delivery of assignments will require you to use [Git](https://en.wikipedia.org/wiki/Git) version management system, and have an account on [GitHub](https://github.com), which is an online service built to store Git repositories in the cloud. We will discuss both git and GitHub in detail, but you should open an account on GitHub unless you already have one as soon as you read these sentences.


## Course Data Package

Hands-on exercises and assignments throughout the course will make use of a previously prepared data package. You need to download this data package on your computer, and uncompress it. You will do it only once, and then every time you open your terminal, you will go into the data pack directory since all commands will assume that you are in that directory.

Since at this point everyone have their terminals ready, all you need to do is to open a terminal, and paste the following commands in it:

```sh
## download the package on your computer using 'curl'
curl -L 'https://www.dropbox.com/scl/fi/kjmkzv35fx4vv4pqo6qhl/PFLS-DATA-PACKAGE.tar.gz?rlkey=brll5t615ubsjx2bbdvjg9gme' -o PFLS-DATA-PACKAGE.tar.gz

## unpack the archive
tar -zxvf PFLS-DATA-PACKAGE.tar.gz

## go into the data pack directory (the part that you will do over
## and over again when you open a new terminal with the intention
## to follow the course content
cd PFLS-DATA-PACKAGE
```

If you got an error in any of the lines above, please do not continue before addressing the issue.

## Final Checks

This is time for Meren to make sure that every participant,

* Has a access to a properly setup computer (with WSL for Windows users, etc),
* Has a working terminal environment that runs BASH,
* Has the data package downloaded and ready to go,
* Is able to run `git`,
* Has a GitHub account,
* Is enthusiastic about this course and put their war paint on to deal with whatever it will bring into their life.

If we are all good to go, we can start now.

## Course Content

### The Shell

The purpose of this section is to familiarize you with your terminal and some of the features of the shell environment.

At the end of this section, you will have an understanding of the power of the command line environment and how to use multiple programs in tandem to get things done interactively.

#### What is shell, exactly?

I think this is a good place to start because the shell is one of the most powerful yet underappreciated programs on your computer. In simplest terms, the shell is a command-line interface that sits between you and the operating system kernel —the hidden core of your system responsible for running program processes, executing their instructions on the CPU, managing memory, and handling storage and peripherals. At any given time, your computer runs a gazillion programs; from those that respond to your inputs through a keyboard and mouse to those that literally show where your cursor is on your screen or keeps track of all the windows you have opened, and all of them go through the kernel which allows your computer to run more programs than the number of available CPU cores, gives all programs the illusion of unlimited memory despite its physical limits, and takes care of priorities among them so all programs big and small can use the hardware resources on your computer in harmony to serve you. The shell is another program running on your computer. Like many others, it operates in what we call 'user space' and communicates with core processes running in 'kernel space' through system calls, ensuring that our commands and requests integrate smoothly into the operating system's workings without disrupting its stability. In many ways, the shell is your command center that gives you direct control over your system’s most powerful functions; in a way, it is the computer equivalent of an airplane cockpit or a nuclear reactor control room. On top of its abilities to run programs, modern shells come with scripting languages and control structures for us to run things even more efficiently.

What happens when you type something in your terminal is a fascinating story, and have many many layers. Once you press Enter, the terminal sends your input to the shell, which parses the command, expands special characters and variables, and determines whether it is a built-in command or an external program. If it’s an external program, the shell searches the relevant directories to locate it, and if the program is found and has the right permissions, the shell creates a new process by invoking the kernel, and puts the target program into the driver seat, upon which the kernel loads the program's code into memory, registers it with the scheduler, and begins execution. Throughout its runtime the program interacts with the kernel via system calls until it either terminates normally (successful or not) or the kernel forcefully kills it due to a myriad of reasons such as unhandled fatal signals, resource exhaustion, or illegal memory access. Regardless of how it ends, the kernel cleans up the kitchen and notifies the shell, which retrieves the exit status and resumes its eternal wait for your next command, ready to start the cycle all over again.

Shells are complex, and you don't need to learn any of these things, obviously. But as an undergraduate student of computer science I *did* implement a shell from scratch. My shell was neither as good nor as talented as any of the modern shells. But doing that, and forcing myself to go through that suffering created many new synapses in my brain that eventually afforded me the framework I needed to understand and solve more complex problems elsewhere. Which, in a big plot twist, brings me to the use of AI. One may argue that almost none of the things we will discuss throughout this course are necessary to learn. Indeed, you can solve the vast majority of programming challenges, including those that we will cover during our exercises, using popular LLM clients such as ChatGPT or DeepSeek almost instantaneously. But there are no cheat codes to life -- without going through the pain of truly understanding something, there is no way to achieve mastery, lead with confidence, or create something original.

<blockquote>
There is a tide in the affairs of people,<br/>
Which, taken at the flood, leads on to fortune;<br/>
Omitted, all the voyage of their life<br/>
Is bound in shallows and in miseries.

<div class="blockquote-author">William Shakespeare (Julius Caesar, Act 4, Scene 3)</div>
</blockquote>

I am asking you to not use any AI tools unless I ask you to throughout this course. Do not prevent yourself from taking advantage of the tide that is rising in your academic life due to this course. Once you have some fundamental understanding of the topics we cover here, you will be much more efficient using AI-assistance for your dealings with shell and any other coding tasks.

Now, please open your terminals and let's start with the simplest shell functions and most commonly used UNIX programs in it that help us navigate through files and directories, study their contents, and special characters that make life easier.

#### Navigating files and directories

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
                          # directory :))  
```

---

`pwd` – Print working directory

Displays the absolute path of the current directory.

```sh
pwd  # Show full path of the current directory
```

---

`mkdir` – Create a new directory

Creates a new empty directory.

```sh
mkdir my_folder       # Create a directory named my_folder
mkdir -p parent/child # Create nested directories
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

`touch` – Create an empty file

Creates a new empty file.

```sh
touch newfile.txt  # Create an empty file called newfile.txt
```

---

`find` – Search for files and directories

Finds files in a directory hierarchy.

```sh
find /home -name "document.txt"  # Search for document.txt in /home
find /var -size +100M            # Find files larger than 100MB in /var
find . -type f -name "*.log"     # Find all .log files in current directory
```

The last command uses `*`, a special character (so-called 'wildcard') that is extremely useful to target multiple files that match to a particular pattern. This is not the only special character, and most shells will process user input commands and interpret a series of special characters when they are found. Before we continue with more fun programs, let's take a look at a list of commonly used special charaters first.

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

#### Special Characters

`$` – Variable substitution and command substitution

Used to reference variables and execute commands. Extremely important character that we will use many many times.

```sh
echo $HOME       # Prints the home directory
echo $(date)     # Runs the date command and prints the result
```

---

`#` – Comment

Everything after `#` on a line is ignored by BASH. It is mostly useful when writing shell scripts.

```sh
## This is a comment
echo "Hello World"  # Prints Hello World
```

---

`*` – Wildcard (matches multiple characters)

Matches all files and directories in a given location.

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
mkdir newdir && cd newdir  # Creates and moves into newdir if successful
```

---

`||` – Logical OR (run second command if first fails)

Runs the second command **only if the first one fails**.

```sh
mkdir mydir || echo "Directory creation failed"  # Prints message if mkdir fails
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

In a nutshell, AWK scans a file line by line, splits input into fields based on a separator, enables pattern-based filtering, and allows users to perform actions on matching lines, and help produce highly formatted reports.

It is really difficult to demonstrate the utility of AWK without a few examples, so I put together the following file of all chancellors of Germany where the columns indicate (1) the name of the chancellor, (2) their education, (3) the age at which they became a chancellor of Germany, (4) the number of year they served at this position, and (5) the year they assumed this position. You can save the contents of this file as `german_chancellors.txt` in your working directory, and follow the examples below:

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

Calculate each chancellor’s age today (assuming we are still in 2025 by the time you're seeing this):

```sh
awk '{print $1, "is (or would have been)", 2025 - $5, "years old today"}' german_chancellors.txt
```

Find chancellors born after women gained voting rights in Germany:

```sh
awk '$5 > 1918 {print $1, "was born in", $5}' german_chancellors.txt
```

---

Let's say we want to print the first letter of the first name and last name of each chancellor who studied law given this file. How can we do that?

First we can filter this file for chancellors who studied law. There are multiple ways to do it. For instance, this is one way to do it:

```sh
awk '{if ($2=="Law") print }' german_chancellors.txt
```

but a more stylish way to do it in AWK would make use of patterns:

```sh
awk '/Law/{print}' german_chancellors.txt
```

where `/Law/` is a pattern and means "run the following statemnt on every line that contains the word Law".

To get the first letter of the first name and last name of each chancellor, we first need to get their names. So we can do it by printing the exact column where names appear:

```sh
awk '/Law/{print $1}' german_chancellors.txt
```

We see where the first and the last names are, but how to access them separately? Well, one way to do it is to 'split' this field into fields. A special keyword in AWK is `FS`, which means field separator. By default, the field separator is whitespace, but we can replace it with anything BEFORE running AWK commands on anything by simply mentioning it to AWK at the very beginning of the process.

Here is an example:

```sh
awk '/Law/{print $1}' german_chancellors.txt | awk 'BEGIN{FS="_"} {print $1 ", " $2}'
```

Cool. Now we can access name and lastname as separate columns. But how to get the first letter of each of these? Well, if you were to search on Google "how to get the first letter of a variable in AWK", you would learn that there is an AWK 'function' called `substr`, and we can use that to get what we want:

```sh
awk '/Law/{print $1}' german_chancellors.txt | awk 'BEGIN{FS="_"} {print(substr($1, 1, 1) substr($2,1,1))}'
```

---

OK. Here is a relatively difficult one: find the average age at which chancellors took office. This will take care of it, and let's brake it down to its individual components to discuss what we are looking at here:

```sh
awk '{sum_age += $3; count++} END { print "Average age when becoming chancellor:", sum_age/count}' german_chancellors.txt
```

### EXC-001

Let's do an exercise using a few of the things we have learned so far. For this exercise, please go into the relevant exervise directory in the course data pack:

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

Try your best, and it is OK if you can't answer each one of them. If you try your hardest, the solutions will make much more sense even if you fail.

Once you are done, we will review the [solutions](solutions/EXC-001) together.

## Shell Scripting

The last section focused on how the command line environment and the common tools that are accessible to us in the UNIX shell can empower its users to perform tasks that would have taken much longer to do manually. For instance, one could find how many genes in a given FASTA file by literally going throught it line by line in their text editor, but the ability to perform this task with a single command is a life saver. Using individual lines of instructions interactively is very powerful, but not suitable to complete repetitive tasks or implement complex ideas that require multiple instructions to be run one after another.

Let's talk more about what shell scripting is and where it comes handy through a realistic example. Consider the following.

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

You could indeed indeed use `grep` and `wc` to put together a command that gives you the number of contigs in a given genome, and run those command one by one on each one of these genomes to generate a report to make decisions later. Here how that looks like on my terminal:

```sh
meren $ grep '>' GCA_018224395.1_genomic.fna | wc -l
1
meren $ grep '>' GCA_019061405.1_genomic.fna | wc -l
30
meren $ grep '>' GCA_022836975.1_genomic.fna | wc -l
192
meren $ grep '>' GCA_023052945.1_genomic.fna | wc -l
1
meren $ grep '>' GCA_902636535.1_genomic.fna | wc -l
108
meren $ grep '>' GCF_000008385.1_genomic.fna | wc -l
1
meren $ grep '>' GCF_000306885.1_genomic.fna | wc -l
1
meren $ grep '>' GCF_000376585.1_genomic.fna | wc -l
1
(...)
```

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

## capture the output of a command, and put it into a variable:
my_birthday="May 14 2002"
seconds=$(date -d "$my_birthday" +%s)
echo "It has been about $seconds seconds since I was born!"

## you can also define a variable that stores lists of things
fasta_files=*.fna
echo $fasta_files

## the placeholder nature of variables can lead to creative
## applications! what do you think will happen when we run
## the next two lines here?
x=ls
$x
```
{% include CODEBLOCKFILENAME filename="variables-user-defined.sh" %}

The second class of variables are **environmental varaibles**: previously defined variables by various processes, including those that are set everytime you open a terminal. You can see all of these variables by simply typing the command `env` in your terminal, and you can access to any of these variables from within your BASH scripts:

```sh
echo "My username is '$USER'. There are $(ls $HOME | wc -l) files in my home folder, which is at '$HOME'."
```
{% include CODEBLOCKFILENAME filename="variables-environmental.sh" %}

The third class of variables are the **built-in variables** that are set everytime a command or script is run. If you have the following shell script,

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

OK. Going back to our exammple in shell, the items we wish to iterate over may come from another command. For instance, the following notation in BASH will give you numbers from 1 to 5 as a sequence (you can copy paste it in your terminal to see):

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

Great. Do you remember our glorius first shell script we wrote to count the number of sequences in all FASTA files in that directory?

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

OK. There is one more topic to cover, and after that we will come back to these genomes, which will be your first assignment!

### Conditionals

So far we discussed how to define and make use of variables, and how to build loops using `for` in our shell scripts. Our `for` loops run on all items, without having to make any decisions. But real-world tasks often require some sort of decision making and performing an operation only if a certain condition is, or a few of them are, met (or do other things if they don't!). Those conditions could include a variety of considerations such as doing something only if a file exists or absent, comparing variables and taking action depending on whether they are equal or not, and thus controlling the flow of our script based on our expectations from it.

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

Here the *condition* is `$number -gt 2`, and the `-gt` operator, which means '*greater than*', does the heavy lifting of the entire opretion. When a number is truly greater than 2, that *condition* becomes a true statement, and it meets the criterion for the `if` statement to continue.

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

What if we didn't want to use `!`? What is it we would need to do to get the numbers that are not greater than 2?

This is a good time to have a look at all the operators that one can use wen they are forming their conditions:

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


One can also extend `if` satements with additional conditions and `else` statements:

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

You have 34 Wolbachia genomes downloaded from the GTDB, and a text file that shows which genome matches to which host organism Wolbachia infects. Please write a shell script that does the following tasks when it is run in the `EXC-002` directory:

* Creates a directory called `WOLBACHIA-GENOMES`,
* Identifies Wolbachia genomes with less than x number of contigs, where x is sent to the shell script as a parameter,
* Creates a copy of each genome that match to the maximum number of contigs condition in `WOLBACHIA-GENOMES`,
* But uses the `wolbachia-hosts.txt` to rename each genome to match the host name from which the Wolbachia was recovered.

We will go trough the [solution](solutions/EXC-002) together once you have given this exercise your best shot.


## Working with Git

### Introduction to Version Control

Probably everyone who is going through this document is familiar with the fact that scientific writing is not a linear process with a clear endpoint. If you have ever had multiple versions of the same document on your disk named like `final.doc`, `final_final.doc`, `really_final.doc`, `final_Jan292025,doc`, `final_Jan292025_final.doc`, you already know what I am talking about with firsthand experience. We often have multiple versions of the same document like that instead of changing the same file is simple: we want to make sure that if we ever need to go back in time to recall a version of the text we are working on, those earlier versions of the text will not be overwritten with the newer versions of it and be forever lost. But of course this leads to a lot of confusion and redundancy.

The same happens when we write code, where the number of changes that occur in a single file can be very high. For instance, here is the top ten files in the anvi'o code based on how many times they were changed during the past 10 years.

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

This is precisely what version control systems do for you. They offer you means to record your changes to files over time, allowing you to keep track of everything, go back to a specific version of your file, or recall what was actually changed between two versions of it. This also enhances transparency, allowing others to track your work’s evolution and collaborate while keeping a clear record of each change.

As tracking changes through version control systems is not limited to code, it is particularly useful for scientists whose work requires them to write things, whether those things are experimental protocols, code, or text, and whether their priority is to track the evolution of their work, collaborate with others, and ensure reproducibility.

In the past few sections of this course, we focused heavily on text-based communication with your computer. You wrote commands in your terminal, put a few of them in a file to run tasks in batch, implemented comprehensive shell scripts. Now that we have explored the benefits of version control in managing changes, let's look at how it integrates with text-based workflows and cloud services.


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

Luckily, Git comes as the clear choice for most of us as a version control system as a free and general purpose version control system supported by cloud solutions and suitable for academoc use. We all should thank [Linus Torvalds](https://en.wikipedia.org/wiki/Linus_Torvalds) for it.

Knowing about Git and being able to use it is already enough to manage your files, track your modifications, and maintain an organized history of your work, not only to manage workflows or shell scripting, but also for your papers, reproducible data analyses, and so on. If you combine Git with Markdown, there is nothing you can't do! For instance, this entire course was prepared using Git and [Markdown](https://www.markdownguide.org/basic-syntax/), a simple mark-up language to create well structured documents that can be ported into anything -- from web pages to PDF documents to Microsoft Word files using a tool like [pandoc](https://pandoc.org/) in the comfort of your terminal environment :)


### Introduction to GitHub

Git is an open-source, stand-alone program that enables you to create local repositories for version control -- local repositories you store on your own computer or the University servers you can connect with your terminal.

GitHub, on the other hand, is a platform that enables you to store your Git repositories on the 'cloud.' It is a proprietary and for-profit developer platform that has over 1 billion dollars in revenue while itself is not open-source. Unlike Git, which is one of us, GitHub is one of *them*, if you will. But even though Linus Torvalds [hates](https://news.ycombinator.com/item?id=36123124) GitHub for various technical reasons, and I personally hate it for the gargantuan enterprise it represents, GitHub has successfully captured the attention of millions of developers worldwide as it made software development more accessible, manageable, interactive, and fun by adding key features to their cloud-based hosting service such as (1) graphical user interfaces for repository management and code reviews, (2) a service for 'pull requests' that help team members to [review and discuss code changes](https://github.com/merenlab/anvio/pull/2155) before merging them, (3) providing a platform to [report issues and discuss them](https://github.com/merenlab/anvio/issues/1248), (4) [action and workflow](https://github.com/merenlab/anvio/actions) support for continuous integration or testing, and many, many more. GitHub is free unless you wish to have private repositories, in which case you are asked to pay a fee (which is about 4 Euros as of 2025).

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
## will kep all your Git repositoies. I usually like to
## keep them in a directory called `github` under my home
## directory, so let's create one for you:
mkdir -p ~/github

## and enter into that directory:
cd ~/github

## now create a new project directory, which will keep
## all our files that will be a part of this projet.
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

You will use this command over and over again. Just so you know, Git keeps system-wide settings in another hidden file in your home directory at `~/.gitconfig`. You can edit this file to add aliases for some commands. Here is how mine looks like jost for your reference:

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

Except where the blanks are, of course, which should be replaced with actual numbers your script will calculate for a given FASTA file. The GC-content

Once you are done with your script, please commit it to a GitHub repository called `PFLS`, and put it in a directory called `EXC-003`, with the file name `fasta-file-processor.sh`.

This is extremely important, because I will test your solutions by cloning your repository on my computer the following way,

```sh
git clone https://github.com/_______/PFLS.git
```

and testing it on a FASTA file of my choosing using a command like this:

```sh
bash PFLS/EXC-003/fasta-file-processor.sh test.fa
```

You can use this information to make sure things will work on my end (if `cd ~/github && ls PFLS/EXC-003/fasta-file-processor.sh` does not produce an error, you're good). You can test your program using any FASTA file, incuding those that we have in the data package, but you do not know the *FASTA files* I will be using to test your script on my end.

If `git clone` step and the step of running your script on a FASTA file with works with the expected out works without a problem, you have a full grade for this exercise, *even if the numbers are not correct*.

---

A small tip! You can use the AWK function `gsub` to calculate the G and C bases in a given sequence:

```sh
echo "ATCGATCGCG" | awk '{gc_count += gsub(/[GgCc]/, "", $1)} END {print gc_count}'
```

But of course what you want is GC-content, which essentiall is the percentage of GC bases to _all_ bases.

The solution is [here](solutions/EXC-003), and usual, we will go through it together once you are done.

Good luck!

### EXC-004

We are rapidly moving to another exercise to immerse you into the joy of writing shell scripts for complex data tasks.

This is the most realistic problem we are going to be working on so far, and Sarahi Garcia has all the background on this problem that we had to resolve last year for a study. She can help us better undertand its context, but here I will provide a very brief description of what is going on, and what needs to be done.

Please go into the `EXC-004` directory in your data pack, and take a look at the contents of it. You will realize that

* There are multiple directories in `RAW-DATA` directory that goes like `DNA57`, `DNA58`, `DNA64`, etc. Each of these directories contains the results of a genome-resolved analysis of the metagenomic sequencing of a single culture generated from Lake Erken samples. Since each culture was started with just a few cells from the environment, the researches who conducted this study was able to recover one or more metagenomic bins. The directory names match to library names rather than sample names, and the actual names of these cultures are stored in the file `sample-translation.txt`. Each of these directories have identical structure, and contain the following two files and another directory:
  - `bins/` -- Contains the FASTA files that emerged from the automatic binning of each sample.
  - `checkm.txt` -- Completion and redundancy estimates for the tentative genomes represented by individual FASTA files in the `bins/` directory based on bacterial and archaeal single-copy core genes.
  - `gtdb.gtdbtk.tax` -- Taxonomy for each of them based on Genome Taxonomy Database.

One of the most critical next step is to estimate the actual abundances of individual genomes. But this particular organization of these data is not very useful for such downstream analyses, and we need a shell script that can help us combine these data into a more meaningful representation. The script needs to do the following things:

* Generate a new directory at the same level of `RAW-DATA` called `COMBINED-DATA`.
* Process each FASTA file in individual directories to (1) copy each `bin-unbinned.fasta` into the to `COMBINED-DATA` directory as `XXX_UNBINNED.fa` and (2) copy every other FASTA file into the to `COMBINED-DATA` directory as `XXX_YYY_ZZZ.fa` where,
  * `XXX` is the culture name recovered from `sample-translation.txt`,
  * `YYY` is `MAG` if the completion is 50% or more and contamination is 5% less according to the information in stored in the relevant `checkm.txt` file, otherwise `YYY` is `BIN`,
  * `ZZZ` is `001`, `002`, and so on, where each `MAG` and each `BIN` for individual cultures have sequential numbering.
* Ensure that each sequence in each FASTA file has a unique defline associated with the culture (`XXX`) -- pro tip: you can use {% include PROGRAM name="anvi-script-reformat-fasta" %} for this.
* Copy `checkm.txt` and `gtdb.gtdbtk.tax` files in individual `bins/` directories into `COMBINED-DATA` as `XXX-CHECKM.txt` and `XXX-GTDB-TAX.txt`.

Based on these instructions, when you run your script in the directory `EXC-004`,

```sh
bash generate-combined-data.sh
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

Once you are done, please commit your script to your GitHub repository for `PFLS`, and put it in a directory called `EXC-004`, with the file name `generate-combined-data.sh`.

## Course Responsibilities

The evaluation in this course will be based on **five parts** of a portfolio that we can divide into **two major components**:

* The first major component is '**class citizenship**' emails, described below in the section **Part I** (20% of your grade).
* The second major component is '**programming exercises**', described below in sections for **Part II, III, IV, and V** (80% of your grade).


### Part I (20% of your grade) – Class Citizenship

Class citizenship emails **will track your attendance and engagement to the course** and will help the course director to have an overall understanding of the evolution of the course.

The class citizenship demands every participant to send a class citizenship email at the end of each day to _meren@hifmb.de_ **and** _sarahi.garcia@uol.de_ (*and* as in, sending the same email to both Meren and Sarahi).

The class citizenship email must be composed of two parts:

1. A brief summary of the main concepts discussed during the day, interpreted by the attendee in their own words.
2. A short question that is relevant to a concept or idea discussed during the lecture.

The last 10 minutes of every course day will be dedicated to writing the class citizenship emails, therefore the attendees will leave the class without having to remember doing it later. **The class citizenship emails that are sent after the end of the class will not be taken into consideration as a mark of attendance**.

The title of the class citizenship email **must follow this pattern word-by-word** where you will need to replace DD/MM/YY with the date, month, and the year (despite the simplicity of this request many students have failed to follow these instructions, so you are our last hope):

> PFLS Citizenship: DD/MM/YY

The best class citizenship emails are those that are brief, genuine, and insightful. In an ideal world the emails should be no less than 50 words, and no more than 150 words. Please do not send notes you take throughout the class -- your notes are for you, not for us. You should use the last 10 minutes of the day to gather your thoughts, and come up with a summary of what you can remember.

Here is an example class citizenship email:

> Summary: Today we discussed what is phylogenomics, how phylogenomic trees are built, and why single-copy core genes are suitable for building phylogenomics trees. We also discussed the relationship between phylogenetics, phylogenomics, and pangenomics with respect to the fraction of
genome used and the evolutionary distance that they can cover.
>
> Question: Since phylogenomics and pangenomics are both useful for inferring evolutionary distances, it seems to me that integrating both methods in a systematic way would yield a more reliable tree. But it looks like the field only uses phylogenomics and pangenomics separately, is there a reason for that?

### Part II, III, and IV (30% of your grade)

This part will be composed of three mini programming exercises that you will have to implement and return. Each programming exercise will provide you with explicit instructions regarding the nature of the data and question, and what the program is expected to achieve. You will use your learnings in the course to complete the programming tasks and submit the resulting source code.

### Part V (50% of your grade)

This part will be composed of a single large programming exercise that will require you to orchestrate multiple programming tools and languages.

## Grading

Grading scale:

|**Grade**|**Threshold**|
|:--|:--|
|1.0|95%|
|1.3|90%|
|1.7|85%|
|2.0|80%|
|2.3|75%|
|2.7|70%|
|3.0|65%|
|3.3|60%|
|3.7|55%|
|4.0|50%|

For all grading related questions, please consult with Sarahi Garcia.

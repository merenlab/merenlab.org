#!/bin/bash

# the purpose of this script is to get information from our internal
# github repository to under merenlab.org

set -e

web_repo_dir=~/github/merenlab.org
courses_repo_dir=~/github/courses

process () {
    # send this function three parameters in the following order:
    #
    #   - the target output directory for contents to be created
    #   - the markdown file for the course syllabus
    #   - the PDF output for the course syllabus
    #
    #
    if [ ! -d $1 ]; then
        mkdir -p $1
    fi

    echo '---' > $1/index.md
    echo 'layout: page' >> $1/index.md
    echo 'noleftpanel: true' >> $1/index.md

    cat $2 | \
        sed 1,1d | \
        sed 's/^#### /##### /g' | \
        sed 's/^### /#### /g' | \
        sed 's/^## /### /g' | \
        sed 's/^# /## /g' | \
        grep -v pagebreak >> $1/index.md

    if [ -n "$3" ]; then
        cp $3 $1/syllabus.pdf
    fi
}

########################################################################################################################
# Introduction to Popular `Omics Strategies
#######################################################################################################################
course_token="courses/ITPOS"
course_input_dir="$courses_repo_dir/introduction-to-popular-omics-strategies/00_SYLLABUS"
course_web_dir="$web_repo_dir/$course_token"

echo "Rendering $course_token ..."

process $course_web_dir \
        $course_input_dir/introduction-to-popular-omics-strategies-syllabus.md \
        $course_input_dir/introduction-to-popular-omics-strategies-syllabus.pdf

######################################################################################################################
# Applied Microbial `Omics
#######################################################################################################################
course_token="courses/AMO"
course_input_dir="$courses_repo_dir/applied-microbial-omics/00_SYLLABUS"
course_web_dir="$web_repo_dir/$course_token"

echo "Rendering $course_token ..."

process $course_web_dir \
        $course_input_dir/applied-microbial-omics.md \
        $course_input_dir/applied-microbial-omics.pdf

#######################################################################################################################
# Ecology of Marine Microbes
#######################################################################################################################
course_token="courses/EMM"
course_input_dir="$courses_repo_dir/ecology-of-marine-microbes/00_SYLLABUS"
course_web_dir="$web_repo_dir/$course_token"

echo "Rendering $course_token ..."

process $course_web_dir \
        $course_input_dir/ecology-of-marine-microbes.md \
        $course_input_dir/ecology-of-marine-microbes.pdf

#######################################################################################################################
# Programming for Life Scientists
#######################################################################################################################
course_token="courses/PFLS"
course_input_dir="$courses_repo_dir/programming-for-life-scientists/00_SYLLABUS_AND_CONTENT"
course_web_dir="$web_repo_dir/$course_token"

echo "Rendering $course_token ..."

# Process the index
process $course_web_dir \
        $course_input_dir/programming-for-life-scientists.md \
        $course_input_dir/programming-for-life-scientists.pdf

# NEXT, we process exercises by reading clear text, and encryptig them on the fly.
# ~/github/courses/programming-for-life-scientists/00_SYLLABUS_AND_CONTENT/EXERCISES/README.md
# explains how this blackmagic works.

exc_paths=`ls -d ~/github/courses/programming-for-life-scientists/00_SYLLABUS_AND_CONTENT/EXERCISES/EXC*`

# This is the FIRST LOOP to get all exercise solutions to be processed
# and rendered. Since the encryption has to wait for the rendered files,
# and Jekyll will take its sweet time to do it, we will update everything
# to kick-start the rendering, and then wait for rendered files in a
# second loop to encrypt them.
for exc_path in $exc_paths
do
    exc_name=`basename $exc_path`
    exc_solution_dir="$course_web_dir/solutions/$exc_name"
    rendered_output="$web_repo_dir/_site/$course_token/solutions/$exc_name/index.html"

    mkdir -p $exc_solution_dir

    # this is a bit cray. here we will remove the rendered output,
    # process the non-rendered origincal index.md from the courses
    # repository for the exercise into its place, and this way
    # make sure when there is a file at $rendered_output, it is
    # the newly rendered file ready for encryption
    rm -rf $rendered_output
    process $exc_solution_dir \
            $exc_path/index.md
done

# This is the SECOND LOOP to wait all exercises to be rendered.
echo -n "  - waiting for Jekyll to render all exercise files "
for exc_path in $exc_paths
do
    exc_name=`basename $exc_path`
    rendered_output="$web_repo_dir/_site/$course_token/solutions/$exc_name/index.html"

    while [ ! -e $rendered_output ]
    do
        echo -n "."
        sleep 1
    done
    echo " done! :)"
done

# This is the THIRD LOOP to encrypt the resulting files.
echo -n "  - encrypting files for "
for exc_path in $exc_paths
do
    exc_name=`basename $exc_path`
    exc_name_callback=`echo $exc_name | awk '{print tolower($0)}'`

    echo -n "$exc_name .. "

    exc_solution_dir="$course_web_dir/solutions/$exc_name"
    rendered_output_dir="$web_repo_dir/_site/$course_token/solutions/$exc_name"

    # what follows is some sort of a cray hack. first, learn the password the instructor
    # set in the courses repo.
    exc_password=`cat $exc_path/password.txt`

    # go into the rendered output directory, encrypt the index.html,
    # and go back to where you came from
    cd $rendered_output_dir
    staticrypt index.html \
        -p $exc_password --short \
        -t $web_repo_dir/_misc/staticrypt-template.html \
        --template-title "Solution for $exc_name" \
        --template-instructions "You shall not pass! Well, unless you have the password. Then you pass .. or <a href=\"/$course_token/#$exc_name_callback\">go back to where you came from</a>."
    cd - > /dev/null 2>&1

    # now remove the non-rendered index.md from the courses output directory
    rm -rf $rendered_output_dir/index.md

    # and replace it with the encrypted index.html so that in the web page that will be
    # committed to GitHub actually has the password protected pages.
    mv $rendered_output_dir/encrypted/index.html $exc_solution_dir/index.html
done
echo "done! :)"

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




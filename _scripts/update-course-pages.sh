#!/bin/bash

# the purpose of this script is to get information from our internal
# github repository to under merenlab.org

set -e

process () {
    # send this function three parameters in the following order:
    #
    #   - the target output directory for contents to be created
    #   - the markdown file for the course syllabus
    #   - the PDF output for the course syllabus
    #

    echo '---' > $1/index.md
    echo 'layout: page' >> $1/index.md
    echo 'noleftpanel: true' >> $1/index.md

    cat $2 | \
        sed 1,1d | \
        sed 's/^### /#### /g' | \
        sed 's/^## /### /g' | \
        sed 's/^# /## /g' | \
        grep -v pagebreak >> $1/index.md

    cp $3 $1/syllabus.pdf
}

# Introduction to Popular `Omics Strategies
process ~/github/merenlab.org/courses/ITPOS \
        ~/github/courses/introduction-to-popular-omics-strategies/00_SYLLABUS/introduction-to-popular-omics-strategies-syllabus.md \
        ~/github/courses/introduction-to-popular-omics-strategies/00_SYLLABUS/introduction-to-popular-omics-strategies-syllabus.pdf

# Applied Microbial `Omics
process ~/github/merenlab.org/courses/AMO \
        ~/github/courses/applied-microbial-omics/00_SYLLABUS/applied-microbial-omics.md \
        ~/github/courses/applied-microbial-omics/00_SYLLABUS/applied-microbial-omics.pdf

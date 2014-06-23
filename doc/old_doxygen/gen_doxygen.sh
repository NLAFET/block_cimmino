#!/usr/bin/env bash

function watch_dir {
    cd $1
    dir=$PWD
    while inotifywait -e modify ./; do
        echo "$1 has changed, working..."
        doxygen &>> /tmp/log_doxygen_abcd
        #cp latex_stuffs/* ../doc_output/doxy/latex
        #cd ../doc_output/doxy/latex
        #make
        #cd $dir        
    done
}

watch_dir ../src &
watch_dir ../include &
# watch_dir ../doc


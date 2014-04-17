#!/usr/bin/bash

function watch_dir {
    cd $1
    while inotifywait -e modify ./; do
        echo "$1 has changed, working..."
        doxygen &>> /tmp/log_doxygen_abcd
    done
}

watch_dir ../doc &
watch_dir ../src &
watch_dir ../include &


#!/usr/bin/env bash

function watch_dir {
    while inotifywait -e modify $1; do
        echo "$1 has changed, working..."
        doxygen 
        # ./gen_sphx.sh pdf
        ./gen_sphx.sh
    done
}

# watch_dir ../src &
watch_dir ./source &
watch_dir ../include


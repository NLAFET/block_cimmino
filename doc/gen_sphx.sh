#!/usr/bin/env bash

export SB="/usr/bin/env sphinx-build"

function do_job {
  $SB -b html -E -j2 source build/html
  $SB -b latex -E -j2 source build/latex
  # cd build/latex
  # make
  # cd -
}

do_job

while inotifywait -e modify ./source; do
    echo "has changed, working..."
    do_job
done

#!/usr/bin/env bash

export SB="/usr/bin/env sphinx-build"

function do_job {
  $SB -b html -E -j2 source build/html
  $SB -b latex -E -j2 source build/latex
  cd build/latex
  sed -n '
  1h
  1!H
  $ {
    g
    s/\\begin{fulllineitems}\n*\\end{fulllineitems}//g
    p
  }' ABCDSolver.tex > temp.tex
  rm ABCDSolver.tex
  mv temp.tex ABCDSolver.tex
  make clean
  echo "Building pdf file"
  make &> /dev/null
  cd -
}

do_job $1

if [ "$1" == "update" ]
      then
        echo "Uploading to abcd.enseeiht.fr"
        rsync -raz --progress --delete \
          build/latex/ABCDSolver.pdf \
          build/html/* \
          abcd.enseeiht.fr:/var/www/abcd/
fi

# while inotifywait -e modify ./source; do
    # echo "has changed, working..."
    # do_job $1
# done

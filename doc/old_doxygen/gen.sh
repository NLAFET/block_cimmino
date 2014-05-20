#!/usr/bin/sh

function comp(){
    pandoc -s --mathjax -V geometry:margin=1in\
        -f markdown+footnotes+definition_lists+example_lists+raw_tex\
        --highlight-style=tango -N\
        -c linear/docco.css --toc\
        -o ../doc_output/abcd_usage.pdf usage.md

    pandoc -s --mathjax \
        -f markdown+footnotes+definition_lists+example_lists+raw_tex\
        --highlight-style=tango\
        -c linear/docco.css --toc\
        -o ../doc_output/abcd_usage.html usage.md
}

function syncToHyde(){
    dir=$PWD
    cp -r ../doc_output/abcd_usage.pdf ../doc_output/abcd_usage.html linear ${HOME}/Programming/html/perso_n7/content/abcd/
    cd ${HOME}/Programming/html/perso_n7/
    hyde gen
    rsync -rz --delete deploy/* zenadi@s.zeapo.com:
    cd $dir
}

comp

if [ "$1" == "once" ]; then
    exit
else
    if [ "$1" == "sync" ]; then
        while inotifywait -e modify usage.md; do
            comp
            syncToHyde
        done
    else
        while inotifywait -e modify usage.md; do
            comp
        done
    fi
fi


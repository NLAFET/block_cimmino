#!/usr/bin/sh

function comp(){
    pandoc -s --mathjax -V geometry:margin=1in -f markdown+footnotes+definition_lists+example_lists\
        --highlight-style=tango\
        -c linear/docco.css --toc\
        -o abcd_usage.pdf usage.md

    pandoc -s --mathjax -f markdown+footnotes+definition_lists+example_lists\
        --highlight-style=tango\
        -c linear/docco.css --toc\
        -o abcd_usage.html usage.md
}

comp

while inotifywait -e modify usage.md; do
    comp
done

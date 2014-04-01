#!/usr/bin/sh


while inotifywait -e modify usage.md; do
    pandoc -s --mathjax -V geometry:margin=1in -f markdown+footnotes+definition_lists+example_lists\
        --highlight-style=tango\
        -o abcd_usage.pdf usage.md
    pandoc -s --mathjax -f markdown+footnotes+definition_lists+example_lists\
        --highlight-style=tango\
        -o abcd_usage.html usage.md
done

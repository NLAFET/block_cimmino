(TeX-add-style-hook
 "header"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "twoside")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("xcolor" "table") ("hyperref" "pdfstartview=FitH")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "calc"
    "doxygen"
    "makeidx"
    "textcomp"
    "xcolor"
    "hyperref"
    "amsmath"
    "tikz"
    "pgfplots")))


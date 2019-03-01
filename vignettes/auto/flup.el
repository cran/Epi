(TeX-add-style-hook
 "flup"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("report" "a4paper" "dvipsnames" "twoside" "12pt")))
   (TeX-run-style-hooks
    "latex2e"
    "topreport"
    "report"
    "rep12")
   (TeX-add-symbols
    "Title"
    "Tit"
    "Version"
    "Dates"
    "Where"
    "Homepage"
    "Faculty")
   (LaTeX-add-labels
    "fig:fu2"
    "fig:Lexis-diagram"
    "fig:Ins-noIns")))


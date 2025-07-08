rem vignettes appear in alphabetical order of filename
rem this script moves .R and .pdf file to inst.doc
rem => no vignette links on CRAN pckage website 

call rw 01flup
call rt 01flup
call bl 01flup
move 01flup.R   ..\inst\doc\
move 01flup.pdf ..\inst\doc\

call rt 02addLexis
call rw 02addLexis
call bl 02addLexis
move 02addLexis.R   ..\inst\doc\
move 02addLexis.pdf ..\inst\doc\

call rt 03crisk
call rw 03crisk
call bl 03crisk
move 03crisk.R   ..\inst\doc\
move 03crisk.pdf ..\inst\doc\

call rt 04simLexis
call rw 04simLexis
call bl 04simLexis
move 04simLexis.R   ..\inst\doc\
move 04simLexis.pdf ..\inst\doc\

call rt 05yll
call rw 05yll
call bl 05yll
move 05yll.R   ..\inst\doc\
move 05yll.pdf ..\inst\doc\

call klean
rem del *.pdf

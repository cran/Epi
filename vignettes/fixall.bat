rem vignettes appear in alphabetical order of filename
rem this script moves .R and .pdf file to inst.doc
rem => no vignette links on CRAN package website 

call rw 01flup
call rt 01flup
call bl 01flup
copy 01flup.R   ..\inst\doc\
copy 01flup.pdf ..\inst\doc\

call rt 02addLexis
call rw 02addLexis
call bl 02addLexis
copy 02addLexis.R   ..\inst\doc\
copy 02addLexis.pdf ..\inst\doc\

call rt 03crisk
call rw 03crisk
call bl 03crisk
copy 03crisk.R   ..\inst\doc\
copy 03crisk.pdf ..\inst\doc\

call rt 04simLexis
call rw 04simLexis
call bl 04simLexis
copy 04simLexis.R   ..\inst\doc\
copy 04simLexis.pdf ..\inst\doc\

call rt 05yll
call rw 05yll
call bl 05yll
copy 05yll.R   ..\inst\doc\
copy 05yll.pdf ..\inst\doc\

call klean
rem del *.pdf

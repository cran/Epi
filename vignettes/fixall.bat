rem vignettes appear in alphabetical order of filename
rem this script moves .R and .pdf file to inst.doc
rem => no vignette links on CRAN pckage website 

call rw aaflup
call rt aaflup
call bl aaflup
move aaflup.R   ..\inst\doc\
move aaflup.pdf ..\inst\doc\

call rt addLexis
call rw addLexis
call bl addLexis
move addLexis.R   ..\inst\doc\
move addLexis.pdf ..\inst\doc\

call rt crisk
call rw crisk
call bl crisk
move crisk.R   ..\inst\doc\
move crisk.pdf ..\inst\doc\

call rt simLexis
call rw simLexis
call bl simLexis
move simLexis.R   ..\inst\doc\
move simLexis.pdf ..\inst\doc\

call rt yll
call rw yll
call bl yll
move yll.R   ..\inst\doc\
move yll.pdf ..\inst\doc\

call klean
rem del *.pdf

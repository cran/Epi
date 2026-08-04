rem vignettes appear in alphabetical order of filename
rem this script moves .R and .pdf file to inst.doc
rem => no vignette links on CRAN package website 

call rw n01-flup
call rt n01-flup
call bl n01-flup
rem copy n01-flup.R   ..\inst\doc\
rem copy n01-flup.pdf ..\inst\doc\

call rt n02-addLexis
call rw n02-addLexis
call bl n02-addLexis
rem copy n02-addLexis.R   ..\inst\doc\
rem copy n02-addLexis.pdf ..\inst\doc\

call rt n03-crisk
call rw n03-crisk
call bl n03-crisk
rem copy n03-crisk.R   ..\inst\doc\
rem copy n03-crisk.pdf ..\inst\doc\

call rt n04-simLexis
call rw n04-simLexis
call bl n04-simLexis
rem copy n04-simLexis.R   ..\inst\doc\
rem copy n04-simLexis.pdf ..\inst\doc\

call rt n05-yll
call rw n05-yll
call bl n05-yll
rem copy n05-yll.R   ..\inst\doc\
rem copy n05-yll.pdf ..\inst\doc\

rem move to website
call 2gp n01-flup.R   Epi
call 2gp n01-flup.pdf Epi
call 2gp n02-addLexis.R   Epi
call 2gp n02-addLexis.pdf Epi
call 2gp n03-crisk.R   Epi
call 2gp n03-crisk.pdf Epi
call 2gp n04-simLexis.R   Epi
call 2gp n04-simLexis.pdf Epi
call 2gp n05-yll.R   Epi
call 2gp n05-yll.pdf Epi

call klean
rem del *.pdf

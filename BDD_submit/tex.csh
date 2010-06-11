#!/bin/tcsh
echo "removing old versions of printable copies"
rm paper.ps
rm paper.pdf
echo "Running latex of paper"
latex paper
latex paper
latex paper
bibtex paper
latex paper
echo "Collating errors into errors.log"
unalias grep
latex paper |grep "LaTeX Warning" >errors.log
echo "creating postcript of paper"
dvips paper
echo "creating PDF of paper"
ps2pdf paper.ps 
echo "---------------------------------------------------------"
echo " Current word count is:"
wc -w paper.tex
echo "---------------------------------------------------------"

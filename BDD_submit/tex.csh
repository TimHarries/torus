#!/bin/tcsh
echo "removing old versions of printable copies"
rm paper.ps
rm paper.pdf
echo "Running latex of paper"
pdflatex paper
pdflatex paper
pdflatex paper
bibtex paper
pdflatex paper
echo "Collating errors into errors.log"
unalias grep
pdflatex paper |grep "LaTeX Warning" >errors.log
echo "---------------------------------------------------------"
echo " Current word count is:"
wc -w paper.tex
echo "---------------------------------------------------------"

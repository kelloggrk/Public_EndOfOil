#!/bin/sh

#-------------------------------------------------------------------------------
# Name:        EndOfOil_paper_script.sh
# Purpose:     Compiles the paper. LaTex code inputs figures, tables, and single
#	       number tex files from subfolders in EndOfOil/output
#
# Created:     8 November, 2024
#-------------------------------------------------------------------------------

# Define path to paper subfolder based on $REPODIR variable
# exported from EndOfOil_bash_file.sh
CODEDIR=$REPODIR/paper


# COMPILE PAPER
cd $CODEDIR
pdflatex -output-directory=$CODEDIR EndOfOil.tex
bibtex EndOfOil
pdflatex -output-directory=$CODEDIR EndOfOil.tex
pdflatex -output-directory=$CODEDIR EndOfOil.tex

# COMPILE APPENDIX
pdflatex -output-directory=$CODEDIR EndOfOil_appx.tex
bibtex EndOfOil_appx
pdflatex -output-directory=$CODEDIR EndOfOil_appx.tex
pdflatex -output-directory=$CODEDIR EndOfOil_appx.tex

# GET CROSS REFERENCES BETWEEN PAPER AND APPX
pdflatex -output-directory=$CODEDIR EndOfOil.tex
pdflatex -output-directory=$CODEDIR EndOfOil_appx.tex

#clean up log files
rm *.log

exit

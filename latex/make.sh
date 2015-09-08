#!/bin/bash
#

rm -f *~
latex manuscript
bibtex manuscript
latex manuscript
latex manuscript
dvips -o manuscript.ps manuscript.dvi
ps2pdf14 manuscript.ps manuscript.pdf

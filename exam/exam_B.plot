#!/usr/bin/env gnuplot

# set the output as png
set term pdfcairo
# Input file contains comma-separated values fields
set output "output_data/B.pdf"
set title "B plot"
set xlabel "rho"
set ylabel "B"

set grid
plot "output_data/B.dat" using 1:2:3 with yerrorbars notitle
# write to file
set output

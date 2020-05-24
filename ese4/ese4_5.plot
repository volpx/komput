#!/usr/bin/env gnuplot
# set the output as png
set term pdfcairo
# Input file contains comma-separated values fields
set output "output_data/cvv.pdf"
set title "CVV"
set xlabel "Time"
set ylabel "cvv"
# set log y
set xrange [0:1]
# set yrange [-100:400]
set grid
plot "output_data/cvv.dat" using 1:2 with lines notitle
# write to file
set output

#!/usr/bin/env gnuplot

# set the output as png
set term pdfcairo
# Input file contains comma-separated values fields
set output "data/Q.pdf"
set title "Q fluctuations T=1 rho=0.75"
set xlabel "Time"
set ylabel "Q"
# set log y
# set xrange [4:4.3]
# set yrange [-100:400]
set grid
plot "data/data_evolution.dat" using 2:6 with lines notitle
# write to file
set output

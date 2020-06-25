#!/usr/bin/env gnuplot

# set the output as png
set term pdfcairo
# Input file contains comma-separated values fields
set output "output_data/Q.pdf"
set title "Q fluctuations"
set xlabel "Time"
set ylabel "Q"
# set log y
# set xrange [4:4.3]
# set yrange [-100:400]
set grid
plot "output_data/data_evolution.dat" using 2:6 with lines notitle
# write to file
set output

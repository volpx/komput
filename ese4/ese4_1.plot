#!/usr/bin/env gnuplot
# set the output as png
set term pngcairo size 1024,768
# Input file contains comma-separated values fields
set output "output_data/energy.png"
set title "Energy fluctuations"
set xlabel "Time"
set ylabel "Energy"
set log y
# set xrange [10:87]
# set yrange [-0:1.01]
set grid
plot "output_data/energy.dat" using 1:2 with lines notitle
# write to file
set output

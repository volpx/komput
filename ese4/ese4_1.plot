#!/usr/bin/env gnuplot
# set the output as png
set term pdfcairo
# Input file contains comma-separated values fields
set output "output_data/energy.pdf"
set title "Energy fluctuations"
set xlabel "Time"
set ylabel "Energy"
# set log y
# set xrange [10:87]
# set yrange [-0:1.01]
set grid
plot "output_data/energy.dat" using 2:5 with lines title "E=T+V"\
 , "output_data/energy.dat" using 2:3 with lines title "T"\
 ,"output_data/energy.dat" using 2:4 with lines title "V"
# write to file
set output

#!/usr/bin/env gnuplot
# set the output as png
set term pngcairo size 1024,768
# Input file contains comma-separated values fields
set output "output_data/pos.png"
set title "Positions"
set xlabel "X"
set ylabel "Y"
set zlabel "Z"
# set range [0:12.7]
set ticslevel 0
splot "output_data/positions.dat" using 2:3:4 with points pointtype 7
# write to file
set output

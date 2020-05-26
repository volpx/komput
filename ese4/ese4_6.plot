#!/usr/bin/env gnuplot
# set the output as png
set term pdfcairo
# Input file contains comma-separated values fields
set output "output_data/cvv_comp.pdf"
set title "CVV comparison"
set xlabel "Time"
set ylabel "cvv"
# set log y
# set xrange [0:1]
# set yrange [0:1]
set grid
psn=0
plot "output_data/rho0.06E1e5/cvv.dat" using 1:2 with linespoints ps psn title "rho 0.06"\
,"output_data/rho0.2/cvv.dat" using 1:2 with linespoints ps psn title "rho 0.2"\
,"output_data/rho0.55/cvv.dat" using 1:2 with linespoints ps psn title "rho 0.55"\
,"output_data/rho0.85/cvv.dat" using 1:2 with linespoints ps psn title "rho 0.85"\
,"output_data/rho1/cvv.dat" using 1:2 with linespoints ps psn title "rho 1"
# write to file
set output

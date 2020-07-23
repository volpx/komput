#!/usr/bin/env gnuplot

# set the output as png
set term pdfcairo
# Input file contains comma-separated values fields
set output "data/Qac.pdf"
set title "Q autocorrelation"
set xlabel "Time"
set ylabel "Qac"
# set log y
set xrange [0:0.2]
# set yrange [-100:400]
l0=1
l1=1/(7.168*0.004)
l2=-0.0
set grid
plot "data/Qautocorr.dat" using 2:3 with lines notitle\
, l0*exp(-l1*x)+l2 with lines notitle
# plot "output_data/Qautocorr.dat" using 2:3 with lines notitle
# write to file
set output

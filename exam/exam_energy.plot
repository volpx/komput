#!/usr/bin/env gnuplot

# set the output as png
set term pdfcairo
# Input file contains comma-separated values fields
set output "data/energy.pdf"
set title "Energy fluctuations T=1 rho=0.75"
set xlabel "Time"
set ylabel "Energy"
# set log y
# set xrange [4:4.3]
# set yrange [-100:400]
set grid
plot "data/data_evolution.dat" using 2:5 with lines title "E=K+V+Nose"\
   , "data/data_evolution.dat" using 2:3 with lines title "K"\
   , "data/data_evolution.dat" using 2:4 with lines title "V"\
, "data/data_evolution.dat" using 2:(column(3)+column(4)) with lines title "K+V"
# write to file
set output

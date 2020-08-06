#!/usr/bin/env gnuplot

# set the output as png
set term pdfcairo
# Input file contains comma-separated values fields
set output "data/temp.pdf"
set title "Temp fluctuations"
set xlabel "Time"
set ylabel "Temp,s,vs"
# set log y
# set xrange [0:1]
# set yrange [-100:400]
set grid
plot "data/data_evolution.dat" using 2:7 with lines title "Temp"\
   , "data/data_evolution.dat" using 2:8 with lines title "s"\
    , "data/data_evolution.dat" using 2:(column(9)/10) with lines title "vs"
# write to file
set output

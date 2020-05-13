#!/usr/bin/env gnuplot
# set the output as png
set term pngcairo size 1024*2,768*2
# Input file contains comma-separated values fields
set output "output_data/pos.png"
set title "Positions"
set xlabel "X"
set ylabel "Y"
set zlabel "Z"
set xrange [0:12.7]
set yrange [0:12.7]
set zrange [0:12.7]

set ticslevel 0
set grid
set style line 8 \
    linecolor rgb '#0060ad' \
    linetype 0 linewidth 2 \
    pointtype 1 pointsize 1.5
set style line 1 \
    linecolor rgb '#dd181f' \
    linetype 0 linewidth 2 \
    pointtype 1 pointsize 1.5
set style line 2 \
    linecolor rgb '#29ff29' \
    linetype 0 linewidth 2 \
    pointtype 1 pointsize 1.5
set style line 3 \
    linecolor rgb '#ff61ea' \
    linetype 0 linewidth 2 \
    pointtype 1 pointsize 1.5
set style line 4 \
    linecolor rgb '#fff693' \
    linetype 0 linewidth 2 \
    pointtype 1 pointsize 1.5
set style line 5 \
    linecolor rgb '#01ffb7' \
    linetype 0 linewidth 2 \
    pointtype 1 pointsize 1.5
set style line 6 \
    linecolor rgb '#738dff' \
    linetype 0 linewidth 2 \
    pointtype 1 pointsize 1.5
set style line 7 \
    linecolor rgb '#ff5500' \
    linetype 0 linewidth 2 \
    pointtype 1 pointsize 1.5

splot "output_data/positions0.dat" using 2:3:4 with linespoints linestyle 8
# write to file
set output

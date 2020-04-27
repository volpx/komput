#!/usr/bin/env gnuplot
# Set output format
set term pngcairo size 720,640
# Input file contains comma-separated values fields
set datafile separator ","
set hidden3d
set pm3d
set grid

# Where to write to
set output "output_data/heat.png"
set title "Diffusion"
set ylabel "Time"
set xlabel "Space"
set zlabel "Temperature?"
# set xrange [0:7.5]
# set yrange [-0:1.01]
# set log z
set view 60,360-50-90
set ticslevel 0
splot "output_data/res_LU.dat" matrix nonuniform notitle
# Write to file
set output

#!/usr/bin/env gnuplot
# Set the output as latexcode
# set terminal epslatex standalone
# Input file contains comma-separated values fields
set datafile separator ","
# set grid
set format xy "$%g$"

# Where to write to
# set output "output_data/nice_neutron.tex"
set title "$\\theta$ for stable values fo $n$"
set xlabel "$x$"
set ylabel "$\\theta_n(x)$"
# set xrange [0:7.5]
# set yrange [-0:1.01]
set hidden3d
set pm3d
# set log z
splot "output_data/res.dat" matrix
# Write to file
# set output
# Make the pdf
# !pdflatex -output-directory=output_data output_data/nice_neutron.tex

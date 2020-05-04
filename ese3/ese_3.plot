#!/usr/bin/env gnuplot
# Set output format
set term pngcairo size 1024,768
# Input file contains comma-separated values fields
set datafile separator ","
set hidden3d
set pm3d
set grid

E=10
r=0.1
sigma=0.4
T=0.25
a=2*r/sigma**2-1
b=-2*r/sigma**2
q=a/2

# # Where to write to
# set output "output_data/plot1.png"
# set title "Diffusion"
# set ylabel "Time"
# set xlabel "Space"
# set zlabel "Temperature?"
# # set xrange [0:7.5]
# # set yrange [-0:1.01]
# # set log z
# set view 60,360-50
# set ticslevel 0
# splot "output_data/res_LU.csv" matrix nonuniform notitle
# # Write to file
# set output

# Second plot
reset
set datafile separator ","
set output "output_data/plot2.png"
set title "Comparison"
set ylabel "V"
set xlabel "S"
set xrange [0:3*E]
set yrange [-1:*]
phi(x) = 1.0/2*(1+erf(x/sqrt(2)))
f(x) = E*exp(-(a/2+1)**2*sigma**2/2*T)*( \
	exp(sigma**2/8*T*(a/2+1)**2)*x/E*phi( (log(x/E)+T*sigma**2/2*(a/2+1))/sqrt(T*sigma**2/2) ) \
	+ exp(sigma**2/8*T*(a/2)**2)*phi( (log(x/E)+T*sigma**2/2*(a/2))/sqrt(T*sigma**2/2) ) )
plot f(x) title "analitic solution" w line \
, "output_data/initialstatus.csv" title "initial status" w line\
, "output_data/finalstatusexplicit.csv" title "final status explicit" w line \
, "output_data/finalstatusLU.csv" title "final status LU" w line

set grid
set output

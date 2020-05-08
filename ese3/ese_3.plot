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
t=0

# Where to write to
set output "output_data/plot1.png"
set title "Numerical solution"
set ylabel "t"
set xlabel "S"
set zlabel "V"
# set xrange [0:7.5]
# set yrange [-0:1.01]
set zrange [0:*]
# set log z
set view 60,360-50
set ticslevel 0
splot "output_data/res_LU.csv" matrix nonuniform notitle w l
# Write to file
set output

# Second plot
reset
set grid
set datafile separator ","
set output "output_data/plot2.png"
set title "Comparison"
set ylabel "V"
set xlabel "S"
set xrange [0:3*E]
set yrange [-1:*]
set key left top
phi(x) = 1.0/2*(1+erf(x/sqrt(2)))
#f(x) = E*exp(-(a/2+1)**2*sigma**2*T)*( \
#	exp(sigma**2*T*(a/2+1)**2)*x/E*phi( (log(x/E)+T*sigma**2*(a/2+1))/sqrt(T*sigma**2) ) \
#	- exp(sigma**2*T*(a/2)**2)*phi( (log(x/E)+T*sigma**2*(a/2))/sqrt(T*sigma**2) ) )

f(x)=E*exp(-(a/2+1)**2*sigma*sigma*(T-t)/2)*(exp(sigma*sigma*(T-t)/2*(a/2+1)**2)*x/E*(1.0/2*(1+erf((log(x/E)+(T-t)*sigma*sigma*(a/2+1))/sqrt((T-t)*sigma*sigma)/sqrt(2))))- exp(sigma*sigma*(T-t)/2*(a/2)**2)*(1.0/2*(1+erf((log(x/E)+(T-t)*sigma*sigma*(a/2))/sqrt((T-t)*sigma*sigma)/sqrt(2)))))

plot f(x) title "analitic solution" w line \
, "output_data/initialstatus.csv" title "initial status" w line\
, "output_data/finalstatusexplicit.csv" title "final status explicit" w line \
, "output_data/finalstatusexplicit_more.csv" title "final status explicit finer grid" w line \
, "output_data/finalstatusLU.csv" title "final status LU" w line \
, "output_data/finalstatusLU_more.csv" title "final status LU finer grid" w line

set grid
set output



set output "output_data/plot3.png"
set title "Comparison, difference from analitical"
set ylabel "V"
set xlabel "S"
set xrange [0:3*E]
set yrange [*:*]
set key left bottom
plot 0 title "analitic solution" w line \
, "output_data/finalstatusexplicit.csv" u 1:($2-f($1)) title "final status explicit" w line \
, "output_data/finalstatusexplicit_more.csv" u 1:($2-f($1)) title "final status explicit finer grid" w line \
, "output_data/finalstatusLU.csv" u 1:($2-f($1)) title "final status LU" w line \
, "output_data/finalstatusLU_more.csv" u 1:($2-f($1)) title "final status LU finer grid" w line

set grid
set output

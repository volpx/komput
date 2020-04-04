#!/usr/bin/env gnuplot
# set the output as latexcode
set terminal epslatex standalone
# Input file contains comma-separated values fields
set datafile separator ","
set grid
set format xy "$%g$"

set output "output_data/nice_neutron.tex"
set title "$\\theta$ for stable values fo $n$"
set xlabel "$x$"
set ylabel "$\\theta_n(x)$"
set xrange [0:7.5]
set yrange [-0:1.01]
plot "output_data/neutron_n_1.50.csv" using 1:2 with lines title "$n$: $1.50$",\
    "output_data/neutron_n_1.75.csv" using 1:2 with lines title "$n$: $1.75$",\
    "output_data/neutron_n_2.00.csv" using 1:2 with lines title "$n$: $2.00$",\
    "output_data/neutron_n_2.25.csv" using 1:2 with lines title "$n$: $2.25$",\
    "output_data/neutron_n_2.50.csv" using 1:2 with lines title "$n$: $2.50$",\
    "output_data/neutron_n_2.75.csv" using 1:2 with lines title "$n$: $2.75$",\
    "output_data/neutron_n_3.00.csv" using 1:2 with lines title "$n$: $3.00$",
# write to file
set output
# make the pdf
!pdflatex -output-directory=output_data output_data/nice_neutron.tex
# reset the session for other plots
# reset

set output "output_data/notnice_neutron.tex"
set title "$\\theta$ for un-stable values fo $n$"
set xrange [0:40]
set yrange [-0.4:1]
set xlabel "$x$"
set ylabel "$\\theta_n(x)$"
plot "output_data/neutron_n_1.00.csv" using 1:2 with lines title "$n$: $1.00$",\
   "output_data/neutron_n_4.50.csv" using 1:2 with lines title "$n$: $4.50$",
set output
!pdflatex -output-directory=output_data output_data/notnice_neutron.tex

set output "output_data/mass_radius_relation.tex"
set title "$M(R)$ relation with $n=3/2$"
set xrange [20:55]
set yrange [0:4]
set xlabel "$R$ [km]"
set ylabel "$M(R)$ [M$_\\odot$]"
plot "output_data/mass_radius_1_5.csv" using (column(1)/1000):2 with lines title "$M(R)$ relation",
set output
!pdflatex -output-directory=output_data output_data/mass_radius_relation.tex

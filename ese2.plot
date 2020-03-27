#!/usr/bin/env gnuplot

# Input file contains comma-separated values fields
set datafile separator ","

set terminal pdf
set output "output_data/nice_neutron.pdf"
set title "theta nice values"
set xrange [0:7.5]
set yrange [0:1]
plot "output_data/neutron_n_1.50.csv" using 1:2 with lines title "n: 1.50",\
    "output_data/neutron_n_1.75.csv" using 1:2 with lines title "n: 1.75",\
    "output_data/neutron_n_2.00.csv" using 1:2 with lines title "n: 2.00",\
    "output_data/neutron_n_2.25.csv" using 1:2 with lines title "n: 2.25",\
    "output_data/neutron_n_2.50.csv" using 1:2 with lines title "n: 2.50",\
    "output_data/neutron_n_2.75.csv" using 1:2 with lines title "n: 2.75",\
    "output_data/neutron_n_3.00.csv" using 1:2 with lines title "n: 3.00",\

set terminal pdf
set output "output_data/notnice_neutron.pdf"
set title "theta not nice values"
set xrange [0:40]
set yrange [-0.4:1]
plot "output_data/neutron_n_1.00.csv" using 1:2 with lines title "n: 1.00",\
    "output_data/neutron_n_4.50.csv" using 1:2 with lines title "n: 4.50",

#!/usr/bin/env gnuplot

# Input file contains comma-separated values fields
set datafile separator ","

set title "GL,LL integration"
plot "output_data/GL-LL.csv" using 1:2 with lines title "Ng_{dumb}",\
    "output_data/GL-LL.csv" using 1:3 with lines title "Nl_{dumb}",\
    "output_data/GL-LL.csv" using 1:4 with lines title "Ng",\
    "output_data/GL-LL.csv" using 1:5 with lines title "Nl"
#set terminal pdf



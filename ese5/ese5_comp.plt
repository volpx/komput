set terminal pdfcairo
set output "output_data/corr_comp.pdf"

set title "cVV"
set ylabel "cVV"
set xlabel "i"

# set xrange [0:500]

plot exp(-x/1600)\
,"output_data/cVV.dat" i 0 u 1:2 w lines title columnheader(1)\
, "output_data/cVV.dat" i 4 u 1:2 w lines title columnheader(1)\
, "output_data/cVV.dat" i 1 u 1:2 w lines title columnheader(1)\
, "output_data/cVV.dat" i 3 u 1:2 w lines title columnheader(1)\
, "output_data/cVV.dat" i 2 u 1:2 w lines title columnheader(1)


set output

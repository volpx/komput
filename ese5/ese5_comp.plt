set terminal pdfcairo
set output "output_data/corr_comp.pdf"

set title "cVV"
set ylabel "cVV"
set xlabel "i"

set xrange [0:20]

plot exp(-x/1)\
, "output_data/cVV.dat" i 3 u 1:2 w lines title columnheader(1)
# ,"output_data/cVV.dat" i 0 u 1:2 w lines title columnheader(1)\
# , "output_data/cVV.dat" i 1 u 1:2 w lines title columnheader(1)\
# , "output_data/cVV.dat" i 2 u 1:2 w lines title columnheader(1)\
# , "output_data/cVV.dat" i 5 u 1:2 w lines title columnheader(1)\
# , "output_data/cVV.dat" i 6 u 1:2 w lines title columnheader(1)\
# , "output_data/cVV.dat" i 7 u 1:2 w lines title columnheader(1)\
# , "output_data/cVV.dat" i 8 u 1:2 w lines title columnheader(1)\
# , "output_data/cVV.dat" i 9 u 1:2 w lines title columnheader(1)\
# , "output_data/cVV.dat" i 4 u 1:2 w lines title columnheader(1)

set output

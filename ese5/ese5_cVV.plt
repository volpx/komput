set terminal pdfcairo
set output "output_data/cVV.pdf"

set title "cVV"
set ylabel "cVV"
set xlabel "i"
plot "output_data/cVV.dat" u 1:2 w lines \
, exp(-x/25)

set output

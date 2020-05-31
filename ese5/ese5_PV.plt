set terminal pdfcairo
set output "output_data/PV.pdf"

set title "P-V"
set ylabel "P"
set xlabel "V"
plot "output_data/PV.dat" u 2:1 w linespoints

set output

set terminal pdfcairo
set output "output_data/Prho.pdf"

set title "P-rho"
set ylabel "P"
set xlabel "rho"
set logscale x
set grid
set yrange [-1:1]
# set xrange [:0.1]

set errorbars  linecolor "#0000ff"
plot "output_data/PV_total.dat" u 4:2:3 w yerror pointtype 7 pointsize 0.3

set output

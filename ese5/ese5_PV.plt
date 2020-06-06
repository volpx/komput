set terminal pdfcairo
# set output "output_data/PV_normal.pdf"
set output "output_data/PV.pdf"

set title "P-V"
set ylabel "P"
set xlabel "V"
set logscale x
set grid
set yrange [-1:1]
# set xrange [1000:]

set errorbars  linecolor "#0000ff"
plot "output_data/PV_total.dat" u 1:2:3 w yerror pointtype 7 pointsize 0.3

set output

set terminal pdfcairo
set output "output_data/Cvrho.pdf"

set title "Cv-rho"
set ylabel "Cv"
set xlabel "rho"
# set logscale x
set grid
# set yrange [0:0.28]
set xrange [:0.2]

set errorbars  linecolor "#0000ff"
plot "output_data/Cv_total.dat" u 1:(column(2)/125):(column(3)/125) w yerror pointtype 7 pointsize 0.3

set output

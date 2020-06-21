set terminal pdfcairo
set output "output_data/Cvrho.pdf"
# set output "output_data/CvV.pdf"

set title "Cv-rho"
# set title "Cv-V"
set ylabel "Cv N^{-1} K_b^{-1}"
set xlabel "rho"
# set xlabel "V"
# set logscale x
set grid
# set yrange [0:0.28]
# set xrange [:0.2]

set errorbars  linecolor "#0000ff"
plot "output_data/Cv_total.dat" u 1:(column(2)/125):(column(3)/125) title "Cv" w yerror pointtype 7 pointsize 0.3
# plot "output_data/Cv_total.dat" u (125/column(1)):(column(2)/125):(column(3)/125) title "Cv" w yerror pointtype 7 pointsize 0.3

set output

set terminal pdfcairo
# set output "output_data/PV_normal.pdf"
set output "output_data/PV.pdf"

set x2tics nomirror
set x2label "rho"
set format x2 "%4.2e"
set format x "%4.1e"
set xtics nomirror
set title "P-V"
set ylabel "P"
set xlabel "V"
# set logscale x
# set logscale x2
set grid
set yrange [*:2]
set x2range [*:*] reverse
# set xrange [800:1.3e5]

set link x2 via 125/x inverse 125/x

set errorbars  linecolor "#0000ff"
plot "output_data/PV_total.dat" u 1:2:3 w yerror pointtype 7 pointsize 0.3 \
 ,"output_data/PV_total.dat" u 4:2 pointtype 0 notitle axes x2y1

set output

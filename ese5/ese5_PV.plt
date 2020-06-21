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

# set errorbars  linecolor "#0000ff"
# set yerrorbars 1 linecolor "#0000ff"

# plot "output_data/PV_total.dat" u 1:2:3 w yerror pointtype 7 pointsize 0.3 title "P"\
#  ,"output_data/PV_total.dat" u 4:2 pointtype 0 notitle axes x2y1 \
#   ,"output_data/PV_total.dat" u 1:(125/column(1)) w points linecolor "#00ff00" title "PV=NRT"\
#  ,"output_data/PV_total.dat" u 1:(-1./3/column(1)*column(5)):(1./3/column(1)*sqrt(column(6))) w yerror pointtype 6 pointsize 0.3 linecolor "#00ff00" title "W"

plot "output_data/PV_total.dat" u 1:2:3 w yerrorbars pointtype 7 pointsize 0.3 title "P"\
 ,"output_data/PV_total.dat" u 4:2 pointtype 0 notitle axes x2y1 \
  ,"output_data/PV_total.dat" u 1:(125/column(1)) w points linecolor "#00ff00" title "PV=NRT" \
 ,"output_data/PV_total.dat" u 1:(-1./3/column(1)*column(5)):(1./3/column(1)*sqrt(column(6))) w yerrorbars  pointtype 6 pointsize 0.3 title "W"

set output

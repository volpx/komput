set terminal gif animate delay 0.1 size 1024*2,768*2
set output 'output_data/animation.gif'
stats 'output_data/positions.dat' nooutput
set title "Positions"
set xlabel "X"
set ylabel "Y"
set zlabel "Z"
set xrange [0:12.7]
set yrange [0:12.7]
set zrange [0:12.7]

do for [i=1:int(STATS_blocks):100] {
	splot 'output_data/positions.dat' u 2:3:4 index i w points pointtype 1 notitle
	print i
}

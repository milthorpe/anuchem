set term post

set yrange[0:2000]
set ytics 0,500,2000
set ylabel 'time (ms)'

set xlabel 'loop size'

set border 3
set xtics nomirror
set ytics nomirror

set output 'scaling.ps'
set title "Scaling of parallel loop with size for different activity generation approaches"

plot "scaling.txt" using 1:2 with linespoints pt 5 title "standard", "scaling.txt" using 1:3 with linespoints pt 5 title "iterator", "scaling.txt" using 1:4 with linespoints pt 5 title "bisect", "scaling.txt" using 1:5 with linespoints pt 5 title "chunked"




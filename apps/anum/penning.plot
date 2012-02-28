set title "Ion in a penning trap"

set output 'penning.ps'
set term post

set xrange[0:40000]
set xtics 0,5000,40000
set xlabel 'time(fs)'

set ylabel 'displacement'

plot 'penning.dat' using 1:2 with linespoints title "z" lt rgb "red"

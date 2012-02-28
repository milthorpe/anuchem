set title "Ion in a penning trap"

set output 'penning.ps'
set term post

set xrange[0:40]
set xtics 0,5,40
set xlabel 'time(fs)'

set ylabel 'displacement'

plot 'penning.dat' using 1:3 with linespoints title "y" lt rgb "dark-blue", 'penning.dat' using 1:4 with linespoints title "y" lt rgb "red"

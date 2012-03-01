set title "Ion in a penning trap"

set output 'penning.ps'
set term post

set xlabel 'time(fs)'
set ylabel 'displacement'

plot 'penning.dat' using 1:2 with lines title "y" lt rgb "red", 'penning.dat' using 1:3 with lines title "z" lt rgb "blue"

set title "Acetaldehyde cation in a Penning trap"

set output 'penning.ps'
set term post

set xlabel 'time(us)'
set ylabel 'displacement (mm)'

plot 'penning.dat' using 1:($2/1000000) with lines title "x" lt rgb "red", 'penning.dat' using 1:($3/1000000) with lines title "y" lt rgb "blue"

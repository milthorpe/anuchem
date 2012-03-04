set title "Acetaldehyde cation in a Penning trap"

set output 'penning.ps'
set term post enhanced

set xlabel 'time({/Symbol m}s)'
set ylabel 'displacement (mm)'

plot 'penning.dat' using ($1/1000):($2/1000) with lines title "x" lt rgb "red", 'penning.dat' using ($1/1000):($3/1000) with lines title "y" lt rgb "blue"

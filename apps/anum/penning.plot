set title "Acetaldehyde cation (m/q = 43.04) in a Penning trap (B = 0.7646T V_T = 1V)"

set output 'penning.ps'
set term post enhanced

set xlabel 'time({/Symbol m}s)'
set ylabel 'displacement (mm)'

plot 'penning.dat' using ($1/1000):($2) with lines title "x" lt rgb "red", 'penning.dat' using ($1/1000):($3) with lines title "y" lt rgb "blue"

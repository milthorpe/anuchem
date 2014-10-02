set title "FTICR-MS: packet of Cs^+ (m/q = 132.9) and x^+ (m/q = 150) B = 7.0 V_T = 1V"

set output 'penning.ps'
set term post enhanced

set xlabel 'time({/Symbol m}s)'
set ylabel 'displacement (mm)'

set xzeroaxis

plot 'penning.dat' using ($1/1000):($3) with lines title "x" lt rgb "red", 'penning.dat' using ($1/1000):($4) with lines title "y" lt rgb "blue"

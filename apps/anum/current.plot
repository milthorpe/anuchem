set title "FTICR-MS induced current: packet of Cs^+ (m/q = 132.9) and x^+ (m/q = 150)"

set output 'current.ps'
set term post enhanced

set xlabel 'time({/Symbol m}s)'
set ylabel 'current (pA)'

set xzeroaxis

plot 'penning.dat' using ($1/1000):($9) with lines notitle lt rgb "blue"

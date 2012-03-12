set title "FTICR-MS induced current: single CH_3CO and HCO ions"

set output 'current.ps'
set term post enhanced

set xlabel 'time({/Symbol m}s)'
set ylabel 'current (aA)'

plot 'penning.dat' using ($1/1000):($8) with lines notitle lt rgb "blue"

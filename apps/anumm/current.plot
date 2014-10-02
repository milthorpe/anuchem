set title "FTICR-MS induced current"

set output 'current.ps'
#set term png enhanced font 'arial,22' size 1000,700
set term post enhanced

set xlabel 'time({/Symbol m}s)'
set ylabel 'current (pA)'

set xrange [1:800]

set xzeroaxis

plot 'current_2500000_3.dat' using 1 with lines notitle lt rgb "blue"

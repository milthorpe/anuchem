set title "Harmonic and Morse oscillators: bond length over time"

set output 'oscillator.ps'
set term post

set xrange[0:40]
set xtics 0,5,40
set xlabel 'time(fs)'

set yrange[0.07:0.13]
set ytics 0.08,0.01,0.13
set ylabel 'H-F distance (nm)'

equilibriumBondLength = 0.09169

plot 'harmonic.dat' using 1:2 with linespoints title "Harmonic" lt rgb "dark-blue", 'morse.dat' using 1:2 with linespoints title "Morse" lt rgb "red", equilibriumBondLength lt -1 title "Equilibrium bond length"

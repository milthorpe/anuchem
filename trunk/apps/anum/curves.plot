set title "Morse and harmonic curves"

set output 'curve_comparison.ps'
set term post

set xrange[0:0.25]
set xtics 0,0.05,0.25
set xlabel 'r(nm)'

set yrange[-100:1400]
set ytics -100,200,1400
set ylabel 'V(kJ/mol)'

f(x) = 0.5 * 582000 * (x-0.09169)**2
ff(x) = 582000 * (x-0.09169)
g(x) = 569.87 * (1 - exp(-(2/0.09169) * (x-0.09169)))**2
#2D(a^2(r-b)*e^(-a(r-b)))
gg(x) = 2 * 569.87 * (2/0.09169)**2 * (x-0.09169) * exp(-(2/0.09169) * (x-0.09169))

plot f(x) title "Harmonic" lt rgb "light-blue", ff(x) title "Harmonic derivative" lt rgb "dark-blue", g(x) title "Morse" lt rgb "forest-green", gg(x) title "Morse derivative" lt rgb "green"

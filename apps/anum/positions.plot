set output 'positions.ps'
set term post enhanced

set xrange [-23.5:23.5]
set yrange [-23.5:23.5]

set size square

set xzeroaxis
set yzeroaxis
set key left top

plot 'positions.dat' using 4:5 with points pt 0 notitle

set output 'positions.ps'
set term post enhanced

set xrange [-10:10]
set yrange [-10:10]

set size square

set xzeroaxis
set yzeroaxis
set key left top

plot 'positions_2500000_3.dat' using 4:5 with points pt 0 notitle

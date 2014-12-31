set output 'positions_25000_yz.eps'
set term postscript eps 22  enhanced

set xrange [-5:5]
set yrange [-5:5]

set size square
set palette model RGB defined ( 0 '#D95F02', 1 '#1B9E77' )
unset key
unset colorbox

set xzeroaxis
set yzeroaxis
set key left top

plot 'positions_25000.dat' using 4:5:2 with points pt 0 palette notitle

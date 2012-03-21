set output 'positions_stills.ps'
set term post enhanced

set xrange [-23.5:23.5]
set yrange [-23.5:23.5]

set noxtic
set noytic
set rmargin 0
set lmargin 0
set tmargin 0
set bmargin 0
set xzeroaxis
set yzeroaxis
set key left top

set multiplot
set size square 0.35,0.35
set origin 0,0.5
plot 'positions_0.dat' using 4:5 with points pt 0 title "0{/Symbol m}s"
set origin 0.33,0.5
plot 'positions_8000.dat' using 4:5 with points pt 0 title "8{/Symbol m}s"
set origin 0.66,0.5
plot 'positions_16000.dat' using 4:5 with points pt 0 title "16{/Symbol m}s"
set origin 0,0.1
plot 'positions_24000.dat' using 4:5 with points pt 0 title "24{/Symbol m}s"
set origin 0.33,0.1
plot 'positions_32000.dat' using 4:5 with points pt 0 title "32{/Symbol m}s"
set origin 0.66,0.1
plot 'positions_40000.dat' using 4:5 with points pt 0 title "40{/Symbol m}s"

set nomultiplot


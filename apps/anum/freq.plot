set title "FTICR-MS frequency spectrum: packet of Cs^+ (m/q = 132.9) and x^+ (m/q = 150)"

set output 'freq.png'
#set term post enhanced
set term png enhanced font 'arial,24' size 1000,700

set xlabel 'frequency (Hz)'
set ylabel 'amplitude (arbitrary units)'

set logscale y 10

set xrange [700000:820000]

plot 'penning.mz' using 1:2 with lines notitle lt rgb "blue"

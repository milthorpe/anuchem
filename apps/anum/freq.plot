set title "FTICR-MS frequency spectrum: packet of Cs^+ (m/q = 132.9) and x^+ (m/q = 150)"

set output 'freq.ps'
set term post enhanced

set xlabel 'frequency (Hz)'
set ylabel 'amplitude (arbitrary units)'

set xrange [600000:900000]

plot 'penning.mz' using 1:2 with lines notitle lt rgb "blue"

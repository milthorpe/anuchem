set title "FTICR-MS frequency spectrum: glutamine (147.07698) and lysine (147.11336)"

set output 'freq_aminos.ps'
set term post enhanced

set xlabel 'frequency (Hz)'
set ylabel 'amplitude (arbitrary units)'

set xrange [489000:490000]

plot 'penning.mz' using 1:($2/1000) with lines notitle lt rgb "blue"

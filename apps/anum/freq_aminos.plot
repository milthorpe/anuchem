set title "FTICR-MS frequency spectrum:\nglutamine (147.07698) and lysine (147.11336)"

set output 'freq_aminos.ps'
set term post enhanced
#set term png enhanced font 'arial,22' size 1000,700

set xlabel 'frequency (Hz)'
set ylabel 'amplitude (arbitrary units)'

set xrange [489000:490050]

set label "489406.59 Hz" at 489250,78
set label "489540.10 Hz" at 489350,114
plot 'penning.mz' using 1:($2/1000) with lines notitle lw 2 lt rgb "blue"


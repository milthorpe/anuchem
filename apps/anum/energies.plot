set title "Kinetic energy distribution in ion packet after 0.25ms, p=3"

set output 'energies_3.ps'
set terminal postscript enhanced 22

#to put an empty boundary around the
#data inside an autoscaled graph.

binwidth=0.00001
bin(x,width)=width*floor(x/width)

set boxwidth binwidth*0.9
set style fill solid 0.5	#fillstyle
set tics out nomirror

set xlabel 'energy (fJ)'
set ylabel 'frequency'

mean = 0.010448
stddev = 0.0006678

set arrow from mean,0 to mean,100 nohead lt 1 lc rgb "red"
set arrow from mean-stddev,0 to mean-stddev,100 nohead lt 2 lc rgb "green" 
set arrow from mean+stddev,0 to mean+stddev,100 nohead lt 2 lc rgb "green" 

plot "./energies_250000_3.dat" using (bin($4,binwidth)):(1.0) smooth freq with boxes lc rgb "blue" notitle



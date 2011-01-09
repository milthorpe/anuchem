set terminal png size 1479,986 crop
set output 'mandelbrot.png'

set notics
set nokey
unset colorbox


set pm3d map
set palette defined (0 "black", 1 "navy", 9 "blue", 25 "dark-cyan", 59 "green", 255 "white")
splot 'mandelbrot.dat' u 2:1:3 matrix with pm3d



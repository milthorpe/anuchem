set terminal png size 1479,986 crop
set output 'mandelbrot.png'

set notics
set nokey
unset colorbox

set palette defined (0 "black", 1 "navy", 16 "blue", 36 "dark-cyan", 81 "green", 255 "white")
plot 'mandelbrot.dat' u 2:1:3 matrix with image



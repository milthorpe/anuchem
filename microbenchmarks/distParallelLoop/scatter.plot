set term post

set yrange[-0.5:3.5]
set ytics 0,1,3
set xtics 0,1000,8000

set border 3
set xtics nomirror
set ytics nomirror

set pointsize 0.5

set output 'standard.ps'
set title "Array elements processed by each worker thread: X10_NTHREADS=4 : standard for/async loop"
plot "standard.dat" using 1:2 with linespoints pt 5 notitle

set output 'iterator.ps'
set title "Array elements processed by each worker thread: X10_NTHREADS=4 : iterator loop"
plot "iterator.dat" using 1:2 with linespoints pt 5 notitle

set output 'bisect.ps'
set title "Array elements processed by each worker thread: X10_NTHREADS=4 : loop bisection"
plot "bisect.dat" using 1:2 with linespoints pt 5 notitle

set output 'part_bisect.ps'
set title "Array elements processed by each worker thread: X10_NTHREADS=4 : partial loop bisection"
plot "part_bisect.dat" using 1:2 with linespoints pt 5 notitle

set output 'chunked.ps'
set title "Array elements processed by each worker thread: X10_NTHREADS=4 : chunked loop"
plot "chunked.dat" using 1:2 with linespoints pt 5 notitle




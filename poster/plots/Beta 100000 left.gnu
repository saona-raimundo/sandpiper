imagefile = 'images\\beta_100000_left'

unset key
set title "Beta = 100000"
set xlabel "selection, s"
set ylabel "dominance, h"
set logscale x 10
# set logscale y 10
set xtics 
set ytics 
set yrange [-0.05:1.05]

set xrange [*:*] reverse
set format x "-%f" # "-%.0s*10^{%T}"

set terminal svg font 'Verdana,10' # size 550,262 
set output sprintf("%s.svg", imagefile)

plot "data\\Beta 100000.txt" using 1:(1.0 - $2) with linespoints dashtype 1 linecolor rgb 'black', 

system(sprintf("inkscape --export-area-drawing --export-filename=%s.png --export-dpi=300 %s.svg", imagefile, imagefile))
imagefile = 'images\\fixed_left'

set title "Negative selection"
set xlabel "selection, s"
set ylabel "polymorphisms"
set logscale x 10
set logscale y 10
set xtics 
set ytics 

set xrange [*:*] reverse
set format x "-%f"

set terminal svg font 'Verdana,10' # size 550,262 
set output sprintf("%s.svg", imagefile)

plot "data\\fixed dominance_4.txt" using 1:2 with linespoints title "h = 0.0" dashtype 1, \
"data\\fixed dominance_3.txt" using 1:2 with linespoints title "h = 0.2" dashtype 2, \
"data\\fixed dominance_2.txt" using 1:2 with linespoints title "h = 0.5" dashtype 3, \
"data\\fixed dominance_1.txt" using 1:2 with linespoints title "h = 0.8" dashtype 4, \
"data\\fixed dominance_0.txt" using 1:2 with linespoints title "h = 1.0" dashtype 5, 

# pause -1

system(sprintf("inkscape --export-area-drawing --export-filename=%s.png --export-dpi=900 %s.svg", imagefile, imagefile))
imagefile = 'images\\fixed_right'

set title "Positive selection"
set xlabel "selection, s"
set ylabel ""
set logscale x 10
set logscale y 10
set xtics 
unset ytics 

set terminal svg font 'Verdana,10' # size 550,262 
set output sprintf("%s.svg", imagefile)

plot "data\\fixed dominance_0.txt" using 1:2 with linespoints title "h = 0.0" dashtype 1, \
"data\\fixed dominance_1.txt" using 1:2 with linespoints title "h = 0.2" dashtype 2, \
"data\\fixed dominance_2.txt" using 1:2 with linespoints title "h = 0.5" dashtype 3, \
"data\\fixed dominance_3.txt" using 1:2 with linespoints title "h = 0.8" dashtype 4, \
"data\\fixed dominance_4.txt" using 1:2 with linespoints title "h = 1.0" dashtype 5, 

system(sprintf("inkscape --export-area-drawing --export-filename=%s.png --export-dpi=300 %s.svg", imagefile, imagefile))

# pause -1
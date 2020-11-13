imagefile = 'images\\beta_5000_right'

unset key
set title "Beta = 5000"
set xlabel "selection, s"
set ylabel ""
set logscale x 10
# set logscale y 10
set xtics 
set ytics 
set yrange [-0.05:1.05]

set terminal svg font 'Verdana,10' # size 550,262 
set output sprintf("%s.svg", imagefile)

plot "data\\Beta 5000.txt" using 1:2 with linespoints dashtype 1 linecolor rgb 'black', 

system(sprintf("inkscape --export-area-drawing --export-filename=%s.png --export-dpi=300 %s.svg", imagefile, imagefile))
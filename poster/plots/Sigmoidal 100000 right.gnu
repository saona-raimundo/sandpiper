rate = "100000"
side = "right"
imagefile = 'images\\fixed_right'


# unset key
set title sprintf("Beta = %s", rate)
set xlabel "selection, s"
set ylabel ""
set logscale x 10
set logscale y 10
set yrange [0.00000001:0.1]
set xtics 
unset ytics 

set terminal svg font 'Verdana,10' # size 550,262 
set output sprintf("images\\sigmoidal_%s_%s.svg", rate, side)

plot sprintf("data\\Sigmoidal %s.txt", rate) using 1:2 with linespoints title "Sigmoidal" dashtype 1 linecolor rgb "black" , \
	"data\\fixed dominance_0.txt" using 1:2 with linespoints title "h = 0.0" dashtype 2 pointtype 2 linecolor rgb "grey" , \
	"data\\fixed dominance_4.txt" using 1:2 with linespoints title "h = 1.0" dashtype 2 pointtype 2 linecolor rgb "grey" ,

# pause -1

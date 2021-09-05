unset key
set title ""
set xlabel ""
set ylabel ""
set xtics 
set ytics 
set logscale y 10
set pm3d map
splot "examples\\ploting_results\\redneck.txt" using 1:2:3
pause -1
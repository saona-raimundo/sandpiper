g(s, beta) = 1 / (1 + exp(-beta * s))

set xrange [-0.00014:0.00002]
set xlabel "selection, s"
set ylabel "dominance, h"

set title "Dominance coefficient for various beta"

plot g(x, 0) title "beta = 0",\
	g(x, 10) title "beta = 10^1",\
	g(x, 100) title "beta = 10^2",\
	g(x, 1000) title "beta = 10^3",\
	g(x, 10000) title "beta = 10^4",\
	g(x, 100000) title "beta = 10^5",\

pause -1
# Sandpiper
N = 2 * 10**4
U = 1.2 * 10**(-8)

h(s, beta) = 1 / (1 + exp(-beta * s))

g(x, s, beta) = exp(2 * N * s * ( x*x + 2 *h(s, beta) *x *(1 - x) - (s >= 0))) * x**(4 * N * U - 1) * (1 - x)**(4 * N * U - 1)


set xrange [0:1]
set xlabel "frequency, x"
set ylabel "unnormalized density"

set title sprintf("Sandpiper, beta = %.f", beta)

s = -1e-4

plot g(x, s, 0) title "beta = 0",\
	g(x, s, 10) title "beta = 10",\
	g(x, s, 100) title "beta = 100",\
	g(x, s, 1000) title "beta = 1000",\
	g(x, s, 10000) title "beta = 10000",\
	g(x, s, 1000000) title "beta = 1000000"

pause -1
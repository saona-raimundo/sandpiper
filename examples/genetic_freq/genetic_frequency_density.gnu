# Sandpiper
N = 2 * 10**4
U = 1.2 * 10**(-8)
beta = 0

h(s) = 1 / (1 + exp(-beta * s))

g(x, s) = exp(2 * N * s * ( x*x + 2 *h(s) *x *(1 - x) - (s >= 0))) * x**(4 * N * U - 1) * (1 - x)**(4 * N * U - 1)


set xrange [0:1]
set xlabel "frequency, x"
set ylabel "unnormalized density"

set title sprintf("Sandpiper, beta = %.f", beta)

plot g(x, 0.001) title "s = 0.001"
#g(x, 0) title "s = 0",\
#	g(x, 0.001) title "s = 0.001",\
#	g(x, -0.001) title "s = -0.001"

pause -1
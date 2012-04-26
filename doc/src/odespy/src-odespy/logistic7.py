def f(y, x):
    return a*y*(1 - y/R)

a = 2;  R = 1E+5;  A = 1

import odespy, numpy
solver = odespy.RK4(f)
solver.set_initial_condition(A)

L = 10  # end of x domain
N = 30  # no of time steps
x_points = numpy.linspace(0, L, N+1)
y, x = solver.solve(x_points)

from matplotlib.pyplot import *
plot(x, y, 'r-')
xlabel('x'); ylabel('y')
show()

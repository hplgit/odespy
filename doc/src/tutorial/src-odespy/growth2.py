"""
As expgr1.py, but class for right-hand side function and
comparison with exact solution.
"""

class ExponentialGrowth:
    def __init__(self, c=1, A=1):
        self.c, self.A = c, A

    def __call__(self, u, t):
        return self.c*u

f = ExponentialGrowth(c=0.1, A=1.5)

import odespy
solver = odespy.RK4(f)
solver.set_initial_condition(f.A)

import numpy
N = 30  # no of time steps
time_points = numpy.linspace(0, 40, N+1)
u, t = solver.solve(time_points)

u_exact = f.A*numpy.exp(f.c*t)
error = numpy.abs(u_exact - u).max()
print 'Max deviation of numerical solution:', error

from matplotlib.pyplot import *
plot(t, u, 'r-', t, u_exact, 'bo')
legend(['RK4, N=%d' % N, 'exact'])
show()

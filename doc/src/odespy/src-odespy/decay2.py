"""As exde1.py, but y and x as names instead of u and t."""

def myrhs(y, x):
    return -y

import odespy
solver = odespy.RK4(myrhs)
solver.set_initial_condition(1)

import numpy
# Make sure integration interval [0, L] is large enough
N = 50
L = 10
x_points = numpy.linspace(0, L, N+1)

def terminate(y, x, stepnumber):
    tol = 0.001
    return True if abs(y[stepnumber]) < tol else False

y, x = solver.solve(x_points, terminate)

print 'Final y(x=%g)=%g' % (x[-1], y[-1])

from matplotlib.pyplot import *
plot(x, y)
show()

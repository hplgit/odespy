"""As exde1.py, but no terminate function and AdaptiveResidual as solver."""
c = -1
A = 1

def f(u, t):
    return c*u

import odespy
solver = odespy.AdaptiveResidual(f, solver='Euler', atol=1E-3)
solver.set_initial_condition(A)

import numpy
# Make sure integration interval [0, T] is large enough
N = 20
T = 4
time_points = numpy.linspace(0, T, N+1)

u, t = solver.solve(time_points)

print 'Final u(t=%g)=%g after %d steps' % (t[-1], u[-1], len(u)-1)

from matplotlib.pyplot import *
plot(t, u)
show()

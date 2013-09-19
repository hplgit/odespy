"""Simple example of exponential growth ODE."""
c = 0.1
A = 1.5

def f(u, t):
    return c*u

import odespy
solver = odespy.RK4(f)
solver.set_initial_condition(A)

import numpy
N = 30  # no of time steps
time_points = numpy.linspace(0, 40, N+1)
u, t = solver.solve(time_points)

from matplotlib.pyplot import *
plot(t, u)
show()

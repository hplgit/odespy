"""Oscillating pendulum."""

import odespy, numpy
from math import sin, pi, sqrt

c = 1
Theta = pi/4

def f(u, t):
    theta, omega = u
    return [omega, -c*sin(theta)]

solver = odespy.Heun(f)
solver.set_initial_condition([Theta, 0])
freq = sqrt(c)      # frequency of oscillations when Theta is small
period = 2*pi/freq  # the period of the oscillations
T = 10*period       # final time
N_per_period = 20   # resolution of one period
N = N_per_period*period
time_points = numpy.linspace(0, T, N+1)

u, t = solver.solve(time_points)

theta = u[:,0]
omega = u[:,1]

from matplotlib.pyplot import *
plot(t, theta, 'r-')
savefig('tmppng'); savefig('tmp.pdf')
show()


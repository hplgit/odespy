from numpy import linspace

def f(u, t, a=1, R=1):
    return a*u*(1 - u/R)

A = 1

import odespy
solver = odespy.RK4(f, f_kwargs=dict(a=2, R=1E+5))
solver.set_initial_condition(A)

T = 10  # end of simulation
N = 30  # no of time steps
time_points = linspace(0, T, N+1)
u, t = solver.solve(time_points)

from matplotlib.pyplot import *
plot(t, u, 'r-')
savefig('tmppng'); savefig('tmp.pdf')
show()

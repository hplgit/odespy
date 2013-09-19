def f(u, t, a=1, R=1):
    return a*u*(1 - u/R)

A = 1

import odespy
solver = odespy.CashKarp(f, f_kwargs=dict(a=2, R=1E+5),
                         first_step=9., verbose=1)
solver.set_initial_condition(A)

from numpy import linspace
T = 10  # end of simulation
N = 3   # no of time steps
time_points = linspace(0, T, N+1)
u, t = solver.solve(time_points)

from matplotlib.pyplot import *
plot(solver.t_all, solver.u_all, 'r-o')
savefig('tmppng'); savefig('tmp.pdf')
show()

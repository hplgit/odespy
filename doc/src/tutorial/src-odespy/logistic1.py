def f(u, t):
    return a*u*(1 - u/R)

a = 2
R = 1E+5
A = 1

import odespy
solver = odespy.RK4(f)
solver.set_initial_condition(A)

from numpy import linspace, exp
T = 10  # end of simulation
N = 30  # no of time steps
time_points = linspace(0, T, N+1)
u, t = solver.solve(time_points)

def u_exact(t):
    return R*A*exp(a*t)/(R + A*(exp(a*t) - 1))

from matplotlib.pyplot import *

plot(t, u, 'r-',
     t, u_exact(t), 'bo')
savefig('tmppng'); savefig('tmp.pdf')
show()

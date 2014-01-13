"""Gaussian function for testing adaptivity."""

from odespy.problems import Logistic
from numpy import exp, linspace

import odespy
problem = Logistic(a=2, R=1, U0=0.0001)
solver = odespy.RKFehlberg
#solver = odespy.lsoda_scipy
#solver = odespy.Lsoda
solver = odespy.RK4
solver = solver(problem.f, verbose=2, first_step=1.100,
                u_exact=problem.u_exact, atol=1E-6, rtol=1E-6)
solver.set_initial_condition(problem.u_exact(0))

T = 10
dt = 0.5
N = round(int(T/float(dt)))
time_points = linspace(0, T, N+1)
u, t = solver.solve(time_points)
print u

from matplotlib.pyplot import *
plot(t, u, 'r-o',
     t, problem.u_exact(t), 'x')
try:
    hold('on')
    plot(solver.t_all, solver.u_all, 'o')
except:
    pass
show()

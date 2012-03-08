"""Gaussian function for testing adaptivity."""

from odespy.problems import Problem
from numpy import exp, linspace

class Gaussian(Problem):
    def __init__(self, mu, sigma):
        self.mu, self.sigma = float(mu), float(sigma)

    def u_exact(self, t):
        return 0.1 + 1/sqrt(2*pi)*exp(-(t-self.mu)**2/self.sigma**2)

    def f(self, u, t):
        return -2*(t-self.mu)/self.sigma**2*u

class Decay(Problem):
    """u(t) = 1 - exp(-a*t):  u' = a*(1 - u)."""
    def __init__(self, a):
        self.a = a

    def f(self, u, t):
        return self.a*(1 - u)

    def u_exact(self, t):
        return 1 - exp(-self.a*t)

import odespy
problem = Decay(1)
solver = odespy.CashKarp(problem.f, verbose=2, first_step=10.0,
                     u_exact=problem.u_exact, atol=1E-6, rtol=1E-6)
solver.set_initial_condition(problem.u_exact(0))

T = 10
dt = 1.0
N = round(int(T/float(dt)))
time_points = linspace(0, T, N+1)
u, t = solver.solve(time_points)
print u

from matplotlib.pyplot import *
plot(solver.t_all, solver.u_all, 'o',
     t, u, 'b-')
#plot(t, u, 'r-o')
show()

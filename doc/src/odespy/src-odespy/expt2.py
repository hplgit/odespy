"""Gaussian function for testing adaptivity."""

from problems import Problem
from numpy import sqrt, exp, pi, linspace

class Gaussian(Problem):
    def __init__(self, mu, sigma):
        self.mu, self.sigma = float(mu), float(sigma)

    def u_exact(t):
        return 1/sqrt(2*pi)*exp(-(t-self.mu)**2/self.sigma**2)

    def __call__(self, u, t):
        return -2*(t-self.mu)/self.sigma**2*u

import odespy
solver = odespy.RK4(f)
solver.set_initial_condition(A)

problem = Gaussian(mu=5, sigma=0.5)
T = problem.mu + 4*problem.sigma
dt = 3*problem.sigma/100.
N = round(int(T/dt))
time_points = linspace(0, T, N+1)
u, t = solver.solve(time_points)

from matplotlib.pyplot import *
plot(t, u)
show()

from numpy import linspace, exp

class Logistic:
    def __init__(self, a, R, A):
        self.a = a
        self.R = R
        self.A = A

    def f(self, u, t):
        a, R = self.a, self.R  # short form
        return a*u*(1 - u/R)

    def u_exact(self, t):
        a, R, A = self.a, self.R, self.A  # short form
        return R*A*exp(a*t)/(R + A*(exp(a*t) - 1))

import odespy
problem = Logistic(a=2, R=1E+5, A=1)
solver = odespy.RK4(problem.f)
solver.set_initial_condition(problem.A)

T = 10  # end of simulation
N = 30  # no of time steps
time_points = linspace(0, T, N+1)
u, t = solver.solve(time_points)

from matplotlib.pyplot import *

plot(t, u, 'r-',
     t, problem.u_exact(t), 'bo')
savefig('tmppng'); savefig('tmp.pdf')
show()

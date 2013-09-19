"""Stochastic oscillator driven by white noise."""


class WhiteNoiseOscillator:
    def __init__(self, b, c, sigma=1):
        self.b, self.c, self.sigma = b, c, sigma

    def connect_solver(self, solver):
        """Solver is needed for time step number and size."""
        self.solver = solver

    def f(self, u, t):
        if not hasattr(self, 'N'):  # is self.N not yet computed?
            # Compute N(t) for all time intervals
            import numpy
            numpy.random.seed(12)
            t = self.solver.t
            dW = numpy.random.normal(loc=0, scale=1, size=len(t)-1)
            dt = t[1:] - t[:-1]
            self.N = self.sigma*dW/numpy.sqrt(dt)

        x, v = u
        N = self.N[self.solver.n]
        return [v, N -self.b*v -self.c*x]

from numpy import pi, linspace
from matplotlib.pyplot import *
import odespy

problem = WhiteNoiseOscillator(b=0.1, c=pi**2, sigma=1)
solvers = [odespy.Heun(problem.f), odespy.RK4(problem.f),
           odespy.ForwardEuler(problem.f)]
for solver in solvers:
    f.connect_solver(solver)
    solver.set_initial_condition([0,0])  # start from rest
    T = 60   # with c=pi**2, the period is 1
    u, t = solver.solve(linspace(0, T, 10001))

    x = u[:,0]
    plot(t, x)
    hold(True)

legend([str(s) for s in solvers])
savefig('tmppng'); savefig('tmp.pdf')
show()


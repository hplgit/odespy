"""Stochastic oscillator driven by white noise."""


class WhiteNoiseOscillator:
    def __init__(self, b, c, sigma=1):
        self.b, self.c, self.sigma = b, c, sigma

    def connect_solver(self, solver):
        """Solver is needed for time step number and size."""
        self.solver = solver

    def __call__(self, u, t):
        if not hasattr(self, 'N'):
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

f = WhiteNoiseOscillator(b=0.1, c=pi**2, sigma=1)
methods = [odespy.Heun(f), odespy.RK4(f),
           odespy.ForwardEuler(f)]
for method in methods:
    f.connect_solver(method)
    solver.set_initial_condition([0,0])  # start from rest
    T = 60   # with c=pi**2, the period is 1
    u, t = solver.solve(linspace(0, T, 10001))

    x = u[:,0]
    plot(t, x)
    hold(True)

legend([str(m) for m in methods])
savefig('tmp.png')
show()


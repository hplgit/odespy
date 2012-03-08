"""As exos1.py, but testing several solvers and setting sin(theta) to theta."""

from math import pi, sqrt

class Problem:
    def __init__(self, c, Theta):
        self.c, self.Theta = float(c), float(Theta)

        self.freq = sqrt(c)
        self.period = 2*pi/self.freq

    def __call__(self, u, t):
        theta, omega = u;  c = self.c
        return [omega, -c*theta]

problem = Problem(c=1, Theta=pi/4)

import odespy
solvers = [
    odespy.ThetaRule(problem, theta=0),   # Forward Euler
    odespy.ThetaRule(problem, theta=0.5), # Midpoint
    odespy.ThetaRule(problem, theta=1),   # Backward Euler
    odespy.RK4(problem),
    odespy.MidpointIter(problem, max_iter=2, eps_iter=0.01),
    odespy.LeapfrogFiltered(problem),
    ]

N_per_period = 20
T = 3*problem.period   # final time
import numpy
import matplotlib.pyplot as mpl
legends = []

for solver in solvers:
    solver_name = str(solver)
    print solver_name

    solver.set_initial_condition([problem.Theta, 0])
    N = N_per_period*problem.period
    time_points = numpy.linspace(0, T, N+1)

    u, t = solver.solve(time_points)

    theta = u[:,0]
    legends.append(solver_name)
    mpl.plot(t, theta)
    mpl.hold('on')
mpl.legend(legends)
mpl.savefig(__file__[:-3] + '.png')
mpl.show()


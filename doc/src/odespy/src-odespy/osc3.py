"""As exos2.py, but testing several methods."""

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

import odesolvers
solvers = [
    odesolvers.ThetaRule(problem, theta=0),   # Forward Euler
    odesolvers.ThetaRule(problem, theta=0.5), # Midpoint
    odesolvers.ThetaRule(problem, theta=1),   # Backward Euler
    odesolvers.RK4(problem),
    odesolvers.RK2(problem),
    odesolvers.MidpointIter(problem, max_iter=5, eps_iter=0.01),
    odesolvers.Leapfrog(problem),
    odesolvers.LeapfrogFiltered(problem),
    ]

theta_exact = lambda t: problem.Theta*numpy.cos(sqrt(problem.c)*t)

import sys
try:
    num_periods = int(sys.argv[1])
except IndexError:
    num_periods = 8  # default

T = num_periods*problem.period       # final time
results = {}
resolutions = [10, 20, 40, 80, 160]  # intervals per period
import numpy

for method in solvers:
    method_name = str(method)
    results[method_name] = {'dt': [], 'error': []}

    solver.set_initial_condition([problem.Theta, 0])

    for N_per_period in resolutions:
        N = N_per_period*problem.period
        time_points = numpy.linspace(0, T, N+1)

        u, t = solver.solve(time_points)

        theta = u[:,0]
        error = numpy.abs(theta_exact(t) - theta)
        error_L2 = sqrt(numpy.sum(error**2)/N)
        if not numpy.isnan(error_L2):  # drop nan
            results[method_name]['dt'].append(t[1] - t[0])
            results[method_name]['error'].append(error_L2)

# Print
for method in results:
    print '%-20s %s' % (method, ' '.join(['%.1E' % e
                        for e in results[method]['error']]))

# Analyze convergence
from scitools.convergencerate import OneDiscretizationPrm
pairwise_rates = OneDiscretizationPrm.pairwise_rates  # short form

print '\n\nConvergence results for %d periods' % num_periods
for method_name in results:
    rates, C = pairwise_rates(results[method_name]['dt'],
                              results[method_name]['error'])
    rates = ', '.join(['%.1f' % rate for rate in rates])
    print '%-20s r: %s E_min=%.1E' % \
          (method_name, rates, min(results[method_name]['error']))

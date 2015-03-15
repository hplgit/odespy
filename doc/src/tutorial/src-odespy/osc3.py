"""As osc2.py, but testing several solvers."""

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
    odespy.RK2(problem),
    odespy.MidpointIter(problem, max_iter=5, eps_iter=0.01),
    odespy.Leapfrog(problem),
    odespy.LeapfrogFiltered(problem),
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

for solver in solvers:
    solver_name = str(solver)
    results[solver_name] = {'dt': [], 'error': []}

    solver.set_initial_condition([problem.Theta, 0])

    for N_per_period in resolutions:
        N = N_per_period*num_periods
        time_points = numpy.linspace(0, T, N+1)

        u, t = solver.solve(time_points)

        theta = u[:,0]
        error = numpy.abs(theta_exact(t) - theta)
        error_L2 = sqrt(numpy.sum(error**2)/N)
        if not numpy.isnan(error_L2):  # drop nan (overflow)
            results[solver_name]['dt'].append(t[1] - t[0])
            results[solver_name]['error'].append(error_L2)

# Print
for solver in results:
    print '%-20s %s' % (solver, ' '.join(['%.1E' % e
                        for e in results[solver]['error']]))

from math import log
print '\n\nConvergence results for %d periods' % num_periods
for solver_name in results:
    r_h = results[solver_name]['dt']
    r_E = results[solver_name]['error']
    rates = [log(r_E[i]/r_E[i-1])/log(r_h[i]/r_h[i-1]) for i
             in range(1, len(r_h))]
    # Reformat rates with 1 decimal for rate
    rates = ', '.join(['%.1f' % rate for rate in rates])
    print '%-20s r: %s E_min=%.1E' % \
          (solver_name, rates, min(results[solver_name]['error']))

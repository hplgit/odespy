"""As exos3.py, but testing many more methods for a fixed resolution."""

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

atol = 1E-7
rtol = 1E-5
adams_or_bdf = 'bdf'
import odesolvers
solvers = [
    odesolvers.AdamsBashMoulton2(problem),
    odesolvers.AdamsBashMoulton3(problem),
    odesolvers.AdamsBashforth2(problem),
    odesolvers.AdamsBashforth3(problem),
    odesolvers.AdamsBashforth4(problem),
    odesolvers.AdaptiveResidual(problem, solver='Euler'),
    odesolvers.Backward2Step(problem),
    odesolvers.BackwardEuler(problem),
    odesolvers.Dop853(problem, rtol=rtol, atol=atol),
    odesolvers.Dopri5(problem, rtol=rtol, atol=atol),
    odesolvers.Euler(problem),
    odesolvers.Heun(problem),
    odesolvers.Leapfrog(problem),
    odesolvers.LeapfrogFiltered(problem),
    odesolvers.MidpointImplicit(problem),
    odesolvers.MidpointIter(problem, max_iter=10, eps_iter=1E-7),
    odesolvers.RK2(problem),
    odesolvers.RK3(problem),
    odesolvers.RK4(problem),
    odesolvers.RKFehlberg(problem, rtol=rtol, atol=atol),
    odesolvers.SymPy_odefun(problem),
    odesolvers.ThetaRule(problem),
    odesolvers.Trapezoidal(problem),
    odesolvers.Vode(problem, rtol=rtol, atol=atol,
                    adams_or_bdf=adams_or_bdf),
    ]

theta_exact = lambda t: problem.Theta*numpy.cos(sqrt(problem.c)*t)

import sys, time
try:
    num_periods = int(sys.argv[1])
except IndexError:
    num_periods = 30  # default

T = num_periods*problem.period       # final time
results = {}
N_per_period = 20
import numpy

for method in solvers:
    method_name = str(method)
    results[method_name] = {'dt': 0, 'error': numpy.nan, 'cpu': 0}

    solver.set_initial_condition([problem.Theta, 0])

    N = N_per_period*problem.period
    time_points = numpy.linspace(0, T, N+1)

    try:
        t0 = time.clock()
        u, t = solver.solve(time_points)
        cpu_time = time.clock() - t0
    except Exception, e:
        print method_name, 'FAILED!'
        print e
        continue  # continue with next pass in the loop

    theta = u[:,0]
    error = numpy.abs(theta_exact(t) - theta)
    error_L2 = sqrt(numpy.sum(error**2)/N)
    #print method_name, error_L2, cpu_time
    results[method_name]['dt'] = t[1] - t[0]
    results[method_name]['error'] = error_L2
    results[method_name]['cpu'] = cpu_time

# Print
for method_name in results:
    print '%-50s %.1E %.2f' % (method_name,
                               results[method_name]['error'],
                               results[method_name]['cpu'])

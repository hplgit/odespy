"""As osc3.py, but testing many more solvers for a fixed resolution."""

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
import odespy
solvers = [
    odespy.AdamsBashMoulton2(problem),
    odespy.AdamsBashMoulton3(problem),
    odespy.AdamsBashforth2(problem),
    odespy.AdamsBashforth3(problem),
    odespy.AdamsBashforth4(problem),
    odespy.AdaptiveResidual(problem, solver='Euler'),
    odespy.Backward2Step(problem),
    odespy.BackwardEuler(problem),
    odespy.Dop853(problem, rtol=rtol, atol=atol),
    odespy.Dopri5(problem, rtol=rtol, atol=atol),
    odespy.Euler(problem),
    odespy.Heun(problem),
    odespy.Leapfrog(problem),
    odespy.LeapfrogFiltered(problem),
    odespy.MidpointImplicit(problem),
    odespy.MidpointIter(problem, max_iter=10, eps_iter=1E-7),
    odespy.RK2(problem),
    odespy.RK3(problem),
    odespy.RK4(problem),
    odespy.RKFehlberg(problem, rtol=rtol, atol=atol),
    odespy.SymPy_odefun(problem),
    odespy.ThetaRule(problem),
    odespy.Trapezoidal(problem),
    odespy.Vode(problem, rtol=rtol, atol=atol,
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

for solver in solvers:
    solver_name = str(solver)
    results[solver_name] = {'dt': 0, 'error': numpy.nan, 'cpu': 0}

    solver.set_initial_condition([problem.Theta, 0])

    N = N_per_period*problem.period
    time_points = numpy.linspace(0, T, N+1)

    try:
        t0 = time.clock()
        u, t = solver.solve(time_points)
        cpu_time = time.clock() - t0
    except Exception, e:
        print solver_name, 'FAILED!'
        print e
        continue  # continue with next pass in the loop

    theta = u[:,0]
    error = numpy.abs(theta_exact(t) - theta)
    error_L2 = sqrt(numpy.sum(error**2)/N)
    #print solver_name, error_L2, cpu_time
    results[solver_name]['dt'] = t[1] - t[0]
    results[solver_name]['error'] = error_L2
    results[solver_name]['cpu'] = cpu_time

# Print
for solver_name in results:
    print '%-50s %.1E %.2f' % (solver_name,
                               results[solver_name]['error'],
                               results[solver_name]['cpu'])

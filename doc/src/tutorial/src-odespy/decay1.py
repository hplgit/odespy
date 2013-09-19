"""Simple example of decay ODE integrated until asymptotic value is reached."""
c = -0.1
A = 1.5

def f(u, t):
    return c*u

import odespy
solver = odespy.RK4(f)
solver.set_initial_condition(A)

import numpy
# Make sure integration interval [0, T] is large enough
N = 50
T = 150
time_points = numpy.linspace(0, T, N+1)

tol = 0.001  # tolerance for terminating the simulation

def terminate(u, t, step_no):
    # Most recent solution is in u[step_no] at time t[step_no]
    if abs(u[step_no]) <= tol:
        return True
    else:
        return False

u, t = solver.solve(time_points, terminate)

print "Solve u'=%g*u, u(0)=%g, for t in [%g, %g] and u>%g" % \
      (c, A, time_points[0], time_points[-1], tol)
print 'Final u(t=%g)=%g after %d steps' % (t[-1], u[-1], len(u)-1)

from matplotlib.pyplot import *
print plot
plot(t, u)
show()

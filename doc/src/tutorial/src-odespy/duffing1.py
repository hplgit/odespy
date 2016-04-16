import odespy, numpy
from numpy import sin, cos, pi, sqrt

c = 1
Theta = pi/4

def f(u, t):
    x, dxdt = u
    return [dxdt, cos(t)**3 - sin(t) - x - dxdt - x**3]

def u_exact(t):
    return cos(t)

solver = odespy.RK4(f)
solver.set_initial_condition([1, 0])
N_P = 120  # no of periods
T = 2*pi*N_P
Ns_per_period = [40, 80, 160, 320]
for N_per_period in Ns_per_period:
    N = N_per_period*N_P
    time_points = numpy.linspace(0, T, N+1)
    u, t = solver.solve(time_points)

    x = u[:,0]
    v = u[:,1]
    error = u_exact(t) - x
    dt = t[1] - t[0]
    print N_per_period, numpy.abs(error).max(), numpy.abs(error).max()/dt**4

import sys
sys.exit(0)
from matplotlib.pyplot import *
i = -10*N_per_period  # start index for the last 10 periods
plot(t[i:], error[i:], 'r-')
savefig('tmppng'); savefig('tmp.pdf')
show()


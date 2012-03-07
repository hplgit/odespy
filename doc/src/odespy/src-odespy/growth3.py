"""As exgr1.py, but c as argument in f."""
c = 0.1
A = 1.5
import odespy

# f has extra positional argument
def f(u, t, c):
    return c*u

solver = odespy.RK4(f, f_args=[c])

# Alternative: f has extra keyword argument
def f(u, t, c=1):
    return c*u

solver = odespy.RK4(f, f_kwargs={'c': c})
solver.set_initial_condition(A)

import numpy
N = 30
time_points = numpy.linspace(0, 40, N+1)
u, t = solver.solve(time_points)

from matplotlib.pyplot import *
plot(t, u)
show()

def f(u, t):
    return 1j*w*u

import odespy, numpy

w = 2*numpy.pi
solver = odespy.RK4(f, complex_valued=True)
solver.set_initial_condition(1+0j)
u, t = solver.solve(numpy.linspace(0, 6, 101))

from matplotlib.pyplot import *
plot(t, u.real)
show()


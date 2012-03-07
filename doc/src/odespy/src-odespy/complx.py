def f(u, t):
    return 1j*w*u

import odespy
import numpy as np
from matplotlib.pyplot import *

w = 2*np.pi
solver = odespy.RK4(f)
solver.set_initial_condition(1+0j)
u, t = solver.solve(np.linspace(0, 6, 101))
plot(t, u.real)
show()


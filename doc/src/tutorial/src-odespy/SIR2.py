"""SIR model for disease modeling."""

def f(u, t, beta, nu):
    S, I, R = u
    return [-beta(t)*S*I,
            beta(t)*S*I - nu*I,
            nu*I]

import odespy
import numpy as np
import sys

# beta = 0.0005; nu = 0.1; I0 = 1
# beta = 0.0001; nu = 0.1; I0 = 1 and 50
beta = lambda t: 0.0005 if t < 12 else 0.00005
nu = 0.1
I0 = 1
S0 = 1500.

solver = odespy.Euler(f, f_args=(beta, nu))
solver.set_initial_condition([S0, I0, 0])
dt = 0.5  # t counts days
T = 60
N = int(T/dt)
t = np.linspace(0, T, N+1)
u, t = solver.solve(t)
S, I, R = u[:,0], u[:,1], u[:,2]

from matplotlib.pyplot import *
plot(t, S, 'r-',
     t, I, 'b-',
     t, R, 'g-')
legend(['S', 'I', 'R'])
savefig('tmppng'); savefig('tmp.pdf')
show()


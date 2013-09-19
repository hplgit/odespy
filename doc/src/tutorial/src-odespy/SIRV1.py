"""SIRV model for disease modeling with vaccination."""

def f(u, t, beta, nu, p):
    S, I, R, V = u
    return [-beta*S*I - p(t)*S,
            beta*S*I - nu*I,
            nu*I,
            p(t)*S]

import odespy
import numpy as np
import sys

I0   = 1
S0 = 1500.
beta = 0.0005
nu = 0.1
p = lambda t: 0.1 if t >= 6 and t <= 15 else 0
#p = lambda t: 0.1

solver = odespy.Euler(f, f_args=(beta, nu, p))
solver.set_initial_condition([S0, I0, 0, 0])
dt = 0.5  # t counts days
T = 60
N = int(T/dt)
t = np.linspace(0, T, N+1)
u, t = solver.solve(t)
S, I, R, V = u[:,0], u[:,1], u[:,2], u[:,3]

from matplotlib.pyplot import *
plot(t, S, 'r-',
     t, I, 'b-',
     t, R, 'g-',
     t, V, 'y-')
legend(['S', 'I', 'R', 'V'])
savefig('tmppng'); savefig('tmp.pdf')
show()


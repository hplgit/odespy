"""As SIRV1.py, but seasonal change in beta."""

def f(u, t, beta, nu, p, mu, rho):
    S, I, R, V = u
    return [-beta(t)*S*I - p(t)*S + mu*R + rho*V,
            beta(t)*S*I - nu*I,
            nu*I - mu*R,
            p(t)*S - rho*V]

import odespy
import numpy as np
import sys

I0   = 1
S0 = 1500.
period = 360
beta = lambda t: 0.0003*np.sin(np.pi/period*t)**6
nu = 1./10
mu = 1./100
rho = 1./150
p = lambda t: 0.01*np.sin(np.pi/period*(t-50))**6
p = lambda t: 0.0

solver = odespy.Euler(f, f_args=(beta, nu, p, mu, rho))
solver.set_initial_condition([S0, I0, 0, 0])
dt = 0.5  # t counts days
T = 1400
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


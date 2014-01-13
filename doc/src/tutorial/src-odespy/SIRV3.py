"""As SIRV1.py, but computing max I(t) as a function of vaccination period."""

def f(u, t, beta, nu, p):
    S, I, R, V = u
    return [-beta*S*I - p(t)*S,
            beta*S*I - nu*I,
            nu*I,
            p(t)*S]

import odespy
import numpy as np
import sys, time

def simulate(beta, nu, p0, S0, I0, V_T):
    p = lambda t: p0 if 6 <= t <= 6+V_T else 0
    solver = odespy.Euler(f, f_args=(beta, nu, p))
    solver.set_initial_condition([S0, I0, 0, 0])
    dt = 0.5  # t counts days
    T = 60
    N = int(T/dt)
    t = np.linspace(0, T, N+1)
    u, t = solver.solve(t)
    I = u[:,1]
    return I.max()

I0   = 1
S0 = 1500.
beta = 0.0005
nu = 0.1
p0 = 0.1
V_Ts = np.linspace(0, 30, 31)
I_max = np.array([simulate(beta, nu, p0, S0, I0, V_T) for V_T in V_Ts])

tol = 5 # tolerance for comparing I_max[i] with final I_max[-1]
i = 0
print I_max
while I_max[i] - I_max[-1] > tol:
    i += 1
print 'No use in letting V_T > %d' % (V_Ts[i])

import matplotlib.pyplot as plt
plt.plot(V_Ts, I_max, 'b-')
plt.legend(['max I(t)'])
#plt.title('max I(t) as function of vaccination period')
plt.savefig('tmp.png'); plt.savefig('tmp.pdf')
plt.show()




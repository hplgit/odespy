# Author: Liwei Wang

"""
Chemical kinetics, from opkdmain.f in ODEPACK
du1/dt = -0.04*u1+1e4*u2*u3
du2/dt = 0.04*u1-1e4*u2*u3-3e7*u2*u2
0 = u1+u2+u3-1

This example is the typical usage of DAE with
user-supplied functions composed in Python.
"""
from odespy import *
import scitools.std as st
import numpy as np

def f(u, t):
    udot = np.zeros(3, float)
    udot[0] = -.04*u[0] + 1e4*u[1]*u[2]
    udot[1] = .04*u[0] - 1e4*u[1]*u[2] - 3e7*u[1]*u[1]
    udot[2] = sum(u) - 1
    return udot

def mas():
    m = [[1, 0, 0],
         [0, 1, 0],
         [0, 0, 0]]
    return m

def jac(u, t):
    return np.asarray(((-.04, 1e4*u[2], 1e4*u[1]),
                      (.04, -1e4*u[2] - 6e7*u[1], -1e4*u[1]),
                      (1.,1.,1.)))

import sys
try:
    n_points = int(sys.argv[1])    # Read from input
except:
    n_points = 10   # default number of time-steps

t0, tn, u0 = 0., 4.,  [1.,0.,0.]
time_points = np.linspace(t0, tn, n_points)
atol, rtol = [1e-6,1e-8,1e-6], 1e-4

st.figure()
method = Radau5Implicit
exact_final = [9.055142e-1, 2.240418e-5, 9.446344e-2]

# Test case 1: Radau5, with f, mas & jac
m = method(f=f, mas=mas, rtol=rtol, atol=atol, jac=jac)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'b-', title="Radau5 with Python functions",
        legend="with f, mas & jac", hold="on")
print 'Max error for test case 1 is %g' % max(u[-1] - exact_final)

# Test case 2: Radau5, with f, mas
m = method(f=f, mas=mas, rtol=rtol, atol=atol)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'r*', title="Radau5 with Python functions",
        legend="with f, mas", hold="on")
print 'Max error for test case 1 is %g' % max(u[-1] - exact_final)



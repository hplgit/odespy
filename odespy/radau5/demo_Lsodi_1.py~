# Author: Liwei Wang

"""
Chemical kinetics, from opkdmain.f in ODEPACK
du1/dt = -0.04*u1+1e4*u2*u3
du2/dt = 0.04*u1-1e4*u2*u3-3e7*u2*u2
0 = u1+u2+u3-1

This example is the typical usage of Lsodi with 
user-supplied functions composed in Python.
"""
from odespy import *
import scitools.std as st
import numpy as np

def res(u, t, s, ires):
    r = np.zeros(3,float)
    r[0] = -.04*u[0] + 1e4*u[1]*u[2] - s[0]
    r[1] = .04*u[0] - 1e4*u[1]*u[2] - 3e7*u[1]*u[1] - s[1]
    r[2] = sum(u) - 1
    return r,ires

def adda(u, t, p):
    p[0][0] += 1.
    p[1][1] += 1.
    return p

def jac(u, t, s):
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
ydoti = [-.04,.04,0.]      # initial value of du/dt
atol, rtol = [1e-6,1e-8,1e-6], 1e-4

st.figure()
method = Lsodi
exact_final = [9.055142e-1, 2.240418e-5, 9.446344e-2]

# Test case 1: Lsodi, with res, adda, ydoti & jac
m = method(res=res, rtol=rtol, atol=atol, ydoti=ydoti,
           adda_lsodi=adda, jac_lsodi=jac)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'b-', title="Lsodi with Python functions",
        legend="with res, adda, ydoti & jac", hold="on")
print 'Max error for test case 1 is %g' % max(u[-1] - exact_final)

# Test case 2: Lsodi, with res, ydoti & adda
m = method(res=res, rtol=rtol, atol=atol, ydoti=ydoti,
           adda_lsodi=adda)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'r*', title="Lsodi with Python functions",
        legend="with res, adda & ydoti", hold="on")
print 'Max error for test case 1 is %g' % max(u[-1] - exact_final)



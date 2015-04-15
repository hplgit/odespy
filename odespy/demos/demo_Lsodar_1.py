# Author: Liwei Wang

"""
Van Der Pol oscillator
u''=100*(1-u*u)*u'-u, with supplied full jacobian matrix
With one root function:  g = u

This example is the typical usage of Lsodar with 
user-supplied functions composed in Python.
"""

from odespy import *
#import scitools.basics,easyviz as st
import scitools.std as st
import numpy as np

def f(u,t):
    u00, u11 = u
    udot = [u11,100.*(1 - u00*u00)*u11 - u00]
    return udot

def jac(u,t):
    u00, u11 = u
    return [[0.,1.],
            [-200.*u00*u11 - 1.,100.*(1. - u00*u00)]]

def g(u,t):
    u00, u11 = u
    return u00

import sys
try:
    n_points = int(sys.argv[1])    # Read from input 
except:
    n_points = 20   # default number of time-steps    

t0, tn, u0 = 0., 200.,  [2.,0.]
time_points = np.linspace(t0, tn, n_points)
atol, rtol = 1e-4, 1e-4


st.figure()
method = Lsodar

# Test case 1: Lsodar, with f, g & jac
m = method(f, rtol=rtol, atol=atol, jac=jac, g=g)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], title="Lsodar with Python functions",
        legend="with f, g & jac", hold="on")

# Test case 2: Lsodar, with f & g
m = method(f, rtol=rtol, atol=atol, g=g)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Lsodar with Python functions", 
        legend="with f & g", hold="on")


# Author: Liwei Wang

"""
Test for dlsodar(), from demo in ODEPACK
u' = ((2*log(u) + 8)/t - 5)*u,  u(0)=1
With two root functions:
g1 = du/dt = ((2*log(u) + 8)/t - 5)*u,  root = 2.5
g2 = log(u)-2.2491,                     root = 2.47 and 2.53
Exact solution is u(t) = exp(-t**2 + 5*t - 4)

This example is the typical usage of Lsodar with 
user-supplied functions composed in Python.
"""

from odespy import *
#import scitools.basics,easyviz as st
import scitools.std as st
import numpy as np

def f(u,t):
    return ((2.*np.log(u)+8.)/t-5.)*u

def u_solution(t):
    return np.exp(-t ** 2 + 5 * t - 4)

def g(u,t):
    return [((2.*np.log(u) + 8.)/t - 5.)*u,\
            np.log(u) - 2.2491]

import sys
try:
    n_points = int(sys.argv[1])    # Read from input 
except:
    n_points = 200   # default number of time-steps    

t0, tn, u0 = 1., 6.,  1.0
time_points = np.linspace(t0, tn, n_points)
atol, rtol = 1e-8, 1e-8


st.figure()
method = Lsodar

# Test case: Lsodar, with f & g 
m = method(f, rtol=rtol, atol=atol, g=g)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u, 'o', 
        title="Lsodar with Python functions",
        legend="with f & g", hold="on")

# Exact solution:  u(t) = exp(-t**2+5*t-4)
st.plot(t, u_solution(time_points), '-', 
        title="Lsodar with Python functions", 
        legend="Exact solution", hold="on")

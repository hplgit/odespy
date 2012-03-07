# Author: Liwei Wang

"""
Van Der Pol oscillator
u''=3*(1-u*u)*u'-u, with supplied full jacobian matrix

This example is the typical usage of Lsode/Lsoda with 
user-supplied functions composed in Python.
"""

from odespy import *
#import scitools.basics,easyviz as st
import scitools.std as st
import numpy as np

def f(u,t):
    u00,u11 = u
    udot = [u11,3.*(1 - u00*u00)*u11 - u00]
    return udot

def jac(u,t):
    u00,u11 = u
    return [[0.,1.],
            [-6.*u00*u11 - 1.,3.*(1. - u00*u00)]]

import sys
try:
    n_points = int(sys.argv[1])    # Read from input 
except:
    n_points = 10   # default number of time-steps    


t0, tn, u0 = 0., 10.,  [2.,0.]
time_points = np.linspace(t0, tn, n_points)
atol, rtol = 1e-4, 1e-4


st.figure()
method = Lsode

# Test case 1: Lsode, with f & jac
m = method(f, rtol=rtol, atol=atol, jac=jac)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], title="Lsode with Python functions",
        legend="with f & jac", hold="on")

# Test case 2: Lsode, with f
m = method(f, rtol=rtol, atol=atol)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Lsode with Python functions", 
        legend="with f", hold="on")

method = Lsoda
st.figure()
# Test case 3: Lsoda, with f & jac
m = method(f, rtol=rtol, atol=atol, jac=jac)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], title="Lsoda with Python functions",
        legend="with f & jac", hold="on")

# Test case 4: Lsoda, with f
m = method(f, rtol=rtol, atol=atol)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Lsoda with Python functions", 
        legend="with f", hold="on")


st.figure()
# Test case 5: Switch to Lsodar
m.switch_to(Lsodar)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '-', title="Lsodar with Python functions", 
        legend="Switch from Lsoda", hold="on")

# Test case 6: Supplement jac
m.set(jac=jac, iter_method=None)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'o', title="Lsodar with Python functions", 
        legend="Supplement jac", hold="on")

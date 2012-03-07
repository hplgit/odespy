# Author: Liwei Wang

"""
u'= A*y, With supplied banded jacobian matrix

This example is the typical usage of Lsode/Lsoda with 
user-supplied functions composed in Python.
"""


from odespy import *
#import scitools.basics,easyviz as st
import scitools.std as st
import numpy as np

def f(u,t):
    udot = np.zeros(25,float)
    for j in range(5):
        for i in range(5):
            k = i+j*5
            udot[k] = -2.*u[k] + u[k - 1]*(i>0) + u[k - 5]*(j>0)
    return (np.asarray(udot))

def jac_banded(u,t,ml,mu):
    pd = np.zeros(6*25,float).reshape(6,25)
    pd[0,:], pd[1,:], pd[5,:], pd[1,4:24] = -2., 1., 1., 0.
    return pd

import sys
try:
    n_points = int(sys.argv[1])    # Read from input 
except:
    n_points = 10   # default number of time-steps    


t0, tn, u0 = 0., 4.,  [1]+24*[0]
ml, mu = 5, 0
time_points = np.linspace(t0, tn, n_points)
atol, rtol = 1e-2, 1e-2


st.figure()
method = Lsode

# Test case 1: Lsode, with f, ml, mu & jac_banded
m = method(f, rtol=rtol, atol=atol, 
           ml=ml, mu=mu, jac_banded=jac_banded)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], title="Lsode with Python functions",
        legend="with f, ml, mu & jac", hold="on")

# Test case 2: Lsode, with f, ml, mu
m = method(f, rtol=rtol, atol=atol, 
           ml=ml, mu=mu)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Lsode with Python functions", 
        legend="with f & jac", hold="on")

# Test case 3: Lsode, with f
m = method(f, rtol=rtol, atol=atol)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'o', title="Lsode with Python functions", 
        legend="with f", hold="on")

method = Lsoda
st.figure()
# Test case 4: Lsoda, with f, ml, mu & jac_banded
m = method(f, rtol=rtol, atol=atol, 
           ml=ml, mu=mu, jac_banded=jac_banded)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], title="Lsoda with Python functions",
        legend="with f, ml, mu & jac", hold="on")

# Test case 5: Lsoda, with f, ml, mu
m = method(f, rtol=rtol, atol=atol, 
           ml=ml, mu=mu)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Lsoda with Python functions", 
        legend="with f & jac", hold="on")

# Test case 6: Lsoda, with f
m = method(f, rtol=rtol, atol=atol)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'o', title="Lsoda with Python functions", 
        legend="with f", hold="on")

st.figure()
# Test case 7: Switch to Lsodar
m.switch_to(Lsodar)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '-', title="Lsodar with Python functions", 
        legend="Switch from Lsoda", hold="on")

# Test case 8: Supplement ml, mu
m.set(ml=ml, mu=mu)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Lsodar with Python functions", 
        legend="Supplement ml,mu", hold="on")

# Test case 9: Supplement jac_banded
m.set(jac_banded=jac_banded, iter_method=None)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'o', title="Lsodar with Python functions", 
        legend="Supplement jac_banded", hold="on")

# Author: Liwei Wang

"""
u'= A*y, With supplied full or banded jacobian matrix

This example is the typical usage of Lsode/Lsoda with 
user-supplied functions composed in Python.
"""


from odespy import *
#import scitools.basics,easyviz as st
import scitools.std as st
import numpy as np

def f(u,t):
    udot = np.zeros(5,float)
    udot[0] = 0.1*u[0]-0.2*u[1]
    udot[1] = -0.3*u[0]+0.1*u[1]-0.2*u[2]
    udot[2] = -0.3*u[1]+0.1*u[2]-0.2*u[3]
    udot[3] = -0.3*u[2]+0.1*u[3]-0.2*u[4]
    udot[4] = -0.3*u[3]+0.1*u[4]
    return (np.asarray(udot))

def jac_full(u,t):    # full matrix
    pd = np.asarray([[.1,-.2,0.,0.,0.],
                   [-.3,.1,-.2,0.,0.],
                   [0.,-.3,.1,-.2,0.],
                   [0.,0.,-.3,.1,-.2],
                   [0.,0.,0.,-.3,.1]])
    return pd

def jac_banded(u,t,ml,mu):    # banded matrix
    pd = np.array([0.,-.2,-.2,-.2,-.2,
                   0.1,.1,.1,.1,.1,
                   -.3,-.3,-.3,-.3,0.]).reshape(3,5)
    return pd


import sys
try:
    n_points = int(sys.argv[1])    # Read from input 
except:
    n_points = 10   # default number of time-steps    


t0, tn, u0 = 0., 10., [1,2,3,4,5]
ml, mu = 1, 1
time_points = np.linspace(t0, tn, n_points)
atol, rtol = 1e-4, 1e-4


st.figure()
method = Lsode

# Test case 1: Lsode, with f, ml, mu & jac_banded
m = method(f, rtol=rtol, atol=atol, 
           ml=ml, mu=mu, jac_banded=jac_banded)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], title="Lsode with Python functions",
        legend="with f, ml, mu & jac_banded", hold="on")

# Test case 2: Lsode, with f, ml, mu
m = method(f, rtol=rtol, atol=atol, 
           ml=ml, mu=mu)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Lsode with Python functions", 
        legend="with f & jac_banded", hold="on")

# Test case 3: Lsode, with f & jac_full
m = method(f, rtol=rtol, atol=atol, jac=jac_full)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'o', title="Lsode with Python functions", 
        legend="with f & jac_full", hold="on")

# Test case 4: Lsode, with f
m = method(f, rtol=rtol, atol=atol)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'o', title="Lsode with Python functions", 
        legend="with f", hold="on")

method = Lsoda
st.figure()
# Test case 5: Lsoda, with f, ml, mu & jac_banded
m = method(f, rtol=rtol, atol=atol, 
           ml=ml, mu=mu, jac_banded=jac_banded)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], title="Lsoda with Python functions",
        legend="with f, ml, mu & jac_banded", hold="on")

# Test case 6: Lsoda, with f, ml, mu
m = method(f, rtol=rtol, atol=atol, 
           ml=ml, mu=mu)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Lsoda with Python functions", 
        legend="with f & jac_banded", hold="on")

# Test case 7: Lsoda, with f & jac_full
m = method(f, rtol=rtol, atol=atol, jac=jac_full)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'o', title="Lsoda with Python functions", 
        legend="with f & jac_full", hold="on")

# Test case 8: Lsoda, with f
m = method(f, rtol=rtol, atol=atol)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'o', title="Lsoda with Python functions", 
        legend="with f", hold="on")


st.figure()
# Test case 9: Switch to Lsodar
m.switch_to(Lsodar)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '-', title="Lsodar with Python functions", 
        legend="Switch from Lsoda", hold="on")

# Test case 10: Supplement jac_full
m.set(jac=jac_full)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Lsodar with Python functions", 
        legend="Supplement jac_full", hold="on")

# Test case 11: Supplement ml, mu, Remove jac_full
m.set(ml=ml, mu=mu, iter_method=None, jac=None)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Lsodar with Python functions", 
        legend="Supplement ml, mu", hold="on")

# Test case 12: Supplement jac_banded
m.set(jac_banded=jac_banded, iter_method=None)

u,t = m.solve(time_points)
st.plot(t, u[:,0], 'o', title="Lsodar with Python functions", 
        legend="Supplement jac_banded", hold="on")


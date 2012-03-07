# Author: Liwei Wang

"""
Test for Lsodes, from demo of ODEPACK
u'= A*u, where A is assumed to be sparse matrix

This example is the typical usage of Lsodes with 
user-supplied functions composed in Python.
"""

from odespy import *
#import scitools.basics,easyviz as st
import scitools.std as st
import numpy as np

def f(u,t):    # A * u
    u_matrix = np.mat(u).T
    A = np.mat([[-4.,1.,0.,1.,0.,0.,0.,0.,0.],
                [1.,-4.,1.,0.,1.,0.,0.,0.,0.],
                [0.,1.,-4.,0.,0.,1.,0.,0.,0.],
                [0.,0.,0.,-4.,1.,0.,1.,0.,0.],
                [0.,0.,0.,1.,-4.,1.,0.,1.,0.],
                [0.,0.,0.,0.,1.,-4.,0.,0.,1.],
                [0.,0.,0.,0.,0.,0.,-4.,1.,0.],
                [0.,0.,0.,0.,0.,0.,1.,-4.,1.],
                [0.,0.,0.,0.,0.,0.,0.,1.,-4.]])
    return np.array((A*(np.mat(u).T)).T).reshape(9,)

def jac_column(u,t,j,ia,ja):
    jac = np.array([[-4.,1.,0.,0.,0.,0.,0.,0.,0.],
                    [1.,-4.,1.,0.,0.,0.,0.,0.,0.],
                    [0.,1.,-4.,0.,0.,0.,0.,0.,0.],
                    [1.,0.,0.,-4.,1.,0.,0.,0.,0.],
                    [0.,1.,0.,1.,-4.,1.,0.,0.,0.],
                    [0.,0.,1.,0.,1.,-4.,0.,0.,0.],
                    [0.,0.,0.,1.,0.,0.,-4.,1.,0.],
                    [0.,0.,0.,0.,1.,0.,1.,-4.,1.],
                    [0.,0.,0.,0.,0.,1.,0.,1.,-4.]])
    return jac[j]

import sys
try:
    n_points = int(sys.argv[1])    # Read from input 
except:
    n_points = 20   # default number of time-steps    

t0, tn, u0 = 0., 3., np.arange(1.,10.,1.)
time_points = np.linspace(t0, tn, n_points)

atol, rtol = 1e-4, 1e-4
ia = [1,3,6,8,11,15,18,21,25,28]   # Describe the sparse structure
ja = [1,2,1,2,3,2,3,1,4,5,2,4,5,6,3,5,6,4,7,8,5,7,8,9,6,8,9]

st.figure()
method = Lsodes

# Test case 1: Lsodes, with f, ia, ja & jac_column
m = method(f, rtol=rtol, atol=atol, 
           ia=ia, ja=ja, jac_column=jac_column)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], title="Lsodes with Python functions",
        legend="with f, ia, ja & jac_column", hold="on")

# Test case 2: Lsodes, with f, ia, ja 
m = method(f, rtol=rtol, atol=atol, 
           ia=ia, ja=ja)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Lsodes with Python functions", 
        legend="with f, ia & ja", hold="on")

# Test case 3: Lsodes, with f, jac_column
m = method(f, rtol=rtol, atol=atol, jac_column=jac_column, lrw=400)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Lsode with Python functions", 
        legend="with f & jac_column", hold="on")

# Test case 4: Lsodes, with f
m = method(f, rtol=rtol, atol=atol)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'o', title="Lsodes with Python functions", 
        legend="with f", hold="on")


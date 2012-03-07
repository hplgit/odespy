# Author: Liwei Wang

"""
  Demo problem from netlib repository.

  This is a simplification of Example 1 of B.P. Sommeijer, L.F. Shampine, 
  and J.G. Verwer, RKC: an Explicit Solver for Parabolic PDEs that shows
  the use of RKC on a substantial problem.  Semi-discretization of the 
  heat equation in three space variables results in 19**3 = 6859 equations.

  In this test example, F is provide as a Python function. This lead to 
  unavoidable efficiency cost to call-back F from extension module. For users 
  who are familiar with Fortran language, there is a better way to provide F 
  in Fortran language instead of Python. See demo_RKC_2_fortran.py.
  
"""

from odespy import *
import scitools.std as st
import numpy as np

def f(u,t):
    uwork = np.zeros((21, 21, 21), float)   
    udot = np.zeros(6859, float)
    for i in range(1,20):
        for j in range(1,20):
            for k in range(1,20):
                uwork[i][j][k] = u[i - 1 + (j - 1)*19 + (k - 1)*19*19]
    for i in range(1,20):
        for j in range(1,20):
            uwork[i][j][0] = np.tanh(5.0*(i*.05 + j*.1 - .5 - t))
            uwork[i][j][20] = np.tanh(5.0*(i*.05 + j*.1 + 1. - t))
    for i in range(1,20):
        for k in range(1,20):
            uwork[i][0][k] = np.tanh(5.0*(i*.05 + k*.075 - .5 - t))
            uwork[i][20][k] = np.tanh(5.0*(i*.05 + k*.075 + 1.5 - t))
    for j in range(1,20):
        for k in range(1,20):
            uwork[0][j][k] = np.tanh(5.*(j*.1 + k*.075 - .5 - t))
            uwork[20][j][k] = np.tanh(5.*(j*.1 + k*.075 + .5 - t))
    for i in range(1,20):
        for j in range(1,20):
            for k in range(1,20):
                arg = 5.*(i*.05 + .1*j + .075*k - .5 - t)
                sh, ch = np.sinh(arg), np.cosh(arg)
                l = i - 1 + (j-1)*19 + (k-1)*19*19
                udot[l] = (uwork[i-1,j,k]-2.*uwork[i,j,k]+uwork[i+1,j,k])*400 +\
                          (uwork[i,j-1,k]-2.*uwork[i,j,k]+uwork[i,j+1,k])*400 +\
                          (uwork[i,j,k-1]-2.*uwork[i,j,k]+uwork[i,j,k+1])*400 +\
                          (-5.*ch + 362.5*sh)/(ch**3)
    return udot

def spcrad(u,t):
    return 4800.0

import sys
try:
    n_points = int(sys.argv[1])    # Read from input 
except:
    n_points = 20   # default number of time-steps    

t0, tn = 0., .7 
u0 = np.zeros(6859,float)
for i in range(19):
    for j in range(19):
        for k in range(19):
            index = i + j*19 + k*19*19
            u0[index] = np.tanh(5.0*((i+1)*.05 + (j+1)*.1 + \
                                (k+1)*.075 - .5))

time_points = np.linspace(t0, tn, n_points)
atol, rtol = 1e-2, 1e-2
jac_constant = 1

st.figure()
method = RKC
print 'This test will take several minutes, please wait...'

# Test case 1: RKC with f, jac_constant, spcrad
m = method(f, rtol=rtol, atol=atol, spcrad=spcrad, jac_constant=jac_constant)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], title="RKC with Python functions",
        legend="with f, jac_constant, spcrad", hold="on")

# Test case 2: RKC with f, jac_constant
m = method(f, rtol=rtol, atol=atol, jac_constant=jac_constant)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="RKC with Python functions",
        legend="with f, jac_constant", hold="on")

# Test case 3: RKC with f, spcrad
m = method(f, rtol=rtol, atol=atol, spcrad=spcrad)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'go', title="RKC with Python functions",
        legend="with f, spcrad", hold="on")


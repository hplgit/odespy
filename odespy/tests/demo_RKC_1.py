# Author: Liwei Wang

"""
  This example shows how to use RKC.  It solves a system of ODEs that
  arise from semi-discretization of the reaction-diffusion equation

            U  = U   + (1 - U)*U**2    for   t >= 0,   0 <= x <= 10
             t    xx

  Dirichlet boundary conditions specify U(0,t) and U(10,t) for all t >= 0
  and the initial values U(x,0) are specified. 

  These values are taken from an analytical solution that is evaluated in 
  sol(x,t) so that the numerical solution can be compared to a known solution. 

  In this test example, F is provide as a Python function. This lead to 
  unavoidable efficiency cost to call-back F from extension module. For users 
  who are familiar with Fortran language, there is a better way to provide F 
  in Fortran language instead of Python. See demo_RKC_1_fortran.py.
  
"""


from odespy import *
import scitools.std as st
import numpy as np

def f(u,t):
    udot = np.zeros(99, float)
    sol0 = 1.0/(1.0 + np.exp(np.sqrt(.5)*(0. - np.sqrt(.5)*t)))
    udot[0] = (sol0 - 2.*u[0] + u[1])/.01 + (1. - u[0])*u[0]**2
    for i in range(1,98):
        udot[i] = (u[i - 1] - 2.*u[i] + u[i + 1])/.01 + (1. - u[i])*u[i]**2
    sol10 = 1.0/(1.0 + np.exp(np.sqrt(.5)*(10. - np.sqrt(.5)*t)))
    udot[98] = (u[97] - 2.*u[98] + sol10)/.01 + (1. - u[98])*u[98]**2
    return udot

import sys
try:
    n_points = int(sys.argv[1])    # Read from input 
except:
    n_points = 20   # default number of time-steps    

t0, tn = 0., 15. 
u0 = [1.0/(1.0 + np.exp(np.sqrt(.5)*i*0.1)) for i in range(1,100)]
time_points = np.linspace(t0, tn, n_points)



atol, rtol = 1e-4, 1e-4
st.figure()
# Test case 1: RKC
m = RKC(f, rtol=rtol, atol=atol)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], title="Python function f",
        legend="RKC", hold="on")

# Test case 2: Lsode 
m = m.switch_to(Lsode)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Python function f",
        legend="Lsode", hold="on")

# Test case 3: RKFehlberg
m = m.switch_to(RKFehlberg)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Python function f", 
        legend="RKFehlberg", hold="on")



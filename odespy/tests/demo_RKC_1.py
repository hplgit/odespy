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

def f(u, t, udot=None):
    n = u.size
    if udot is None:
        udot = np.zeros(n)

    sol0 = 1.0/(1.0 + np.exp(np.sqrt(.5)*(0. - np.sqrt(.5)*t)))
    udot[0] = (sol0 - 2.*u[0] + u[1])/.01 + (1. - u[0])*u[0]**2
    for i in range(1,n):
        udot[i] = (u[i - 1] - 2.*u[i] + u[i + 1])/.01 + (1. - u[i])*u[i]**2
    sol10 = 1.0/(1.0 + np.exp(np.sqrt(.5)*(10. - np.sqrt(.5)*t)))
    udot[-1] = (u[-2] - 2.*u[-1] + sol10)/.01 + (1. - u[-1])*u[-1]**2
    return udot

import sys
try:
    n_points = int(sys.argv[1])    # Read from input
except:
    n_points = 20   # default number of time-steps

t0, tn = 0., 15.
n = 9  # no of points in x grid
x = np.linspace(0, 10, n)
dx = x[1] - x[0]
u0 = [1.0/(1.0 + np.exp(np.sqrt(.5)*dx)) for i in range(1,n+1)]
time_points = np.linspace(t0, tn, n_points)



atol, rtol = 1e-4, 1e-4
st.figure()

# Needs to debug RKFehlberg for this particular example,
# replace 99 by a variable
raise Exception('Reimplement this test for variable grid size and debug RKFehlberg then.')
m = RKFehlberg(f, rtol=rtol, atol=atol)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'o', title="Python function f",
        legend="RKFehlberg", hold="on")
sys.exit(1)

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
st.plot(t, u[:,0], 'o', title="Python function f",
        legend="RKFehlberg", hold="on")




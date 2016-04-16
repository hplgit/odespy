# Author: Liwei Wang

"""
Same example as in demo_Lsodar_2.py
u' = ((2*log(u) + 8)/t - 5)*u,  u(0)=1
With two root functions:
g1 = du/dt = ((2*log(u) + 8)/t - 5)*u,  root = 2.5
g2 = log(u)-2.2491,                     root = 2.47 and 2.53
Exact solution is u(t) = exp(-t**2 + 5*t - 4)

This example is the typical usage of Lsodar with 
user-supplied functions composed in Python.
"""

from odespy import *
import scitools.std as st
import numpy as np
import os

def u_solution(t):
    return np.exp(-t ** 2 + 5 * t - 4)

f_str = """
      subroutine f_f77(neq, t, u, udot)
Cf2py intent(hide) neq
Cf2py intent(out) udot
      integer neq
      double precision t, u, udot
      dimension u(neq), udot(neq)
      udot(1) = ((2.0d0*log(u(1)) + 8.0d0)/t - 5.0d0)*u(1)
      return
      end
"""

g_str = """
      subroutine g_f77(neq, t, u, ng, groot)
Cf2py intent(hide) neq
Cf2py optional, intent(hide) ng
Cf2py intent(in) t, u
Cf2py intent(out) groot
      integer neq, ng
      double precision t, u, groot
      dimension u(neq), groot(ng)
      groot(1) = ((2.0d0*log(u(1)) + 8.0d0)/t - 5.0d0)*u(1)
      groot(2) = log(u(1)) - 2.2491d0
      return
      end
"""

import sys
try:
    n_points = int(sys.argv[1])    # Read from input 
except:
    n_points = 200   # default number of time-steps    

t0, tn, u0 = 1., 6.,  1.0
time_points = np.linspace(t0, tn, n_points)
atol, rtol, ng = 1e-8, 1e-8, 2
ng = 2     # number of constraint functions in g

# Compile these Fortran subroutines
from numpy import f2py
f2py.compile(f_str+'\n'+g_str, modulename='callback', verbose=False)   
import callback
f_f77, g_f77 = callback.f_f77, callback.g_f77

st.figure()
method = Lsodar
# Test case: Lsodar, with f & g 
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol, g_f77=g_f77, ng=ng)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u, 'o', 
        title="Lsodar with Fortran subroutines",
        legend="with f & g", hold="on")

# Exact solution:  u(t) = exp(-t**2+5*t-4)
st.plot(t, u_solution(time_points), '-', 
        title="Lsodar with Fortran subroutines", 
        legend="Exact solution", hold="on")

os.remove('callback.so')

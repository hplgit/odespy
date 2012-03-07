# Author: Liwei Wang

"""
Same example as in demo_Lsode_Lsoda_1.py:
Van Der Pol oscillator
u''=3*(1-u*u)*u'-u, with supplied full jacobian matrix

This example intends to show users how to improve efficiency with 
call-back functions composed in Fortran language. 
"""
from odespy import *
import scitools.std as st
import numpy as np
import os

f_str = """
      subroutine f_f77(neq, t, u, udot)
Cf2py intent(hide) neq
Cf2py intent(out) udot
      integer neq
      double precision t, u, udot,i,j,k
      dimension u(neq), udot(neq)
      udot(1) = u(2)
      udot(2) = 3d0*(1d0 - u(1)**2)*u(2) - u(1)
      return
      end
"""

jac_full_str = """
      subroutine jac_f77(neq, t, u, ml, mu, pd, nrowpd)
Cf2py intent(hide) neq, ml, mu, nrowpd
Cf2py intent(out) pd
      integer neq, ml, mu, nrowpd
      double precision t, u, pd
      dimension u(neq), pd(nrowpd,neq)
      pd(1,1) = 0d0
      pd(1,2) = 1d0
      pd(2,1) = -6d0*u(1)*u(2) - 1d0
      pd(2,2) = 3d0*(1d0-u(1)**2)
      return
      end
"""

import sys
try:
    n_points = int(sys.argv[1])    # Read from input 
except:
    n_points = 10   # default number of time-steps    

t0, tn, u0 = 0., 10.,  [2.,0.]
time_points = np.linspace(t0, tn, n_points)
atol, rtol = 1e-4, 1e-4

# Compile these Fortran subroutines
string_to_compile = '\n'.join([f_str, jac_full_str])
from numpy import f2py
f2py.compile(string_to_compile, modulename='callback', verbose=False)   
import callback
f_f77, jac_f77 = callback.f_f77, callback.jac_f77


st.figure()
method = Lsode

# Test case 1: Lsode, with f & jac
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol, jac_f77=jac_f77)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '-', title="Lsode with Fortran subroutines", 
        legend="with f & jac_full", hold="on")

# Test case 2: Lsode, with f
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'o', title="Lsode with Fortran subroutines", 
        legend="with f", hold="on")

method = Lsoda
st.figure()
# Test case 3: Lsoda, with f & jac_full
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol, jac_f77=jac_f77)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '-', title="Lsoda with Fortran subroutines", 
        legend="with f & jac_full", hold="on")

# Test case 4: Lsoda, with f
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Lsoda with Fortran subroutines", 
        legend="with f", hold="on")

st.figure()
# Test case 5: Switch to Lsodar
m.switch_to(Lsodar)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '-', title="Lsodar with Fortran subroutines", 
        legend="Switch from Lsoda", hold="on")

# Test case 6: Supplement jac_full
m.set(jac_f77=jac_f77)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Lsodar with Fortran subroutines", 
        legend="Supplement jac_full", hold="on")

os.remove('callback.so')

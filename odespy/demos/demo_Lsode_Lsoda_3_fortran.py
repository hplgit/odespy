# Author: Liwei Wang

"""
Same example as in demo_Lsode_Lsoda_3.py:
u'= A*y, With supplied full or banded jacobian matrix

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
      double precision t, u, udot
      dimension u(neq), udot(neq)
      udot(1) = 0.1*u(1)-0.2*u(2)
      udot(2) = -0.3*u(1)+0.1*u(2)-0.2*u(3)
      udot(3) = -0.3*u(2)+0.1*u(3)-0.2*u(4)
      udot(4) = -0.3*u(3)+0.1*u(4)-0.2*u(5)
      udot(5) = -0.3*u(4)+0.1*u(5)
      return
      end
"""

# banded Jacobian matrix
jac_banded_str = """
      subroutine jac_banded_f77(neq, t, u, ml, mu, pd, nrowpd)
Cf2py intent(hide) neq, ml, mu, nrowpd
Cf2py intent(out) pd
      integer neq, ml, mu, nrowpd
      double precision t, u, pd
      dimension u(neq), pd(nrowpd,neq)
      integer i
      do 10 i = 1, 5
        pd(1,i) = -2d-1
        pd(2,i) = 1d0
 10     pd(3,i) = -3d-1
      pd(1,1) = 0d0
      pd(2,1) = 1d-1
      pd(3,5) = 0d0
      return
      end
"""

# full Jacobian matrix
jac_full_str = """
      subroutine jac_f77(neq, t, u, ml, mu, pd, nrowpd)
Cf2py intent(hide) neq, ml, mu, nrowpd
Cf2py intent(out) pd
      integer neq, ml, mu, nrowpd
      double precision t, u, pd
      dimension u(neq), pd(neq,neq)
      integer i
      pd(1,1) = 1d-1
      do 10 i = 2, 5
        pd(i,i) = 1d-1
        pd(i,i-1) = -3d-1
 10     pd(i-1,i) = -2d-1
      return
      end
"""

import sys
try:
    n_points = int(sys.argv[1])    # Read from input 
except:
    n_points = 10   # default number of time-steps    


t0, tn, u0 = 0., 10., [1,2,3,4,5]
ml, mu = 1, 1
time_points = np.linspace(t0, tn, n_points)
atol, rtol = 1e-4, 1e-4

# Compile these Fortran subroutines
string_to_compile = '\n'.join([f_str, jac_banded_str, jac_full_str])
from numpy import f2py
f2py.compile(string_to_compile, modulename='callback', verbose=False)   
import callback
f_f77, jac_banded_f77, jac_f77 = callback.f_f77, callback.jac_banded_f77, \
    callback.jac_f77

st.figure()
method = Lsode

# Test case 1: Lsode, with f, ml, mu & jac_banded
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol, 
           ml=ml, mu=mu, jac_banded_f77=jac_banded_f77)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], title="Lsode with Fortran subroutines",
        legend="with f, ml, mu & jac", hold="on")

# Test case 2: Lsode, with f, ml, mu
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol, 
           ml=ml, mu=mu)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Lsode with Fortran subroutines", 
        legend="with f & jac_banded", hold="on")

# Test case 3: Lsode, with f & jac_full
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol, jac_f77=jac_f77)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'o', title="Lsode with Fortran subroutines", 
        legend="with f & jac_full", hold="on")

# Test case 4: Lsode, with f
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'o', title="Lsode with Fortran subroutines", 
        legend="with f", hold="on")

method = Lsoda
st.figure()
# Test case 5: Lsoda, with f, ml, mu & jac_banded
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol, 
           ml=ml, mu=mu, jac_banded_f77=jac_banded_f77)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], title="Lsoda with Fortran subroutines",
        legend="with f, ml, mu & jac_banded", hold="on")

# Test case 6: Lsoda, with f, ml, mu
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol, 
           ml=ml, mu=mu)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Lsoda with Fortran subroutines", 
        legend="with f & jac_banded", hold="on")

# Test case 7: Lsoda, with f & jac_full
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol, jac_f77=jac_f77)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'o', title="Lsoda with Fortran subroutines", 
        legend="with f & jac_full", hold="on")

# Test case 8: Lsoda, with f
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'o', title="Lsoda with Fortran subroutines", 
        legend="with f", hold="on")


st.figure()
# Test case 9: Switch to Lsodar
m.switch_to(Lsodar)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '-', title="Lsodar with Fortran subroutines", 
        legend="Switch from Lsoda", hold="on")

# Test case 10: Supplement jac_full
m.set(jac_f77=jac_f77)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Lsodar with Fortran subroutines", 
        legend="Supplement jac_full", hold="on")

# Test case 11: Supplement ml, mu, Remove jac_full
m.set(ml=ml, mu=mu, iter_method=None, jac=None)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Lsodar with Fortran subroutines", 
        legend="Supplement ml, mu", hold="on")

# Test case 12: Supplement jac_banded as string of Fortran code
m.set(jac_banded_f77=jac_banded_str)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'o', title="Lsodar with Fortran subroutines", 
        legend="Supplement jac_banded", hold="on")


os.remove('callback.so')



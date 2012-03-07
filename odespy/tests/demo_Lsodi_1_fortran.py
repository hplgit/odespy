# Author: Liwei Wang

"""
Sample example as in demo_Lsodi_1.py:
Simplified Galerkin Solution of Burgers Equation

This example is the typical usage of Lsodi with 
user-supplied functions composed in Fortran code.

"""
from odespy import *
import scitools.std as st
import numpy as np
import os

res_str = """
      subroutine res_f77(neq, t, u, s, r, ires)
Cf2py intent(hide) neq
Cf2py intent(out) r
Cf2py intent(in, out) ires
      double precision t, u, s, r
      dimension u(neq), s(neq), r(neq)
      r(1) = -.04*u(1) + 1.D4*u(2)*u(3) - s(1)
      r(2) = .04*u(1) - 1.D4*u(2)*u(3) - 3.D7*u(2)*u(2) - s(2)
      r(3) = u(1) + u(2) + u(3) - 1.
      return
      end
"""

adda_str = """
      subroutine adda_lsodi_f77(neq, t, u, ml, mu, pd, nrowpd)
Cf2py intent(in, hide) neq, ml, mu
Cf2py intent(in, hide), depend(pd) nrowpd
Cf2py intent(in, out) pd
      double precision t, u, pd
      dimension u(neq), pd(nrowpd, neq)
      pd(1,1) = pd(1,1) + 1.
      pd(2,2) = pd(2,2) + 1.
      return
      end
"""

jac_str = """
      subroutine jac_lsodi_f77(neq, t, u, s, ml, mu, pd, nrowpd)
Cf2py intent(in, hide) neq, ml, mu, nrowpd
Cf2py intent(out) pd
      integer neq, ml, mu, nrowpd
      double precision t, u, pd, s
      dimension u(neq), s(neq), pd(nrowpd, neq)
      pd(1,1) = -.04
      pd(1,2) = 1.D4*u(3)
      pd(1,3) = 1.D4*u(2)
      pd(2,1) = .04
      pd(2,2) = -1.D4*u(3) - 6.D7*u(2)
      pd(2,3) = -1.D4*u(2)
      pd(3,1) = 1.
      pd(3,2) = 1.
      pd(3,3) = 1.
      return
      end
"""
# compile these Fortran subroutines
string_to_compile = '\n'.join([res_str, adda_str, jac_str])
from numpy import f2py
f2py.compile(string_to_compile, modulename = "callback", verbose=False)
import callback
res_f77, adda_f77, jac_f77 = callback.res_f77, callback.adda_lsodi_f77, \
    callback.jac_lsodi_f77

t0, tn, u0, n_points = 0., 4.,  [1.,0.,0.], 10
time_points = np.linspace(t0, tn, n_points)
ydoti = [-.04,.04,0.]      # initial value of du/dt
atol, rtol = [1e-6,1e-8,1e-6], 1e-4

st.figure()
method = Lsodi
exact_final = [9.055142e-1, 2.240418e-5, 9.446344e-2]

# Test case 1: Lsodi, with res, adda, ydoti & jac
m = method(res_f77=res_f77, rtol=rtol, atol=atol, ydoti=ydoti,
           adda_lsodi_f77=adda_f77, jac_lsodi_f77=jac_f77)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'r-', title="Lsodi with Fortran subroutines",
        legend="with res, adda, ydoti & jac", hold="on")
print 'Max error for test case 1 is %g' % max(u[-1] - exact_final)

# Test case 2: Lsodi, with res, ydoti & adda
m = method(res_f77=res_f77, rtol=rtol, atol=atol, ydoti=ydoti,
           adda_lsodi_f77=adda_f77)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'g*', title="Lsodi with Fortran subroutines",
        legend="with res, adda & ydoti", hold="on")
print 'Max error for test case 2 is %g' % max(u[-1] - exact_final)

os.remove('callback.so')

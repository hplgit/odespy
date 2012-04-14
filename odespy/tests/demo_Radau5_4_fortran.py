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

f_str = """
      subroutine f_f77(neq, t, u, udot)
Cf2py intent(hide) neq
Cf2py intent(out) udot
      integer neq
      double precision t, u, udot
      dimension u(neq), udot(neq)
      udot(1) = -.04*u(1) + 1.D4*u(2)*u(3)
      udot(2) = .04*u(1) - 1.D4*u(2)*u(3) - 3.D7*u(2)*u(2)
      udot(3) = u(1) + u(2) + u(3) - 1.0D0
      return
      end

"""

jac_str = """
      subroutine jac_f77_radau5(neq,t,u,dfu,ldfu,rpar,ipar)
Cf2py intent(hide) neq,rpar,ipar
Cf2py intent(in)   t,u,ldfu
Cf2py intent(out)  dfu
      integer neq,ipar,ldfu
      double precision t,u,dfu,rpar
      dimension u(neq),dfu(ldfu,neq)
      dfu(1,1) = -.04
      dfu(1,2) = 1.D4*u(3)
      dfu(1,3) = 1.D4*u(2)
      dfu(2,1) = .04
      dfu(2,2) = -1.D4*u(3) - 6.D7*u(2)
      dfu(2,3) = -1.D4*u(2)
      dfu(3,1) = 1.
      dfu(3,2) = 1.
      dfu(3,3) = 1.
      return
      end
"""

mas_str = """
      subroutine mas_f77(neq,mas,lmas,rpar,ipar)
Cf2py intent(hide)   neq,rpar,ipar
Cf2py intent(in)     neq,lmas
Cf2py intent(out)    mas
      integer neq,lmas,ipar
      double precision mas(lmas,neq),rpar
      mas(1,1) = 1.
      mas(1,2) = 0.
      mas(1,3) = 0.
      mas(2,1) = 0.
      mas(2,2) = 1.
      mas(2,3) = 0.
      mas(3,1) = 0.
      mas(3,2) = 0.
      mas(3,3) = 0.
      return
      end
"""

# compile these Fortran subroutines
string_to_compile = '\n'.join([f_str, jac_str])
from numpy import f2py
f2py.compile(string_to_compile, modulename = "tmp_callback", verbose=False)
import tmp_callback
f_f77, jac_f77_radau5 = tmp_callback.f_f77, tmp_callback.jac_f77_radau5

t0, tn, u0, n_points = 0., 4.,  [1.,0.,0.], 10
time_points = np.linspace(t0, tn, n_points)
atol, rtol = [1e-6,1e-8,1e-6], 1e-4

st.figure()
method = Radau5Implicit
exact_final = [9.055142e-1, 2.240418e-5, 9.446344e-2]

# Test case 1: Radau5, with f, mas & jac
m = method(None, f_f77=f_str, rtol=rtol, atol=atol,
           jac_f77_radau5=jac_str, mas_f77=mas_str)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'r-', title="Radau5 with Fortran subroutines",
        legend="with f, mas & jac", hold="on")
print 'Max error for test case 1 is %g' % max(u[-1] - exact_final)

# Test case 2: Radau5, with f & mas
m = method(None, f_f77=f_str, rtol=rtol, atol=atol, mas_f77=mas_str)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'g*', title="Radau5 with Fortran subroutines",
        legend="with f & mas", hold="on")
print 'Max error for test case 2 is %g' % max(u[-1] - exact_final)

os.remove('tmp_callback.so')

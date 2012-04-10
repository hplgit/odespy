# Author: Liwei Wang

"""
Same example as in demo_Radau5_2.py:
Van Der Pol oscillator
u''=3*(1-u*u)*u'-u, with supplied full jacobian matrix

This example intends to show users how to improve efficiency with
call-back functions composed in Fortran language.
"""

from odespy import *
import scitools.std as st
import numpy as np

f_str = """
      subroutine f_f77(neq, t, u, udot)
Cf2py intent(hide) neq
Cf2py intent(out) udot
      integer neq
      double precision t, u, udot
      dimension u(neq), udot(neq)
      udot(1) = u(2)
      udot(2) = 3d0*(1d0 - u(1)**2)*u(2) - u(1)
      return
      end
"""
f_str = """
      subroutine f_f77(neq, t, u, udot)
Cf2py intent(hide) neq
Cf2py intent(out) udot
      integer neq
      double precision t, u, udot
      dimension u(neq), udot(neq)
      udot(1) = u(2)
      udot(2) = 3*(1 - u(1)**2)*u(2) - u(1)
      return
      end
"""

jac_full_str = """
      subroutine jac_f77_radau5(neq,t,u,dfu,ldfu,rpar,ipar)
Cf2py intent(hide) neq,rpar,ipar
Cf2py intent(in)   t,u,ldfu
Cf2py intent(out) dfu
      integer neq,ipar,ldfu
      double precision t,u,dfu,rpar
      dimension u(neq),dfu(ldfu,neq),rpar(*),ipar(*)
      dfu(1,1) = 0d0
      dfu(1,2) = 1d0
      dfu(2,1) = -6d0*u(1)*u(2) - 1d0
      dfu(2,2) = 3d0*(1d0-u(1)**2)
      return
      end
"""

jac_full_str = """
      subroutine jac_f77_radau5(neq,t,u,dfu,ldfu,rpar,ipar)
Cf2py intent(hide) neq,rpar,ipar
Cf2py intent(in)   t,u,ldfu
Cf2py intent(out) dfu
      integer neq,ipar,ldfu
      double precision t,u,dfu,rpar
      dimension u(neq),dfu(ldfu,neq),rpar(*),ipar(*)
      dfu(1,1) = 0
      dfu(1,2) = 1
      dfu(2,1) = -6*u(1)*u(2) - 1
      dfu(2,2) = 3*(1d0-u(1)**2)
      return
      end
"""

import sys
try:
    n_points = int(sys.argv[1])    # Read from input
except:
    n_points = 10   # default number of time-steps

# Compile these Fortran subroutines
f_f77, jac_f77 = compile_f77(f_str, jac_full_str)

t0, tn, u0 = 0., 10.,  [2.,0.]
time_points = np.linspace(t0, tn, n_points)
atol, rtol = 1e-4, 1e-4


st.figure()
method = Radau5Explicit

# Test case 1: Radau5, with f & jac
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol, jac_f77_radau5=jac_f77)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], title= "Van der Pol oscillator, with Radau5 & Lsoda",
        legend="Radau5 with f & jac", hold="on")

# Test case 2: Radau5, with f
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', legend="Radau5 with f", hold="on")

method = Lsoda
# Test case 3: Lsoda, with f
m = m.switch_to(Lsoda)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], legend="Lsoda with f & jac", hold="on")

# Test case 4: Lsoda, with f & jac
m.set(jac_f77_radau5=jac_f77)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', legend="Lsoda with f", hold="on")

os.remove('tmp_callback.so')

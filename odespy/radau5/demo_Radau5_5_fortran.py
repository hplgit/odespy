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
      integer i
      udot(1) = -1.25*(u(1)**2)+1.25*(u(2)-2.*u(1))
      do 10 i = 2,8
   10   udot(i) = 1.25*(u(i-1)**2-u(i+1)**2)+1.25*
     1     (u(i+1)-2.*u(i)+u(i-1))
      udot(9) = 1.25*(u(8)**2) + 1.25*(u(8)-2.*u(9))
      return
      end

"""

jac_str = """
      subroutine jac_f77_radau5(neq,t,u,dfu,ldfu,rpar,ipar)
Cf2py intent(hide) neq
Cf2py intent(hide) rpar,ipar
Cf2py intent(in)   t,u,ldfu
Cf2py intent(out)  dfu
      integer neq,ipar,ldfu
      double precision t,u,dfu,rpar
      dimension u(neq),dfu(ldfu,neq),rpar(*),ipar(*)
      integer i
      dfu(1,1) = -2.5
      dfu(1,2) = 1.25-u(2)*2.5
      do 10 i=2,8
        dfu(i,i-1) = 1.25+u(i-1)*2.5
        dfu(i,i  ) = -2.5
   10   dfu(i,i+1) = 1.25 - u(i+1)*2.5
      dfu(9,8) = 1.25+u(8)*2.5
      dfu(9,9) = -2.5
      return
      end
"""

jac_banded_str = """
      subroutine jac_f77_radau5(neq,t,u,dfu,ldfu,rpar,ipar)
Cf2py intent(hide) neq,rpar,ipar
Cf2py intent(in) t,u,ldfu
Cf2py intent(out) dfu
      integer neq,ipar,ldfu
      double precision t,u,dfu,rpar
      dimension u(neq),dfu(ldfu,neq)
      integer i
      dfu(2,1) = -2.5
      dfu(1,2) = 1.25-u(2)*2.5
      do 10 i=2,8
        dfu(3,i-1) = 1.25+u(i-1)*2.5
        dfu(2,i  ) = -2.5
   10   dfu(1,i+1) = 1.25-u(i+1)*2.5
      dfu(3,8) = 1.25+u(9)*2.5
      dfu(2,9) = -2.5
      return
      end
"""

mas_str = """
      subroutine mas_f77(neq,mas,lmas,rpar,ipar)
Cf2py intent(hide) ipar,rpar
Cf2py intent(in)   neq,lmas
Cf2py intent(out)  mas
      integer neq,lmas,ipar
      double precision mas(lmas,neq),rpar
      integer i
      mas(1,1) = 4./6.
      mas(1,2) = 1./6.
      do 10 i = 2,8
        mas(i,i-1) = 1./6.
        mas(i,i)   = 4./6.
   10   mas(i,i+1) = 1./6.
      mas(9,8) = 1./6.
      mas(9,9) = 4./6.
      return
      end
"""

mas_banded_str = """
      subroutine mas_f77(neq,mas,lmas,rpar,ipar)
Cf2py intent(hide)   rpar,ipar
Cf2py intent(in)     neq,lmas
Cf2py intent(out)    mas
      integer neq,lmas,ipar(*)
      double precision mas(lmas,neq),rpar(*)
      integer i
      do 10 i = 1,9
        mas(1,i) = 1./6.
        mas(2,i) = 4./6.
   10   mas(3,i) = 1./6.
      mas(1,1) = 0.
      mas(3,9) = 0.
      return
      end
"""

# compile these Fortran subroutines

string_to_compile = '\n'.join([f_str, jac_str, mas_str])
from numpy import f2py
f2py.compile(string_to_compile, modulename = "tmp_callback", verbose=False)
import tmp_callback
f_f77, mas_f77, jac_f77_radau5 = \
    tmp_callback.f_f77, tmp_callback.mas_f77, tmp_callback.jac_f77_radau5

u0 = np.zeros(9, float)
u0[1], u0[2:5], u0[5] = .5, 1., .5

t0, tn, n_points, ml, mu, mlmas, mumas = 0., .4, 50, 1, 1, 1, 1
time_points = np.linspace(t0, tn, n_points)
atol, rtol = 1e-3, 1e-3

exact_final = [1.07001457e-1, 2.77432492e-1, 5.02444616e-1, 7.21037157e-1,
               9.01670441e-1, 8.88832048e-1, 4.96572850e-1, 9.46924362e-2,
               -6.90855199e-3]
st.figure()
method = Radau5Implicit

# Test case 1: Radau5, with f, mas & jac
m = method(None, f_f77=f_str, mas_f77=mas_f77, rtol=rtol, atol=atol,
           jac_f77_radau5=jac_str)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'b-', title="Radau5 with Fortran subroutines",
        legend="with f, mas & jac_full", hold="on")
print 'Max error for test case 1 is %g' % max(u[-1] - exact_final)

# Test case 2: Radau5, with f, mas
m = method(None, f_f77=f_f77, mas_f77=mas_f77, rtol=rtol, atol=atol)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'r*', title="Radau5 with Fortran subroutines",
        legend="with f, mas", hold="on")
print 'Max error for test case 2 is %g' % max(u[-1] - exact_final)

# compile these Fortran subroutines
os.remove('tmp_callback.so')
string_to_compile = '\n'.join([f_str, jac_banded_str, mas_banded_str])
from numpy import f2py
f2py.compile(string_to_compile, modulename = "tmp_callback2", verbose=False)
import tmp_callback2
f_f77, mas_banded, jac_banded = \
    tmp_callback2.f_f77, tmp_callback2.mas_f77, tmp_callback2.jac_f77_radau5

# Test case 3: Radau5, with f, mas_banded, ml, mu & jac_banded
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol,
           mlmas=mlmas, mumas=mumas, mas_f77=mas_banded,
           ml=ml, mu=mu,
           jac_f77_radau5=jac_banded)

m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'b-', title="Radau5 with Fortran subroutines",
        legend="with f, mas_banded & jac_banded", hold="on")
print 'Max error for test case 3 is %g' % max(u[-1] - exact_final)

os.remove('tmp_callback2.so')


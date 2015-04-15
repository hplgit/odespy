# Author: Liwei Wang

"""
u'= A*y, With supplied banded jacobian matrix

This example is the typical usage of Lsode/Lsoda with
user-supplied functions composed in Python.
"""

from odespy import *
import scitools.std as st
import numpy as np

f_str = """
      subroutine f_f77(neq, t, y, ydot)
Cf2py intent(hide) neq
Cf2py intent(out) ydot
      integer neq
      double precision t, y, ydot
      dimension y(neq), ydot(neq)
      integer i, j, k, ng
      double precision alph1, alph2, d
      data alph1/1.0d0/, alph2/1.0d0/, ng/5/
      do 10 j = 1,ng
      do 10 i = 1,ng
        k = i + (j - 1)*ng
        d = -2.0d0*y(k)
        if (i .ne. 1) d = d + y(k-1)*alph1
        if (j .ne. 1) d = d + y(k-ng)*alph2
 10     ydot(k) = d
      return
      end
"""

jac_banded_str = """
      subroutine jac_f77_radau5(neq,t,u,dfu,ldfu,rpar,ipar)
Cf2py intent(hide) neq,rpar,ipar
Cf2py intent(in)   t,u,ldfu
Cf2py intent(out)  dfu
      integer neq,ipar,ldfu
      double precision t,u,dfu,rpar
      dimension u(neq),dfu(ldfu,neq),rpar(*),ipar(*)
      integer i,j
      do 10 j = 1,neq
        dfu(1,j) = -2.0d0
        dfu(2,j) = 1.0d0
 10     dfu(6,j) = 1.0d0
      do 20 j = 5,neq,5
 20     dfu(2,j) = 0.0d0
      return
      end
"""

# Compile these Fortran subroutines
from numpy import f2py
f2py.compile(f_str+'\n'+jac_banded_str, modulename='tmp_callback', verbose=False)
import tmp_callback
f_f77, jac_banded_f77 = tmp_callback.f_f77, tmp_callback.jac_f77_radau5

import sys
try:
    n_points = int(sys.argv[1])    # Read from input
except:
    n_points = 10   # default number of time-steps


t0, tn, u0 = 0., 4.,  [1]+24*[0]
ml, mu = 5, 0
time_points = np.linspace(t0, tn, n_points)
atol, rtol = 1e-2, 1e-2


st.figure()
method = Radau5Explicit

# Test case 1: Radau5, with f, ml, mu & jac_banded
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol,
           ml=ml, mu=mu, jac_f77_radau5=jac_banded_f77)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], title= "Radau5 & Lsoda",
        legend="Radau5 with f, ml, mu & jac", hold="on")

# Test case 2: Radau5, with f, ml, mu
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol, ml=ml, mu=mu)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*',
        legend="Radau5 with f & ml,mu", hold="on")

# Test case 3: Radau5, with f
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'o',
        legend="with f", hold="on")

method = Lsoda
st.figure()
# Test case 4: Lsoda, with f, ml, mu & jac_banded
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol,
           ml=ml, mu=mu)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0],
        legend="Lsoda with f, ml, mu & jac", hold="on")

# Test case 5: Lsoda, with f
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'o',
        legend="Lsoda with f", hold="on")

os.remove('tmp_callback.so')

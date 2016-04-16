# Author: Liwei Wang

"""
Same example as in demo_Lsode_Lsoda_2.py:
u'= A*y, With supplied banded jacobian matrix

This example intends to show users how to improve efficiency with 
call-back functions composed in Fortran language. 

If there are more complicated mathematical calculation and long loops 
involved in user-supplied functions, efficiency would be improved a 
lot with composing these functions in Fortran Language. 

"""

from odespy import *
import scitools.std as st
import numpy as np
import os

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

jac_str = """
      subroutine jac_banded_f77(neq, t, y, ml, mu, pd, nrowpd)
Cf2py intent(hide) neq, ml, mu, nrowpd
Cf2py intent(out) pd
      integer neq, ml, mu, nrowpd, j, mband, mu1, mu2, ng
      double precision t, y, pd, alph1, alph2
      dimension y(neq), pd(nrowpd,neq)
      data alph1/1.0d0/, alph2/1.0d0/, ng/5/
      mband = ml + mu + 1
      mu1 = mu + 1
      mu2 = mu + 2
      do 10 j = 1,neq
        pd(mu1,j) = -2.0d0
        pd(mu2,j) = alph1
 10     pd(mband,j) = alph2
      do 20 j = ng,neq,ng
 20     pd(mu2,j) = 0.0d0
      return
      end
"""

import sys
try:
    n_points = int(sys.argv[1])    # Read from input 
except:
    n_points = 10   # default number of time-steps    


t0, tn, u0 = 0., 4.,  [1]+24*[0]
ml, mu = 5, 0
time_points = np.linspace(t0, tn, n_points)
atol, rtol = 1e-2, 1e-2

# Compile these Fortran subroutines
from numpy import f2py
f2py.compile(f_str+'\n'+jac_str, modulename='callback', verbose=False)   
import callback
f_f77, jac_banded_f77 = callback.f_f77, callback.jac_banded_f77


st.figure()
method = Lsode

# Test case 1: Lsode, with f, ml, mu & jac_banded
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol, 
           ml=ml, mu=mu, jac_banded_f77=jac_banded_f77)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], title="Lsode with Fortran subroutines",
        legend="with f, ml, mu & jac_banded", hold="on")

# Test case 2: Lsode, with f, ml, mu
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol, 
           ml=ml, mu=mu)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Lsode with Fortran subroutines", 
        legend="with f & jac_banded", hold="on")

# Test case 3: Lsode, with f
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'o', title="Lsode with Fortran subroutines", 
        legend="with f", hold="on")

method = Lsoda
st.figure()
# Test case 4: Lsoda, with f, ml, mu & jac_banded
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol, 
           ml=ml, mu=mu, jac_banded_f77=jac_banded_f77)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], title="Lsoda with Fortran subroutines",
        legend="with f, ml, mu & jac_banded", hold="on")

# Test case 5: Lsoda, with f, ml, mu
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol, 
           ml=ml, mu=mu)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Lsoda with Fortran subroutines", 
        legend="with f & jac_banded", hold="on")

# Test case 6: Lsoda, with f
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'o', title="Lsoda with Fortran subroutines", 
        legend="with f", hold="on")

st.figure()
# Test case 7: Switch to Lsodar
m.switch_to(Lsodar)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '-', title="Lsodar with Fortran subroutines", 
        legend="Switch from Lsoda", hold="on")

# Test case 8: Supplement ml, mu
m.set(ml=ml, mu=mu)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Lsodar with Fortran subroutines", 
        legend="Supplement ml, mu", hold="on")

# Test case 9: Supplement jac_banded as string of Fortran code
m.set(jac_banded_f77=jac_str, iter_method=None)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'o', title="Lsodar with Fortran subroutines", 
        legend="Supplement jac_banded", hold="on")

os.remove('callback.so')

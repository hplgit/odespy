# Author: Liwei Wang

"""

Same example as in demo_Lsodes.py
u'= A*u, where A is assumed to be sparse matrix

This example is the typical usage of Lsodes with 
user-supplied functions composed in Fortran code.
"""

from odespy import *
#import scitools.basics,easyviz as st
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
      integer i, igrid, j, l, m
      data igrid/3/

      do 5 i = 1,neq
 5      udot(i) = 0.0d0
      do 20 m = 1,igrid
        do 10 l = 1,igrid
          j = l + (m - 1)*igrid
          if (m .ne. 1) udot(j-igrid) = udot(j-igrid) + u(j)
          if (l .ne. 1) udot(j-1) = udot(j-1) + u(j)
          udot(j) = udot(j) - 4.0d0*u(j)
          if (l .ne. igrid) udot(j+1) = udot(j+1) + u(j)
 10       continue
 20     continue
      return
      end
"""

jac_str = """
      subroutine jac_column_f77(neq, t, u, j, ia, ja, pd)
Cf2py intent(hide) neq, ia, ja
Cf2py intent(out) pd
      integer neq, j, ia, ja
      double precision t, u, pd
      dimension u(neq), ia(neq + 1), ja(*), pd(neq)
      integer igrid, l, m
      data igrid/3/
      m = (j - 1)/igrid + 1
      l = j - (m - 1)*igrid
      pd(j) = -4.0d0
      if (m .ne. 1) pd(j-igrid) = 1.0d0
      if (l .ne. 1) pd(j-1) = 1.0d0
      if (l .ne. igrid) pd(j+1) = 1.0d0
      return
      end
"""

import sys
try:
    n_points = int(sys.argv[1])    # Read from input 
except:
    n_points = 20   # default number of time-steps    

# Compile these Fortran subroutines
from numpy import f2py
f2py.compile(f_str+'\n'+jac_str, modulename='callback', verbose=False)   
import callback
f_f77, jac_column_f77 = callback.f_f77, callback.jac_column_f77

t0, tn, u0 = 0., 3., np.arange(1.,10.,1.)
time_points = np.linspace(t0, tn, n_points)

atol, rtol = 1e-4, 1e-4
ia = [1,3,6,8,11,15,18,21,25,28]   # Describe the sparse structure
ja = [1,2,1,2,3,2,3,1,4,5,2,4,5,6,3,5,6,4,7,8,5,7,8,9,6,8,9]

st.figure()
method = Lsodes

# Test case 1: Lsodes, with f, ia, ja & jac_column
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol, 
           ia=ia, ja=ja, jac_column_f77=jac_column_f77)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], title="Lsodes with Fortran subroutines",
        legend="with f, ia, ja & jac_column", hold="on")

# Test case 2: Lsodes, with f, ia, ja 
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol, 
           ia=ia, ja=ja)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Lsodes with Fortran subroutines", 
        legend="with f, ia & ja", hold="on")

# Test case 3: Lsodes, with f, jac_column
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol, 
           jac_column_f77=jac_column_f77)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Lsode with Fortran subroutines", 
        legend="with f & jac_column", hold="on")


os.remove('callback.so')

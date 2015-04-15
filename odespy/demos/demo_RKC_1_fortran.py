# Author: Liwei Wang

"""
This example intends to show users of RKC how to improve efficiency with
call-back functions composed in Fortran language.

Same example as in demo_RKC_1.py with call-back function composed
 directly in Fortran.
"""
from odespy import *
import scitools.std as st
import numpy as np
import os

f_str = """
      subroutine f_f77(neqn,t,y,dy)
c---------------------------------------------------------------
c  Semi-discretization of reaction-diffusion equation by central
c  differences.  The analytical solution sol(x,t) is used for
c  Dirichlet boundary conditions at x = 0 and x = 10.
c---------------------------------------------------------------
Cf2py intent(hide)   neqn
Cf2py intent(out)    dy
      integer           neqn
      double precision  t,y(neqn),dy(neqn)
      integer           i
      double precision  sol, sol0, sol10
      dy(1) = (sol(0d0,t) - 2d0*y(1) + y(2))*1d2 +
     &        (1d0 - y(1))*y(1)**2
      do 10 i = 2,neqn-1
         dy(i) = (y(i-1) - 2d0*y(i) + y(i+1))*1d2 +
     &           (1d0 - y(i))*y(i)**2
10    continue
      dy(neqn) = (y(neqn-1)- 2d0*y(neqn) + sol(10d0,t))*1d2 +
     &           (1d0 - y(neqn))*y(neqn)**2
      return
      end

      double precision function sol(x,t)
c------------------------------------------------------------
c  An analytical solution to the reaction-diffusion equation.
c------------------------------------------------------------
      double precision x,t
      double precision v,z
      v = sqrt(0.5d0)
      z = x - v*t
      sol = 1d0/(1d0 + exp(v*z))
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
f2py.compile(f_str, modulename='callback', verbose=False)
import callback
f_f77 = callback.f_f77

t0, tn = 0., 15.
u0 = [1.0/(1.0 + np.exp(np.sqrt(.5)*i*0.1)) for i in range(1,100)]
time_points = np.linspace(t0, tn, n_points)

atol, rtol = 1e-4, 1e-4

st.figure()
# Test case 1: RKC
m = RKC(None, f_f77=f_f77, rtol=rtol, atol=atol)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], title="Fortran subroutine f",
        legend="Rkc", hold="on")

# Test case 2: Lsoda
m = Lsode(None, f_f77=f_f77, rtol=rtol, atol=atol)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Fortran subroutine f",
        legend="Lsoda", hold="on")

# Test case 3: RKFehlberg
m = m.switch_to(RKFehlberg)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Fortran subroutine f",
        legend="RKFehlberg", hold="on")

os.remove('callback.so')

# Author: Liwei Wang

"""
This example intends to show users of RKC how to improve efficiency with 
call-back functions composed in Fortran language. 

Same example as in test_RKC_2.py with call-back function composed
 directly in Fortran. 
"""
from odespy import *
import scitools.std as st
import numpy as np
import os

f_str = """
      subroutine f_f77(neq,t,u,udot)
Cf2py intent(hide)   neq
Cf2py intent(out)    udot
      integer          neq
      double precision t,u(neq),udot(neq)
      integer          i,j,k,l
      double precision v(0:20,0:20,0:20),dx,dy,dz,dxsq,dysq,dzsq,
     &                 arg,sh,ch,sol
      integer          nx,ny,nz
      data             nx/19/, ny/19/, nz/19/
c
      dx = 1d0/(nx+1)
      dy = 1d0/(ny+1)
      dz = 1d0/(nz+1)
      dxsq = dx*dx
      dysq = dy*dy
      dzsq = dz*dz
      do 30 i = 1,nx
        do 20 j = 1,ny
          do 10 k = 1,nz
            v(i,j,k) = u(i + (j-1)*nx + (k-1)*nx*ny)
10        continue
20      continue   
30    continue   
c
      do 50 i = 1,nx
        do 40 j = 1,ny
          v(i,j,0) = sol(i*dx,j*dy,0d0,t)
          v(i,j,nz+1) = sol(i*dx,j*dy,1d0,t)
40      continue
50    continue
c
      do 70 i = 1,nx
        do 60 k = 1,nz
          v(i,0,k) = sol(i*dx,0d0,k*dz,t)
          v(i,ny+1,k) = sol(i*dx,1d0,k*dz,t)
60      continue
70    continue
c
      do 90 j = 1,ny
        do 80 k = 1,nz
          v(0,j,k) = sol(0d0,j*dy,k*dz,t)
          v(nx+1,j,k) = sol(1d0,j*dy,k*dz,t)
80      continue          
90    continue
c
      do 120 i = 1,nx
        do 110 j = 1,ny
          do 100 k = 1,nz
            arg = 5d0*(i*dx + 2d0*j*dy + 1.5d0*k*dz - 0.5d0 - t)
            sh = sinh(arg)
            ch = cosh(arg)
            l = i + (j-1)*nx + (k-1)*nx*ny
            udot(l) = (v(i-1,j,k) - 2d0*v(i,j,k) + v(i+1,j,k))/dxsq +
     &                (v(i,j-1,k) - 2d0*v(i,j,k) + v(i,j+1,k))/dysq +
     &                (v(i,j,k-1) - 2d0*v(i,j,k) + v(i,j,k+1))/dzsq +
     &                (-5d0*ch + 362.5d0*sh)/(ch**3)
100       continue
110     continue
120   continue
c
      return
      end


      double precision function sol(x,y,z,t)
      double precision x,y,z,t
      double precision arg
      arg = 5d0*(x + 2d0*y + 1.5d0*z - 0.5d0 - t)
      sol = tanh(arg)
      return
      end
"""


spcrad_str = """
      double precision function spcrad_f77(neqn,t,y)
Cf2py intent(hide)     neqn
      integer          neqn
      double precision t,y(neqn)
      integer          nx,ny,nz
      data             nx/19/, ny/19/, nz/19/ 
      spcrad_f77 = 4d0*((nx+1)**2 + (ny+1)**2 + (nz+1)**2)
      return
      end

"""

import sys
try:
    n_points = int(sys.argv[1])    # Read from input 
except:
    n_points = 20   # default number of time-steps    

# Compile these Fortran subroutines
string_to_compile = '\n'.join([f_str, spcrad_str])
from numpy import f2py
f2py.compile(string_to_compile, modulename='callback', verbose=False)   
import callback
f_f77, spcrad_f77 = callback.f_f77, callback.spcrad_f77

t0, tn = 0., .7 
u0 = np.zeros(6859,float)
for i in range(19):
    for j in range(19):
        for k in range(19):
            index = i + j*19 + k*19*19
            u0[index] = np.tanh(5.0*((i+1)*.05 + (j+1)*.1 + \
                                (k+1)*.075 - .5))
time_points = np.linspace(t0, tn, n_points)
atol, rtol = 1e-4, 1e-4
jac_constant = 1

st.figure()
method = RKC
# Test case 1: Rkc with f, jac_constant, spcrad
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol, 
           spcrad_f77=spcrad_f77, jac_constant=jac_constant)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], title="Rkc with Fortran subroutines",
        legend="with f, jac_constant, spcrad", hold="on")

# Test case 2: Rkc with f, jac_constant
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol, 
           jac_constant=jac_constant)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], '*', title="Rkc with Fortran subroutines",
        legend="with f, jac_constant", hold="on")

# Test case 3: Rkc with f, spcrad
m = method(None, f_f77=f_f77, rtol=rtol, atol=atol, 
           spcrad_f77=spcrad_f77)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], title="Rkc with Fortran subroutines",
        legend="with f, spcrad", hold="on")

os.remove('callback.so')



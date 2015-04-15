# Author: Liwei Wang

"""
Test for dlsoibt(). Demo from opkdmain.f  in ODEPACK

The Burgers equation.
du/dt  =  - u * du/dx  +  eta * d**2 u/dx**2,   eta = .05,
on -1 .le. x .le. 1.  The boundary conditions are
du/dx = 0  at x = -1 and at x = 1.
The initial profile is a square wave,
u = 1 in ABS(x) .lt. .5,  u = .5 at ABS(x) = .5,  u = 0 elsewhere.
The result is a system A * dy/dt = g(y),
of size NEQ = 41, where y(i) is the approximation to u at x = x(i),
with x(i) = -1 + (i-1)*delx, delx = 2/(NEQ-1) = .05.  The individual
equations in the system are
dy(1)/dt = ( y(3) - 2*y(2) + y(1) ) * eta / delx**2,
dy(NEQ)/dt = ( y(NEQ-2) - 2*y(NEQ-1) + y(NEQ) ) * eta / delx**2,
and for i = 2, 3, ..., NEQ-1,
(1/6) dy(i-1)/dt + (4/6) dy(i)/dt + (1/6) dy(i+1)/dt
  = ( y(i-1)**2 - y(i+1)**2 ) / (4*delx)
    + ( y(i+1) - 2*y(i) + y(i-1) ) * eta / delx**2.
"""
from odespy import *
#import scitools.basics,easyviz as st
import scitools.std as st
import numpy as np

def res(y, t, s, ires):
    r = np.zeros(41, float)
    r[0] = 20.*(y[2] - 2.0*y[1] + y[0]) - s[0]
    r[1:40] = (y[0:39]**2 - y[2:41]**2)*5. + \
        20.*(y[2:41] - 2*y[1:40] + y[0:39]) - \
        (s[0:39] + 4.*s[1:40] + s[2:41])/6.
    r[40] = 20.*(y[38] - 2.*y[39] + y[40]) - s[40]
    return r, ires

def adda(y, t, pa, pb, pc):
    pa[0][0][:] += 4./6.
    pb[0][0][:] += 1./6.
    pc[0][0][:] += 1./6.
    pa[0][0][0] += 2./6.
    pa[0][0][40] += 2./6.
    return pa, pb, pc

def jac(y, t, s):
    pa, pb, pc = np.zeros(41,float), np.zeros(41,float), np.zeros(41,float)
    pa.shape = pb.shape = pc.shape = (1, 1, 41)
    pa[0][0][0], pb[0][0][0], pc[0][0][0] = 20., -40., 20.
    pa[0][0][1:40] = -40. 
    pb[0][0][1:40] = 20. - y[2:41]*10.
    pc[0][0][1:40] = 20. + y[0:39]*10.
    pa[0][0][40], pb[0][0][40], pc[0][0][40] = 20. , 20., -40.
    return pa,pb,pc

import sys
try:
    n_points = int(sys.argv[1])    # Read from input 
except:
    n_points = 5   # default number of time-steps    

t0, tn, u0 = 0., .4, [0.]*10+[.5]+[1.]*19+[.5]+[0.]*10
time_points = np.linspace(t0, tn, n_points)
atol, rtol = 1e-3, 1e-3
mb, nb = 1, 41

st.figure()
method = Lsoibt

# Test case 1: Lsoibt, with res, adda, jac, mb, nb
m = method(rtol=rtol, atol=atol, res=res, adda_lsoibt=adda, 
           mb=mb, nb=nb, jac_lsoibt=jac)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'b-', title="Lsoibt with Python functions",
        legend="with res, jac & adda", hold="on")

# Test case 2: Lsoibt, with res, adda, mb, nb
m = method(rtol=rtol, atol=atol, res=res, adda_lsoibt=adda, 
           mb=mb, nb=nb)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'g*', title="Lsoibt with Python functions",
        legend="with res & adda", hold="on")



# Author: Liwei Wang

"""
Test example from opkddemos.f in ODEPACK.

Simplified Galerkin Solution of Burgers Equation

     Diffusion coefficient is eta =  0.50E-01
     Uniform mesh on interval  -0.100E+01 to    0.100E+01
     Zero boundary conditions
     Time limits: t0 =  0.00000E+00   tlast =  0.40000E+00
     Half-bandwidths ml =  1   mu =  1
     System size neq =   9

    this program solves a semi-discretized form of the Burgers equation,
    u  = -(u*u/2)  + eta * u
     t           x          xx
    for a = -1 .le. x .le. 1 = b, t .ge. 0.
    Here eta = 0.05.
    Boundary conditions: u(-1,t) = u(1,t) = 0.
    Initial profile: square wave
    u(0,x) = 0    for 1/2 .lt. abs(x) .le. 1
    u(0,x) = 1/2  for abs(x) = 1/2
    u(0,x) = 1    for 0 .le. abs(x) .lt. 1/2
    Initial profile:
        0.0000E+00  0.0000E+00  0.5000E+00  0.1000E+01  0.1000E+01  0.1000E+01
        0.5000E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
"""

from odespy import *
import scitools.std as st
import numpy as np

def f(u, t):
    udot = np.zeros(9, float)
    udot[0] = -1.25*(u[0]**2) + 1.25*(u[1] - 2.*u[0])
    udot[1:8] = 1.25*(u[0:7]**2 - u[2:9]**2) + \
        1.25*(u[2:9] - 2.*u[1:8] + u[0:7])
    udot[8] = 1.25*(u[7]**2) + 1.25*(u[7] - 2.*u[8])
    return udot

def mas():
    m = [[4./6., 1./6.,    0., 0., 0., 0., 0., 0., 0.],
         [1./6., 4./6., 1./6., 0., 0., 0., 0., 0., 0.],
         [0., 1./6., 4./6., 1./6., 0., 0., 0., 0., 0.],
         [0., 0., 1./6., 4./6., 1./6., 0., 0., 0., 0.],
         [0., 0., 0., 1./6., 4./6., 1./6., 0., 0., 0.],
         [0., 0., 0., 0., 1./6., 4./6., 1./6., 0., 0.],
         [0., 0., 0., 0., 0., 1./6., 4./6., 1./6., 0.],
         [0., 0., 0., 0., 0., 0., 1./6., 4./6., 1./6.],
         [0., 0., 0., 0., 0., 0., 0.   , 1./6., 4./6.]]
    return m

def mas_banded():
    m = [[   0., 1./6., 1./6., 1./6., 1./6., 1./6., 1./6., 1./6., 1./6.],
         [4./6., 4./6., 4./6., 4./6., 4./6., 4./6., 4./6., 4./6., 4./6.],
         [1./6., 1./6., 1./6., 1./6., 1./6., 1./6., 1./6., 1./6.,    0.]]
    return m

def jac_full(u, t):   # Full jacobian
    p = np.zeros((9,9),float)
    p[0][0], p[0][1] =  -2.5, 1.25 - u[1]*2.5
    for i in range(1, 8):
        p[i][i-1] = 1.25 + u[i-1]*2.5
        p[i][i] = -2.5
        p[i][i+1] = 1.25 - u[i+1]*2.5
    p[8][7], p[8][8] = 1.25 + u[7]*2.5, -2.5
    return p

def jac_banded(u, t, ml, mu):   # Banded jacobian
    p = np.zeros((mu+2, 9),float)
    p[mu][0], p[mu-1][1] = -2.5, 1.25 - u[1]*2.5
    for i in range(1, 8):
        p[mu+1][i-1] = 1.25 + u[i-1]*2.5
        p[mu][i] = -2.5
        p[mu-1][i+1] = 1.25 - u[i+1]*2.5
    p[mu+1][7], p[mu][8] = 1.25 + u[8]*2.5, -2.5
    return p

u0 = np.zeros(9, float)
u0[1], u0[2:5], u0[5] = .5, 1., .5

t0, tn, n_points, ml, mu, mlmas, mumas = 0., .4, 5, 1, 1, 1, 1
time_points = np.linspace(t0, tn, n_points)
atol, rtol = 1e-3, 1e-3

exact_final = [1.07001457e-1, 2.77432492e-1, 5.02444616e-1, 7.21037157e-1,
               9.01670441e-1, 8.88832048e-1, 4.96572850e-1, 9.46924362e-2,
               -6.90855199e-3]
st.figure()
method = Radau5Implicit

# Test case 1: Radau5, with f, mas & jac
m = method(f=f, mas=mas, rtol=rtol, atol=atol, jac=jac_full)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'b-', title="Radau5 with Python functions",
        legend="with f, mas & jac_full", hold="on")
print 'Max error for test case 1 is %g' % max(u[-1] - exact_final)

# Test case 2: Radau5, with f, mas
m = method(f=f, mas=mas, rtol=rtol, atol=atol)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'r*', title="Radau5 with Python functions",
        legend="with f, mas", hold="on")
print 'Max error for test case 2 is %g' % max(u[-1] - exact_final)

# Test case 3: Radau5, with f, mas_banded, ml, mu & jac_banded
m = method(f=f, mas=mas_banded, rtol=rtol, atol=atol,
           ml=ml, mu=mu, jac_banded=jac_banded,
           mlmas=mlmas, mumas=mumas)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'b-', title="Radau5 with Python functions",
        legend="with f, mas & jac_banded", hold="on")
print 'Max error for test case 3 is %g' % max(u[-1] - exact_final)



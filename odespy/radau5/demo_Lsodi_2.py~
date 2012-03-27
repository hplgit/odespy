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
#import scitools.basics,easyviz as st
import scitools.std as st
import numpy as np

def res(u, t, s, ires):
    r = np.zeros(9,float)
    r[0] = -1.25*(u[1]**2) + 1.25*(u[1] - 2.*u[0])
    r[1:8] = 1.25*(u[0:7]**2 - u[2:9]**2) + \
        1.25*(u[2:9] - 2.*u[1:8] + u[0:7])
    r[8] = 1.25*(u[7]**2) + 1.25*(u[7] - 2.*u[8])
    if ires != -1:
        r[0] -= s[0]*2./3. + s[1]/6.
        r[1:8] -= s[0:7]/6. + s[1:8]*2./3. + s[2:9]/6.
        r[8] -= s[7]/6. + s[8]*2./3.
    return r,ires

def adda_full(u, t, p):  # Full form matrix A
    n = (p.shape)[1]
    fac1, fac4 = 1./6., 4./6.
    p[0][1] += fac1
    p[0][0] += fac4
    for i in range(1, n-1):
        p[i][i-1] += fac1
        p[i][i] += fac4
        p[i][i+1] += fac1
    p[n-1][n-1] += fac4
    p[n-1][n-2] += fac1
    return p

def adda_banded(u, t, p, ml, mu): # Banded matrix A
    fac1, fac4 = 1./6., 4./6.
    p[mu-1] += fac1
    p[mu] += fac4
    p[mu+1] += fac1
    return p

def jac_full(u, t, s):   # Full jacobian
    nps = len(u) + 1
    p = np.zeros((nps-1)**2,float)
    p.shape = (nps-1, nps-1)
    eodsq = .05/((2./nps)**2)
    diag = -2.*eodsq
    p[0][0], p[0][1] =  diag, eodsq - u[1]*nps/4.
    for i in range(1, nps-3):
        p[i][i-1] = eodsq + u[i-1]*nps/4.
        p[i][i] = diag
        p[i][i+1] = eodsq - u[i+1]*nps/4.
    p[nps-2][nps-3], p[nps-2][nps-2] = eodsq + u[nps-3]*nps/4., diag
    return p

def jac_banded(u, t, s, ml, mu):   # Banded jacobian
    nps = len(u) + 1
    eodsq = .05/((2./nps)**2)
    diag = -2.*eodsq
    r2d = nps/4.
    p = np.zeros((mu + 2, nps - 1),float)
    p[mu][0], p[mu-1][1] = diag, eodsq - u[1]*r2d
    for i in range(1, nps-2):
        p[mu+1][i-1] = eodsq + u[i-1]*r2d
        p[mu][i] = diag
        p[mu-1][i+1] = eodsq - u[i+1]*r2d
    p[mu+1][nps-3], p[mu][nps-2] = eodsq + u[nps-2]*r2d, diag
    return p

u0 = np.zeros(9, float)
u0[1], u0[2:5], u0[5] = .5, 1., .5

t0, tn, n_points, ml, mu = 0., .4, 5, 1, 1
time_points = np.linspace(t0, tn, n_points)
atol, rtol = 1e-3, 1e-3

exact_final = [1.07001457e-1, 2.77432492e-1, 5.02444616e-1, 7.21037157e-1,
               9.01670441e-1, 8.88832048e-1, 4.96572850e-1, 9.46924362e-2,
               -6.90855199e-3] 

st.figure()
method = Lsodi

# Test case 1: Lsodi, with res, adda_full, & jac_full
m = method(res=res, rtol=rtol, atol=atol, 
           adda_lsodi=adda_full, jac_lsodi=jac_full)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'r-', title="Lsodi with Python functions",
        legend="with res, adda_full & jac_full", hold="on")
print 'Max error with test case 1 is %g' % max(u[-1] - exact_final)


# Test case 2: Lsodi, with res & adda_full
m = method(res=res, rtol=rtol, atol=atol, adda_lsodi=adda_full)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'b*', title="Lsodi with Python functions",
        legend="with res & adda_full", hold="on")
print 'Max error with test case 2 is %g' % max(u[-1] - exact_final)

# Test case 3: Lsodi, with res, adda_banded, ml, mu, jac_banded
m = method(res=res, rtol=rtol, atol=atol, 
           adda_banded_lsodi=adda_banded,
           jac_banded_lsodi=jac_banded,
           ml=ml, mu=mu)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'go', title="Lsodi with Python functions",
        legend="with res, adda_banded, jac_banded, ml, mu", hold="on")
print 'Max error with test case 3 is %g' % max(u[-1] - exact_final)

# Test case 4: Lsodi, with res, adda_banded, ml, mu
m = method(res=res, rtol=rtol, atol=atol, 
           adda_banded_lsodi=adda_banded,
           ml=ml, mu=mu)
m.set_initial_condition(u0)
u,t = m.solve(time_points)
st.plot(t, u[:,0], 'y-', title="Lsodi with Python functions",
        legend="with res, adda_banded, ml, mu", hold="on")
print 'Max error with test case 4 is %g' % max(u[-1] - exact_final)




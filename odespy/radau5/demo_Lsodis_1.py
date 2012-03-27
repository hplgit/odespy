# Author: Liwei Wang

"""        
     Test case for dlsodis(). Demo from opkdmain.f  in ODEPACK
     The Burgers equation.
       du/dt  =  - u * du/dx  +  eta * d**2 u/dx**2,   eta = .05,
     on -1 .le. x .le. 1. The boundary conditions are periodic:
       u(-1,t) = u(1,t)  and  du/dx(-1,t) = du/dx(1,t)
     The initial profile is a square wave,
       u = 1 in ABS(x) .lt. .5,  u = .5 at ABS(x) = .5,  u = 0 elsewhere.

     The PDE is discretized in x by a simplified Galerkin method,
     using piecewise linear basis functions, on a grid of 40 intervals.
     The result is a system A * dy/dt = g(y), of size NEQ = 40,
     where y(i) is the approximation to u at x = x(i), with
     x(i) = -1 + (i-1)*delx, delx = 2/NEQ = .05.

     The individual equations in the system are (in order):
          (1/6)dy(NEQ)/dt+(4/6)dy(1)/dt+(1/6)dy(2)/dt
              = r4d*(y(NEQ)**2-y(2)**2)+eodsq*(y(2)-2*y(1)+y(NEQ))
      for i = 2,3,...,nm1,
     (1/6)dy(i-1)/dt+(4/6)dy(i)/dt+(1/6)dy(i+1)/dt
     = r4d*(y(i-1)**2-y(i+1)**2)+eodsq*(y(i+1)-2*y(i)+y(i-1))
     and finally
     (1/6)dy(nm1)/dt+(4/6)dy(NEQ)/dt+(1/6)dy(1)/dt
     = r4d*(y(nm1)**2-y(1)**2)+eodsq*(y(1)-2*y(NEQ)+y(nm1))
     where r4d = 1/(4*delx), eodsq = eta/delx**2 and nm1 = NEQ-1.
"""
from odespy import *
import numpy as np

def res(y, t, s, ires):
    r = np.zeros(40, float)
    r[0] = 5.*(y[39]**2 - y[1]**2) + 20.*(y[1] - 2.*y[0] + y[39])
    r[1:39] = 5.*(y[0:38]**2 - y[2:40]**2) + \
        20.*(y[2:40] - 2.*y[1:39] + y[0:38])
    r[39] = 5.*(y[38]**2 - y[0]**2) + 20.*(y[0] - 2.*y[39] + y[38])
    r[0] -= (s[39] + 4.*s[0] + s[1])/6.
    r[1:39] -= (s[0:38] + 4.*s[1:39] + s[2:40])/6.
    r[39] -= (s[38] + 4.*s[39] + s[0])/6.
    return r,ires

def adda(y, t, j, ia, ja, p):
    p[j] += 2./3.
    p[(j+1) % 40] += 1./6.
    p[(j-1) % 40] += 1./6.
    return np.asarray(p)

def jac(y, t, s, j, ia, ja):
    p = np.zeros(40, float)
    p[(j-1) % 40] = 10.*y[j] + 20.
    p[j] = -40.
    p[(j+1) % 40] = -10.*y[j] + 20.
    return np.asarray(p)

t0, tn, n_points = 0., .4, 5
u0 =  [0.]*10 + [.5] + [1.]*19 + [.5] + [0.]*9
time_points = np.linspace(t0, tn, n_points)
atol, rtol = 1e-3, 1e-3

method = Lsodis

exact_final = [1.8371e-2, 1.3578e-2, 1.5864e-2, 2.3805e-2, 3.7245e-2,
               5.6630e-2, 8.2538e-2, 1.1538e-1, 1.5522e-1, 2.0172e-1,
               2.5414e-1, 3.1150e-1, 3.7259e-1, 4.3608e-1, 5.0060e-1,
               5.6482e-1, 6.2751e-1, 6.8758e-1, 7.4415e-1, 7.9646e-1,
               8.4363e-1, 8.8462e-1, 9.1853e-1, 9.4500e-1, 9.6433e-1,
               9.7730e-1, 9.8464e-1, 9.8645e-1, 9.8138e-1, 9.6584e-1,
               9.3336e-1, 8.7497e-1, 7.8213e-1, 6.5315e-1, 4.9997e-1,
               3.4672e-1, 2.1758e-1, 1.2461e-1, 6.6208e-2, 3.3784e-2]

# Test case 1: with res, adda & jac
m = method(res=res, adda_lsodis=adda, atol=atol, rtol=rtol, jac_lsodis=jac)
m.set_initial_condition(u0)
y, t = m.solve(time_points)
print 'Max error for test case 1 is %g' % max(y[-1] - exact_final)

# Test case 2:  With res & adda
m = method(res=res, adda_lsodis=adda, atol=atol, rtol=rtol,lrw=4000,liw=100)
m.set_initial_condition(u0)
y, t = m.solve(time_points)
print 'Max error for test case 2 is %g' % max(y[-1] - exact_final)


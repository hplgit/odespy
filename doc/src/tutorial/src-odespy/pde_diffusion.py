"""Temperature evolution in a rod, computed by explicit odespy solvers."""

from numpy import linspace, zeros, linspace, array
import matplotlib.pyplot as plt
import time, sys

def rhs(u, t, L=None, beta=None, x=None):
    N = len(u) - 1
    dx = x[1] - x[0]
    rhs = zeros(N+1)
    rhs[0] = dsdt(t)
    for i in range(1, N):
        rhs[i] = (beta/dx**2)*(u[i+1] - 2*u[i] + u[i-1]) + \
                 f(x[i], t)
    rhs[N] = (beta/dx**2)*(2*u[i-1] - 2*u[i]) + f(x[N], t)
    return rhs

def K(u, t, L=None, beta=None, x=None):
    N = len(u) - 1
    dx = x[1] - x[0]
    K = zeros((N+1,N+1))
    K[0,0] = 0
    for i in range(1, N):
        K[i,i-1] = beta/dx**2
        K[i,i] = -2*beta/dx**2
        K[i,i+1] = beta/dx**2
    K[N,N-1] = (beta/dx**2)*2
    K[N,N] = (beta/dx**2)*(-2)
    return K

def rhs_vec(u, t, L=None, beta=None, x=None):
    N = len(u) - 1
    dx = x[1] - x[0]
    rhs = zeros(N+1)
    rhs[0] = dsdt(t)
    rhs[1:N] = (beta/dx**2)*(u[2:N+1] - 2*u[1:N] + u[0:N-1]) + \
               f(x[1:N], t)
    i = N
    rhs[i] = (beta/dx**2)*(2*u[i-1] - 2*u[i]) + f(x[N], t)
    return rhs

def K_vec(u, t, L=None, beta=None, x=None):
    """Vectorized computation of K."""
    N = len(u) - 1
    dx = x[1] - x[0]
    K = zeros((N+1,N+1))
    K[0,0] = 0
    K[1:N-1] = beta/dx**2
    K[1:N] = -2*beta/dx**2
    K[2:N+1] = beta/dx**2
    K[N,N-1] = (beta/dx**2)*2
    K[N,N] = (beta/dx**2)*(-2)
    return K

def s(t):
    return 423

def dsdt(t):
    return 0

def f(x, t):
    return 0


N = 40
L = 1
x = linspace(0, L, N+1)
f_kwargs = dict(L=L, beta=1, x=x)
u = zeros(N+1)

U_0 = zeros(N+1)
U_0[0] = s(0)
U_0[1:] = 283

import odespy
solvers = {
    'FE': odespy.ForwardEuler(
        rhs, f_kwargs=f_kwargs),
    'BE': odespy.BackwardEuler(
        rhs, f_is_linear=True, jac=K,
        f_kwargs=f_kwargs, jac_kwargs=f_kwargs),
    'B2': odespy.Backward2Step(
        rhs, f_is_linear=True, jac=K,
        f_kwargs=f_kwargs, jac_kwargs=f_kwargs),
    'theta': odespy.ThetaRule(
        rhs, f_is_linear=True, jac=K, theta=0.5,
        f_kwargs=f_kwargs, jac_kwargs=f_kwargs),
    'RKF': odespy.RKFehlberg(
        rhs, rtol=1E-6, atol=1E-8, f_kwargs=f_kwargs),
    'RKC': odespy.RKC(
        rhs, rtol=1E-6, atol=1E-8, f_kwargs=f_kwargs,
        jac_constant=True),
    }

try:
    method = sys.argv[1]
    dt = float(sys.argv[2])
    T = float(sys.argv[3])
except IndexError:
    method = 'FE'
    dx = x[1] - x[0]
    dt = dx**2/(2*beta) # Forward Euler limit
    print 'Forward Euler stability limit:', dt
    T = 1.2

solver = solvers[method]
solver.set_initial_condition(U_0)
N_t = int(round(T/float(dt)))
time_points = linspace(0, T, N_t+1)
u, t = solver.solve(time_points)

# Check how many time steps required by adaptive vs
# fixed-step methods
if hasattr(solver, 't_all'):
    print '# time steps:', len(solver.t_all)
    plt.figure()
    plt.plot(array(solver.t_all[1:]) - array(solver.t_all[:-1]))
    plt.title('Evolution of the time step in %s' %
              solver.__class__.__name__)
    plt.savefig('tmp.png'); plt.savefig('tmp.pdf')
    plt.show()
else:
    print '# time steps:', len(t)

# Make movie
import os
os.system('rm tmp_*.png')  # remove old plot files
plt.figure()
plt.ion()
y = u[0,:]
lines = plt.plot(x, y)
plt.axis([x[0], x[-1], 273, 1.2*s(0)])
plt.xlabel('x')
plt.ylabel('u(x,t)')
counter = 0
for i in range(0, u.shape[0]):
    print t[i]
    lines[0].set_ydata(u[i,:])
    plt.legend(['t=%.3f' % t[i]])
    plt.draw()
    if i % 5 == 0: # plot every 5 steps
        plt.savefig('tmp_%04d.png' % counter)
        counter += 1
    #time.sleep(0.2)

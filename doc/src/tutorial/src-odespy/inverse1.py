"""
Find the inverse of a function f(x).
Common solver: solve y = f(x) wrt x.
Here we differentiate to find x'(y) = 1/f'(x). x(0)=q where f(q)=0.
Solving the ODE for x(y) constructs the inverse of f.
"""

def rhs(x, y, f, h=1E-4):
    dfdx = (f(x+h) - f(x))/h
    return 1./dfdx

from numpy import sqrt, linspace

def f(x):
    return sqrt(x)
    #return 2*x

import odespy
solver = odespy.RungeKutta2(rhs, f_args=[f], f_kwargs={'h': 1E-3})
solver.set_initial_condition(0)
y_points = linspace(0, 4, 41)
x, y = solver.solve(y_points)
from matplotlib.pyplot import *
g_e = lambda y: y**2  # exact inverse function
y_ = linspace(0, 4, 11)
plot(y, x, 'r-',
     x, f(x), 'b-',
     y_, g_e(y_), 'go')
legend(['computed inverse', 'original f(x)', 'exact inverse'])
savefig('tmppng'); savefig('tmp.pdf')
show()




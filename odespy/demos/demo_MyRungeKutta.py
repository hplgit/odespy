# Author: Liwei Wang

"""
This example intends to show users how to apply MyRungeKutta to define 
own RungKutta solvers.
"""
from odespy import *
#import scitools.basics,easyviz as st
import scitools.std as st
import numpy as np

# arrays to test user-defined Runge-Kutta methods
bt = dict(
    FehlBerg = dict(
    table = np.array(\
        [[0., 0., 0., 0., 0., 0., 0.],
         [.25, .25, 0., 0., 0., 0., 0.],
         [.375, .09375, .28125, 0., 0., 0., 0.],
         [.92307692, .87938097, -3.27719618, 3.32089213, 0., 0., 0.],
         [1., 2.03240741,-8., 7.17348928,-.20589669, 0., 0.],
         [.5, -.2962963, 2., -1.38167641, .45297271, -.275, 0.],
         [0., .11574074, 0., .54892788, .53533138, -.2, 0.],
         [0., .11851852, 0., .51898635, .50613149, -.18, .03636364]]),
    order = (4,5)),
    RungeKutta3 =  dict(
    table = np.array(\
        [[0., 0., 0., 0.],
         [.5, .5, 0., 0.],
         [1., -1., 2., 0.],
         [0., .16666667, .66666667, .16666667]]),
    order = 3))

# Test for user-defined methods:
# Sample ODE problem: u = e**-t, T = [0.,.25,.5,.75,.1], 
def f(u,t):
    return -u
u0, t0, tn, n_points = 1., 0., 1., 5
time_points = np.linspace(t0, tn, n_points)

orders = []  # list for calculated orders
for m in bt.keys():
    # user-defined method, without order suplied 
    method = MyRungeKutta(f, butcher_tableau=bt[m]['table'])
    orders += [method.get_order()]
    method.set_initial_condition(u0)
    u,t = method.solve(time_points)
    error = abs((u[-1] - np.exp(-1.))/np.exp(-1.))
    print 'Error is %g with solver %s' % (error, m)





# Author: Liwei Wang

"""
This example intends to show users how to use AdaptivevResidual, i.e. 
to integrate with calculated residual as error-check criteria for 
a specified solver.

Currently, only scalar ODE problem can be used with AdaptiveResidual.

"""

from odespy import *
import numpy as np
import scitools.std as st

f = lambda u,t: -u
u0, t0, tn, n_points = 1., 0., 5., 10
atol, rtol = 1e-1, 1e-1
time_points = np.linspace(t0, tn, n_points)

for solver in ('RKFehlberg', 'RungeKutta2'):
    m = AdaptiveResidual(f, solver=solver, atol=atol, rtol=rtol)
    m.set_initial_condition(u0)
    u,t = m.solve(time_points, print_info=True)
    st.plot(t, u, '-', legend=solver, hold='on')

st.plot(time_points, np.exp(-time_points),'*', 
        legend='Exact', hold='on',
        title='AdaptiveResidual with rtol=atol=0.1')


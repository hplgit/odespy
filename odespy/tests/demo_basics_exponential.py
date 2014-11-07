# Author: Liwei Wang

# Scalar ODE:  Exponential
# u' = - u,   u = exp(-t)

from odespy import *
import numpy as np
import scitools.std as st

solver_no = 0   # number of solvers in figure
st.figure()

exceptions=['Lsodi', 'Lsodis', 'Lsoibt', 'MyRungeKutta',
            'MySolver', 'Lsodes', 'EulerCromer']
solvers = [solver for solver in list_available_solvers() if solver not in exceptions]

f=lambda u,t:-u
u0, time_points = 1., np.linspace(0., 10., 100)

print """Scalar ODE: Exponential u = exp(-t), u' = -u"""

# Loop for all possible solvers
for solver in solvers:
    try:
        method = eval(solver)(f, atol=1e-6)
        method.set_initial_condition(u0)
        u,t = method.solve(time_points)
        if solver_no % 8 == 0:
            # More than 8 solvers in current figure?
            st.figure()           # Initialize new figure.
        st.plot(t, u, hold="on", legend=solver, axis=[0., 10., 0., 1.5])
        solver_no += 1
        print 'Succeed when solver is %s' % solver
    except:
        print 'Failed when solver is %s' % solver

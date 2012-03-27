from Radau5 import *
from odespy import *
import scitools.std as st

st.figure()

exceptions=['RKC', 'Lsodes', 'Leapfrog', 'Lsodi', 'Lsodis', 'Lsoibt', 
                'MyRungeKutta', 'MySolver', 'AdaptiveResidual']

solvers = [solver for solver in list_all_solvers() if solver not in exceptions]

f = lambda (u00,u11),t: [u11, 3.*(1 - u00**2)*u11 - u00]
jac = lambda (u00,u11),t: \
    [[0., 1.], [-6.*u00*u11 - 1., 3.*(1. - u00**2)]]

u0, time_points = [2., 0.], np.linspace(0., 10., 10)

print """Van der Pol oscillator problem:
     u'' = 3*(1 - u**2)*u' - u"""

# Loop for all possible solvers
#for solver in solvers:

for solver in ['Vode']:
    method = eval(solver)(f)
    method.set_initial_condition(u0)
    u,t = method.solve(time_points)
    st.plot(t, u[:,0], '*', hold="on", legend=solver, axis=[0.,10.,-4.,4.])
    print 'Succeed when solver is %s' % solver


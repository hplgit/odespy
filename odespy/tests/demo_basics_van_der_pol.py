# Stiff ODE:  Van der Pol oscillator
# u'' = 3*(1 - u**2)*u' - u
from odespy import *
import scitools.std as st

solver_no = 1    # number of solvers in figure
st.figure()

exceptions=['RKC', 'Lsodes', 'Leapfrog', 'Lsodi', 'Lsodis', 'Lsoibt',
                'MyRungeKutta', 'MySolver', 'AdaptiveResidual', 'EulerCromer']
solvers = [solver for solver in list_available_solvers() if solver not in exceptions]

f = lambda (u00,u11),t: [u11, 3.*(1 - u00**2)*u11 - u00]
jac = lambda (u00,u11),t: \
    [[0., 1.], [-6.*u00*u11 - 1., 3.*(1. - u00**2)]]

u0, time_points = [2., 0.], np.linspace(0., 10., 100)

print """Van der Pol oscillator problem:
     u'' = 3*(1 - u**2)*u' - u"""

# Loop for all possible solvers
for solver in solvers:
    try:
        method = eval(solver)(f)
        method.set_initial_condition(u0)
        u,t = method.solve(time_points)
        if solver_no == 8:        # More than 8 solvers in current figure?
            st.figure()           # Initialize new figure.
            solver_no = 1
        st.plot(t, u[:,0], hold="on", legend=solver, axis=[0.,10.,-4.,4.])
        solver_no += 1
        print 'Succeed when solver is %s' % solver
    except:
        print 'Failed when solver is %s' % solver

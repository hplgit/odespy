# Author: Liwei Wang
# Scalar ODE with complex value 
# u' = 1/(t - 10 + 1j)

from odesolvers import *
import scitools.std as st

solver_no = 1    # number of solvers in figure
st.figure()

complex_solvers = ['AdamsBashMoulton2', 'AdamsBashMoulton3', 'AdamsBashforth2',
                   'AdamsBashforth3', 'AdamsBashforth4', 'AdaptiveResidual',
                   'Backward2Step', 'BackwardEuler', 'ForwardEuler', 
                   'Heun', 'Leapfrog', 'LeapfrogFiltered', 
                   'MidpointImplicit', 'MidpointIter', 'RungeKutta3', 
                   'RungeKutta2', 'RungeKutta4', 'RungeKuttaFehlberg', 
                   'ThetaRule', 'Trapezoidal']

fail_solvers = [name for name in list_solvers() if name not in complex_solvers]

f = lambda u, t: 1./(t - 10. + 1j)
u0, time_points = 0, np.linspace(0., 20., 200)

print """ Scalar ODE with complex value
 u' =  1./(t - 10 + 1j) """

succeed_no = 0
for solver in complex_solvers:
    try:      
        method = eval(solver)(f)
        method.set_initial_condition(u0)
        u,t = method.solve(time_points)
        success = np.allclose(u[-1].imag, -2*np.arctan(10))
        print solver, 'Succeed' if success else 'Fail'
    except:
        print 'Failed when solver is %s' % solver


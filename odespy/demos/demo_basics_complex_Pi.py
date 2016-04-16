# Author: Liwei Wang
# Scalar ODE with complex value
# u' = 1/(t - 10 + 1j)

from odespy import *

solver_no = 1    # number of solvers in figure

complex_solvers = ['AdamsBashMoulton2', 'AdamsBashMoulton3', 'AdamsBashforth2',
                   'AdamsBashforth3', 'AdamsBashforth4', 'AdaptiveResidual',
                   'Backward2Step', 'BackwardEuler', 'Euler',
                   'Heun', 'Leapfrog', 'LeapfrogFiltered',
                   'MidpointImplicit', 'MidpointIter', 'RK3',
                   'RK2', 'RK4', 'RKFehlberg',
                   'ThetaRule', 'Trapezoidal']
fail_solvers = [name for name in list_available_solvers() if name not in complex_solvers]

f = lambda u, t: 1./(t - 10. + 1j)
u0, time_points = 0, np.linspace(0., 20., 200)

print """ Scalar ODE with complex value
 u' =  1./(t - 10 + 1j) """

succeed_no = 0
for solver in complex_solvers:
    try:
        method = eval(solver)(f, atol=1e-6)
        method.set_initial_condition(u0)
        u,t = method.solve(time_points)
        success = np.allclose(u[-1].imag, -2*np.arctan(10))
        print solver, 'Succeed' if success else 'Fail'
    except:
        print 'Failed when solver is %s' % solver

"""
Unit tests for odespy solvers. Compatible with nose and pytest.
All available solvers are tested in four test problems:

 * scalar exponential decay
 * linear 2nd-order ODE for sine oscillations
 * van der Pol oscillator
 * scalar ODE with complex solution

Test method: assume an implementation is correct, record the difference
between its solution and an exact analytical solution, and compare
future simulations with the recorded difference (check that these are
within machine precision, here taken as 1E-14).
"""
# Author: Liwei Wang and Hans Petter Langtangen

import odespy
import numpy as np
import nose.tools as nt


# Dictionaries for each model problem, holding info about
# f, jac, and other parameters needed in the tests

Exponential = dict(
    help="Scalar ODE : Exponential u' = -u, u = exp(-t)",
    f=lambda u,t: -u,
    f_with_args=lambda u,t,a,b: a*u + b,
    f_args=[-1.,0.,],
    f_with_kwargs=lambda u,t,a=1.,b=0.: a*u + b,
    f_kwargs=dict(a=-1., b=0.,),
    # These solvers are not suitable for this ODE problem
    exceptions=['Lsodi', 'Lsodis', 'Lsoibt', 'MyRungeKutta',
                'MySolver', 'Lsodes', 'SymPy_odefun', 'AdaptiveResidual',
                'Radau5', 'Radau5Explicit', 'Radau5Implicit', 'EulerCromer'],
    time_points=np.linspace(0., 1., 10),
    terminate=lambda u,t,step_number: u[step_number] <= 0.4,
    stop_value=0.4,  # as defined in function terminate
    u0=1.,
    atol=1e-4,
    exact=lambda t: np.exp(-t),
    exact_diff=dict(
AdamsBashMoulton2=0.0003980947207471,
AdamsBashMoulton3=0.0005344126203188,
AdamsBashforth2=0.0018831513893969,
AdamsBashforth3=0.0003980947207471,
AdamsBashforth4=0.0005344126203188,
Backward2Step=0.0061598493269601,
BackwardEuler=0.0195709377009968,
BogackiShampine=0.0000229803241930,
CashKarp=0.0000000009049181,
CrankNicolson=0.0003752771263927,
Dop853=0.0000000000000006,
Dopri5=0.0000000012835836,
DormandPrince=0.0000000024189326,
Euler=0.0214400250568239,
Fehlberg=0.0000000900126869,
ForwardEuler=0.0214400250568239,
Heun=0.0008237438335759,
Leapfrog=0.0086696524501236,
LeapfrogFiltered=0.0348582551157020,
Lsoda=0.0000865549875511,
Lsodar=0.0000865549875511,
Lsode=0.0000876264677737,
MidpointImplicit=0.0003752771263927,
MidpointIter=0.0004456997096194,
RK2=0.0008237438335759,
RK3=0.0000229842151632,
RK4=0.0000005126486515,
RKC=0.0002698600815208,
RKF45=0.0000000061711490,
RKFehlberg=0.0000000898517009,
RungeKutta1=0.0214400250568239,
RungeKutta2=0.0008237438335759,
RungeKutta3=0.0000229881145921,
RungeKutta4=0.0000005126448920,
ThetaRule=0.0003752771263926,
Trapezoidal=0.0008237438335759,
Vode=0.0000226336696082,
lsoda_scipy=0.0000865549875511,
odefun_sympy=0.0000000000000000,
        ),
    )

Sine = dict(
    help="Sine, u'' = -u, u = sin(t)",
    f=lambda (u00,u11),t : [u11, -u00],
    jac=lambda (u00,u11),t: [[0., 1.], [- 1., 0.]],
    time_points=np.linspace(0., 1, 10),
    # These solvers are not suitable for this ODE problem
    exceptions=['Lsodi', 'Lsodes', 'Lsodis', 'Lsoibt', 'MyRungeKutta',
                'MySolver', 'AdaptiveResidual',
                'Radau5', 'Radau5Explicit', 'Radau5Implicit'],
    # EulerCromer works, despite the fact that the order of the ODEs is wrong
    terminate=lambda u,t,step_number: u[step_number][0] >= 0.5,
    stop_value=0.5,  # as defined in function terminate
    u0=[0., 1.],
    exact=lambda t: np.sin(t),
    exact_diff=dict(
AdamsBashMoulton2=0.0004527366536171,
AdamsBashMoulton3=0.0006643694031060,
AdamsBashforth2=0.0029177761328268,
AdamsBashforth3=0.0004527366536171,
AdamsBashforth4=0.0006643694031060,
Backward2Step=0.0096965380945701,
BackwardEuler=0.0472845069359756,
BogackiShampine=0.0000005327302367,
CashKarp=0.0000000022329617,
CrankNicolson=0.0005813318416031,
Dop853=0.0000000000000007,
Dopri5=0.0000000023675535,
DormandPrince=0.0000000041139426,
Euler=0.0454289746828882,
EulerCromer=0.0015802187182712,
Fehlberg=0.0000001419622074,
ForwardEuler=0.0454289746828882,
Heun=0.0012657016106976,
Leapfrog=0.0060023738285783,
LeapfrogFiltered=0.0501071455461246,
Lsoda=0.0000003263486278,
Lsodar=0.0000003263486278,
Lsode=0.0000014387203310,
MidpointImplicit=0.0005813318416031,
MidpointIter=0.0006938799224260,
RK2=0.0012657016106976,
RK3=0.0000451552290783,
RK4=0.0000007894844520,
RKC=0.0000116460267607,
RKF45=0.1108826285260906,
RKFehlberg=0.0000001416510547,
RungeKutta1=0.0454289746828882,
RungeKutta2=0.0012657016106976,
RungeKutta3=0.0000451493839151,
RungeKutta4=0.0000007894902545,
ThetaRule=0.0005813318416031,
Trapezoidal=0.0012657016106976,
Vode=0.0000019483075576,
lsoda_scipy=0.0000003263486278,
odefun_sympy=0.0000000000000000,
        ),
)

VanDerPol = dict(
    help="Van der Pol oscillator problem:  u'' = 3*(1 - u**2)*u' - u",
    f=lambda (u00,u11),t: [u11, 3.*(1 - u00**2)*u11 - u00],
    jac=lambda (u00,u11),t: [[0., 1.], [-6.*u00*u11 - 1., 3.*(1. - u00**2)]],
    u0=[2., 0.],
    # These solvers are not suitable for this ODE problem
    exceptions=['RKC', 'Lsodes', 'Leapfrog', 'Lsodi', 'Lsodis', 'Lsoibt',
                'MyRungeKutta', 'MySolver', 'AdaptiveResidual',
                'Radau5', 'Radau5Explicit', 'Radau5Implicit',
                'EulerCromer'],
    time_points=np.linspace(0., 1., 50),
    terminate=lambda u,t,step_number: u[step_number][0] <= 1.8,
    stop_value=1.8,  # as defined in function terminate
    exact_final=[1.7508022, -0.27099777],
    exact_final_diff=dict(
AdamsBashMoulton2=0.0375044598833543,
AdamsBashMoulton3=0.0375044958017816,
AdamsBashforth2=0.0375289430957795,
AdamsBashforth3=0.0375022466046535,
AdamsBashforth4=0.0375043238355941,
Backward2Step=0.0374948094483003,
BackwardEuler=0.0370754207976609,
BogackiShampine=0.0375036894156218,
CashKarp=0.0375036952618713,
CrankNicolson=0.0375047131527575,
Dop853=0.0375036952176235,
Dopri5=0.0375036952001755,
DormandPrince=0.0375036950592158,
Euler=0.0379365976159907,
Fehlberg=0.0375036945602063,
ForwardEuler=0.0379365976159907,
Heun=0.0375041929371673,
LeapfrogFiltered=0.0369386129565421,
Lsoda=0.0375036977160605,
Lsodar=0.0375036977160605,
Lsode=0.0375036965388451,
MidpointImplicit=0.0375047131527575,
MidpointIter=0.0374995994977019,
RK2=0.0375090516330185,
RK3=0.0375036046610746,
RK4=0.0375037080029861,
RKF45=0.0375036952592016,
RKFehlberg=0.0375036943807978,
RungeKutta1=0.0379365976159907,
RungeKutta2=0.0375090516330185,
RungeKutta3=0.0375036020516037,
RungeKutta4=0.0375037080029375,
ThetaRule=0.0374998399893818,
Trapezoidal=0.0375041929371673,
Vode=0.0375037017349515,
lsoda_scipy=0.0375036977160605,
odefun_sympy=0.0375036952176233,
        ),
    )

Complex = dict(
    help=" Scalar ODE with complex value  u' =  1./(t - 1 + 1j)",
    f=lambda u, t: 1./(t - 1. + 1j),
    complex_valued=True,
    time_points=np.linspace(0., 2., 15),
    u0=0.,
    exceptions=odespy.list_not_suitable_complex_solvers(),
    exact_final=(7.27645842122e-08-1.57079632742j),
    exact_final_diff=dict(
AdamsBashMoulton2=0.0003775451466630,
AdamsBashMoulton3=0.0003433008296798,
AdamsBashforth2=0.0090480564618154,
AdamsBashforth3=0.0014544213334521,
AdamsBashforth4=0.0005554178310184,
Backward2Step=0.0013132502301100,
BackwardEuler=0.0714487418239782,
CrankNicolson=0.0008503313468471,
Euler=0.0714488873119148,
ForwardEuler=0.0714488873119148,
Heun=0.0017006724680425,
Leapfrog=0.0034008405626300,
LeapfrogFiltered=0.1060876548770292,
MidpointImplicit=0.0008503313468471,
MidpointIter=0.0017006724680425,
RK2=0.0008503313468471,
RK3=0.0000000728375602,
RK4=0.0000000728375602,
RKFehlberg=0.0000000648779403,
ThetaRule=0.0017006724680423,
Trapezoidal=0.0017006724680425,
Vode=0.0000005669366841,
        ),
    )

def _run_test_problems(problem):
    """Main function for executing a unit test."""
    # Test all available solvers for a problem, except those listed
    # as exceptions in the problem dict
    methods = [method for method in odespy.list_available_solvers()
               if method not in problem['exceptions']]

    # When constructing a method, supply all keys in problem except
    # the ones below as **kwargs arguments to the constructor
    special_names = ('f', 'time_points', 'u0', 'exact', 'exact_final',
                     'help', 'terminate', 'exceptions', 'stop_value',
                     'f_with_args', 'f_args', 'f_with_kwargs',
                     'f_kwargs')
    kwargs = {key: problem[key] for key in problem
              if key not in special_names}
    print problem.get('help', '')
    for method in methods:

        solver = eval('odespy.%s' % method)(problem['f'], **kwargs)
        print 'Testing %s' % method,
        solver.set_initial_condition(problem['u0'])
        u, t = solver.solve(problem['time_points'])

        # Test additional versions of the right-hand side function
        if 'f_with_args' in problem and 'f_args' in problem:
            #print 'Testing %s with f_args' % method
            solver = eval('odespy.%s' % method)(
                problem['f_with_args'], f_args=problem['f_args'], **kwargs)
            solver.set_initial_condition(problem['u0'])
            u, t = solver.solve(problem['time_points'])
        if 'f_with_kwargs' in problem and 'f_kwargs' in problem:
            #print 'Testing %s with f_kwargs' % method
            solver = eval('odespy.%s' % method)(
                problem['f_with_kwargs'], f_kwargs=problem['f_kwargs'],**kwargs)
            solver.set_initial_condition(problem['u0'])
            u, t = solver.solve(problem['time_points'])

        # Compare with exact values
        if 'exact_diff' in problem:
            t = np.asarray(t)
            u = np.asarray(u)
            exact = problem['exact'](t)
            if len(u.shape) == 1:  # make 2d-array for scalar ODE
                u = u.reshape(len(u), 1)
            diff = np.abs(u[:,0] - exact).max()
            if method in problem['exact_diff']:
                nt.assert_almost_equal(diff, problem['exact_diff'][method],
                                       delta=1E-14)
                print '...ok'
            else:
                pass
                print ' no exact diff available for comparison'
        if 'exact_final_diff' in problem:
            u_final = np.asarray(u[-1])
            exact_final = np.asarray(problem['exact_final'])
            diff = np.abs(u_final - exact_final).max()
            if method in problem['exact_final_diff']:
                nt.assert_almost_equal(diff, problem['exact_final_diff'][method],
                                       delta=1E-14)
                print '...ok'
            else:
                print ' no exact final diff available for comparison'

def test_exponentinal():
    _run_test_problems(Exponential)

def test_sine():
    _run_test_problems(Sine)

def test_vanderpol():
    _run_test_problems(VanDerPol)

def test_complex():
    _run_test_problems(Complex)

def test_EulerCromer_1dof():
    """Plain u'' + u = 0."""
    solver = odespy.EulerCromer(lambda (v, x), t: [-x, v])
    #solver = odespy.EulerCromer(lambda u, t: [-u[1], u[0]])
    solver.set_initial_condition([0, 1])
    P = 2*np.pi
    N = 60
    dt = P/N
    num_periods = 8
    T = num_periods*P
    time_points = np.linspace(0, T, num_periods*N+1)
    u, t = solver.solve(time_points)
    x = u[:,1]
    v = u[:,0]
    x_exact = lambda t: np.cos(t)
    x_e = x_exact(t)
    #plot(t, x, t, x_e, legend=('EC', 'exact'))
    print 'Testing EulerCromer',
    # Test difference in final value
    if N == 60:
        diff_exact = 0.0014700112828
        diff_EC = abs(x[-1] - x_e[-1])
        tol = 1E-14
        assert abs(diff_exact - diff_EC) < tol, \
               'diff_exact=%g, diff_EulerCromer=%g' % (diff_exact, diff_EC)
        print '...ok'

def test_terminate():
    problem = Exponential
    solver = odespy.RKFehlberg(problem['f'])
    solver.set_initial_condition(problem['u0'])
    u, t = solver.solve(problem['time_points'], terminate=problem['terminate'])
    exact_diff = 0.0321206486802586
    diff = abs(problem['stop_value']-u[-1])
    print 'Testing RKFehlberg with terminate function',
    nt.assert_almost_equal(diff, exact_diff, delta=1E-14)
    print '...ok'

def test_switch_to():
    problem = Exponential
    solver = odespy.RKFehlberg(problem['f'])
    solver.set_initial_condition(problem['u0'])
    u, t = solver.solve(problem['time_points'])
    solver_new = solver.switch_to(odespy.Lsode)
    u2, t2 = solver_new.solve(problem['time_points'])
    diff = np.abs(u - u2).max()
    exact_diff = 0.0000015284633990
    print 'Testing switch_to from RKFehlberg to Lsode',
    nt.assert_almost_equal(diff, exact_diff, delta=1E-14)
    print '...ok'


if __name__ == '__main__':
    test_exponentinal()
    test_sine()
    test_vanderpol()
    test_complex()
    test_EulerCromer_1dof()
    test_switch_to()
    test_terminate()

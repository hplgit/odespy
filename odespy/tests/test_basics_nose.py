"""Unit test for odespy solvers."""
# Author: Liwei Wang and Hans Petter Langtangen

import odespy
import scitools.std as st
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
    exact=lambda t: np.exp(-t),
    # These solvers are not suitable for this ODE problem
    exceptions=['Lsodi', 'Lsodis', 'Lsoibt', 'MyRungeKutta',
                'MySolver', 'Lsodes', 'SymPy_odefun', 'AdaptiveResidual',
                'Radau5', 'Radau5Explicit', 'Radau5Implicit'],
    time_points=np.linspace(0., 1., 10),
    terminate=lambda u,t,step_number: u[step_number] <= 0.4,
    stop_value=0.4,  # as defined in function terminate
    u0=1.,
    atol=1e-4)

Sine = dict(
    help="Sine, u'' = -u, u = sin(t)",
    f=lambda (u00,u11),t : [u11, -u00],
    jac=lambda (u00,u11),t: [[0., 1.], [- 1., 0.]],
    exact=lambda t: np.sin(t),
    time_points=np.linspace(0., 1, 10),
    # These solvers are not suitable for this ODE problem
    exceptions=['Lsodi', 'Lsodes', 'Lsodis', 'Lsoibt', 'MyRungeKutta',
                'MySolver', 'AdaptiveResidual',
                'Radau5', 'Radau5Explicit', 'Radau5Implicit'],
    terminate=lambda u,t,step_number: u[step_number][0] >= 0.5,
    stop_value=0.5,  # as defined in function terminate
    u0=[0., 1.])

VanDerPol = dict(
    help="Van der Pol oscillator problem:  u'' = 3*(1 - u**2)*u' - u",
    f=lambda (u00,u11),t: [u11, 3.*(1 - u00**2)*u11 - u00],
    jac=lambda (u00,u11),t: [[0., 1.], [-6.*u00*u11 - 1., 3.*(1. - u00**2)]],
    u0=[2., 0.],
    # These solvers are not suitable for this ODE problem
    exceptions=['RKC', 'Lsodes', 'Leapfrog', 'Lsodi', 'Lsodis', 'Lsoibt',
                'MyRungeKutta', 'MySolver', 'AdaptiveResidual',
                'Radau5', 'Radau5Explicit', 'Radau5Implicit'],
    time_points=np.linspace(0., 1., 50),
    terminate=lambda u,t,step_number: u[step_number][0] <= 1.8,
    stop_value=1.8,  # as defined in function terminate
    exact_final=[1.7508022, -0.27099777])

Complex = dict(
    help=" Scalar ODE with complex value  u' =  1./(t - 1 + 1j)",
    f=lambda u, t: 1./(t - 1. + 1j),
    complex_valued=True,
    time_points=np.linspace(0., 2., 15),
    u0=0.,
    exceptions=odespy.list_not_suitable_complex_solvers(),
    exact_final=(7.27645842122e-08-1.57079632742j))

def _run_test_problems(problem):
    method_no = 0      # Number of successful methods
    methods = [method for method in odespy.list_available_solvers()
               if method not in problem['exceptions']]

    special_names = ('f', 'time_points', 'u0', 'exact', 'exact_final',
                     'help', 'terminate', 'exceptions', 'stop_value',
                     'f_with_args', 'f_args', 'f_with_kwargs',
                     'f_kwargs')
    kwargs = {key: problem[key] for key in problem
              if key not in special_names}
    print problem.get('help', '')
    for method in methods:

        solver = eval('odespy.%s' % method)(problem['f'], **kwargs)
        print 'Testing %s' % method
        solver.set_initial_condition(problem['u0'])
        u, t = solver.solve(problem['time_points'])

        if 'f_with_args' in problem and 'f_args' in problem:
            print 'Testing %s with f_args' % method
            solver = eval('odespy.%s' % method)(
                problem['f_with_args'], f_args=problem['f_args'], **kwargs)
            solver.set_initial_condition(problem['u0'])
            u, t = solver.solve(problem['time_points'])
        if 'f_with_kwargs' in problem and 'f_kwargs' in problem:
            print 'Testing %s with f_kwargs' % method
            solver = eval('odespy.%s' % method)(
                problem['f_with_kwargs'], f_kwargs=problem['f_kwargs'],**kwargs)
            solver.set_initial_condition(problem['u0'])
            u, t = solver.solve(problem['time_points'])

        # Compare if exact values are specified
        if 'exact' in problem:
            exact = problem['exact'](np.asarray(t))
            u = np.asarray(u)
            if len(u.shape) == 1:  # make 2d-array for scalar ODE
                u = u.reshape(len(u), 1)
            diff = np.abs(u[:,0] - exact).max()
            nt.assert_almost_equal(diff, 0, delta=0.2)
        elif 'exact_final' in problem:
            u_final, exact_final = u[-1], problem['exact_final']
            u_final = np.asarray(u_final)
            exact_final = np.asarray(problem['exact_final'])
            diff = np.abs(u_final - exact_final).max()
            nt.assert_almost_equal(diff, 0, delta=0.2)

def test_exponentinal():
    _run_test_problems(Exponential)

def test_sine():
    _run_test_problems(Sine)

def test_vanderpol():
    _run_test_problems(VanDerPol)

def test_complex():
    _run_test_problems(Complex)

def test_terminate():
    problem = Exponential
    solver = odespy.RKFehlberg(problem['f'])
    solver.set_initial_condition(problem['u0'])
    u, t = solver.solve(problem['time_points'], terminate=problem['terminate'])
    nt.assert_almost_equal(u[-1], problem['stop_value'], delta=0.5)

def test_switch_to():
    problem = Exponential
    solver = odespy.RKFehlberg(problem['f'])
    solver.set_initial_condition(problem['u0'])
    u, t = solver.solve(problem['time_points'])
    solver_new = solver.switch_to(odespy.Lsode)
    u2, t2 = solver.solve(problem['time_points'])
    diff = np.abs(u - u2).max()
    nt.assert_almost_equal(diff, 0, delta=0.1)


if __name__ == '__main__':
    test_switch_to()
    test_terminate()
    test_exponentinal()
    test_sine()
    test_vanderpol()
    test_complex()

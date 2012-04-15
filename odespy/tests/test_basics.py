"""Unit test for odespy solvers."""
# Author: Liwei Wang

import odespy
import scitools.std as st
import numpy as np

from numpy.testing import assert_array_almost_equal, TestCase, \
    run_module_suite, assert_almost_equal


class TestBasics(TestCase):
    """
    Test basic ODE problems for all possible solvers
    """

    def test_exponentinal(self):
        self._run_test_problems(Exponential)

    def test_sine(self):
        self._run_test_problems(Sine)

    def test_vanderpol(self):
        self._run_test_problems(VanDerPol)

    def test_complex(self):
        self._run_test_problems(Complex)

    def test_switch_to(self):
        for problem in [Exponential, Sine]:
            self._load_problem(problem)
            solver = odespy.RKFehlberg(self.f, **self.kwargs)
            solver.set_initial_condition(self.u0)
            u, t = solver.solve(self.time_points)

            solver_new = solver.switch_to(odespy.Lsode)
            u_new, t_new = solver_new.solve(self.time_points)

            assert_array_almost_equal(\
                u, u_new,
                err_msg='''
       Failed for switch from RKFehlberg to Lsode with problem %s''' \
                    % self.help,
                decimal=2, verbose=False)

    def test_terminate(self):
        for problem in [Exponential, Sine, VanDerPol]:
            self._load_problem(problem)

            solver = odespy.RKFehlberg(self.f, **self.kwargs)
            solver.set_initial_condition(self.u0)
            u, t = solver.solve(self.time_points, terminate=self.terminate)

            u_stop = u[-1][0] if len(u.shape) == 2 else u[-1]
            assert_almost_equal(u_stop, self.stop_value, verbose=True,
                                decimal=1)

    def test_f_args(self):
        self._load_problem(Exponential)

        solver = odespy.RKFehlberg(self.f_with_args,
                                   f_args=self.f_args,
                                   **self.kwargs)
        solver.set_initial_condition(self.u0)
        u, t = solver.solve(self.time_points)

        exact = self.exact(t)
        assert_array_almost_equal(\
            u, exact,
            err_msg='Failed with f_args',
            decimal=2, verbose=True)

    def test_f_kwargs(self):
        self._load_problem(Exponential)

        solver = odespy.RKFehlberg(self.f_with_kwargs,
                                   f_kwargs=self.f_kwargs,
                                   **self.kwargs)
        solver.set_initial_condition(self.u0)
        u, t = solver.solve(self.time_points)

        exact = self.exact(t)
        assert_array_almost_equal(\
            u, exact,
            err_msg='Failed with f_args',
            decimal=2, verbose=True)



    def _load_problem(self, problem):
        """
        Load keys in the `problem` dict into two categories:
        attributes (like ``self.f``, ``self.time_points``) and
        a dict kwargs of optional parameters which can be
        transferred as is to the solver class (a la
        ``solver(self.f, **self.kwargs)``).
        """
        special_names = ('f', 'time_points', 'u0', 'exact', 'exact_final',
                         'help', 'terminate', 'exceptions', 'stop_value',
                         'f_with_args', 'f_args', 'f_with_kwargs',
                         'f_kwargs')
        self.kwargs = {}
        for name in problem:
            if name in special_names:
                setattr(self, name, problem[name])
            else:  # optional parameters
                self.kwargs[name] = problem[name]


    def _run_test_problems(self, problem):
        self._load_problem(problem)

        method_no = 0      # Number of successful methods
        methods = [method for method in odespy.list_available_solvers()
                   if method not in self.exceptions]

        print self.help
        for method in methods:
            print 'Testing %s' % method
            # Start up integration
            solver = eval('odespy.%s' % method)(self.f, **self.kwargs)
            solver.set_initial_condition(self.u0)
            u, t = solver.solve(self.time_points)

            # Compare if exact values are specified
            if hasattr(self, 'exact'):
                exact = self.exact(np.asarray(t))
                u = np.asarray(u)
                if len(u.shape) == 1:  # make 2d-array for scalar ODE
                    u = u.reshape(len(u), 1)
                assert_array_almost_equal(\
                    u[:,0], exact,
                    err_msg='Failed with method %s' % method,
                    decimal=1, verbose=False)
            elif hasattr(self, 'exact_final'):
                u_final, exact_final = u[-1], self.exact_final
                if not np.iterable(u_final):
                    u_final = np.asarray(u_final)
                    exact_final = np.asarray(self.exact_final)
                try:
                    assert_array_almost_equal(\
                        u_final, exact_final,
                        err_msg='Failed with result of method %s' % method,
                        decimal=1, verbose=True)
                except Exception, e:
                    print e
                    print 'Running method', method, 'for', problem
                    raise e

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
                'MySolver', 'Lsodes', 'SymPy_odefun'],
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
                'MySolver', 'AdaptiveResidual'],
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
                'MyRungeKutta', 'MySolver', 'AdaptiveResidual'],
    time_points=np.linspace(0., 1., 50),
    terminate=lambda u,t,step_number: u[step_number][0] <= 1.8,
    stop_value=1.8,  # as defined in function terminate
    exact_final=[1.7508022, -0.27099777])

Complex = dict(
    help=" Scalar ODE with complex value  u' =  1./(t - 1 + 1j)",
    f=lambda u, t: 1./(t - 1. + 1j),
    time_points=np.linspace(0., 2., 15),
    u0=0.,
    exceptions=['BogackiShampine', 'CashKarp', 'Dop853', 'Dopri5',
                'DormandPrince', 'Fehlberg', 'RungeKutta1',
                'Lsoda', 'Lsodar', 'Lsode', 'Lsodes', 'Lsodi', 'Lsodis',
                'Lsoibt', 'MyRungeKutta', 'MySolver', 'RKC', 'RKF45',
                'RungeKutta2', 'RungeKutta3', 'RungeKutta4',
                'odefun_sympy', 'Vode', 'lsoda_scipy',
                'Radau5', 'Radau5Explicit', 'Radau5Implicit'],
    exact_final=(7.27645842122e-08-1.57079632742j))

if __name__ == "__main__":
    run_module_suite()


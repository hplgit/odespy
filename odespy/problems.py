import numpy as np
from solvers import compile_f77

class Problem:
    '''
    The user derives a subclass and implements the right-hand side
    function::

        def f(self, u, t):
            """Python function for the right-hand side of the ODE u'=f(u,t)."""

    Optionally, the Jacobian can be computed::

        def jac(self, u, t):
            """Python function for the Jacobian of the right-hand side: df/du."""
    Some solvers also allow constraings (DAE problems)::

        def constraints(self, u, t):
            """Python function for additional constraints: g(u,t)=0."""

    Some problem classes will also define command-line arguments::

        def define_command_line_arguments(self, parser):
            """
            Initialize an argparse object for reading command-line
            option-value pairs. `parser` is an ``argparse`` object.
            """

    Other functions, for which there are default implementations
    in the superclass, are ``u_exact`` for returning the exact
    solution (if known), ``verify`` for checking that a numerical
    solution is correct (within a tolerance of the analytical solution),

    See the tutorial for examples on subclasses of ``Problem``.
    '''
    stiff = False    # classification of the problem is stiff or not
    complex_ = False # True if f(u,t) is complex valued
    not_suitable_solvers = []  # list solvers that should be be used
    short_description = ''     # one-line problem description

    def __init__(self):
        pass

    def __contains__(self, attr):
        """Return True if attr is a method in instance self."""
        return hasattr(self, attr) and callable(getattr(self,attr))

    def terminate(self, u, t, step_number):
        """Default terminate function, always returning False."""
        return False

    def get_initial_condition(self):
        """
        Return vector of initial conditions.
        Not necessary to implement in subclass if the
        initial condition is stored in self.U0
        (this method then returns that condition).
        """
        if hasattr(self, 'U0'):
            return self.U0
        else:
            raise NotImplementedError(
                'class %s must implement get_initial_condition' %
                self.__class__.__name__)

    def default_parameters(self):
        """
        Compute suitable time_points, atol/rtol, etc. for the
        particular problem. Useful for quick generation of test
        cases, demos, unit tests, etc.
        """
        return {}

    def u_exact(self, t):
        """
        Implementation of the exact solution.
        Return None if no exact solution is available.
        """
        return None

    def verify(self, u, t, atol=None, rtol=None):
        """
        Return True if u at time points t coincides with an exact
        solution within the prescribed tolerances. If one of the
        tolerances is None, return max computed error (infinity norm).
        Return None if the solution cannot be verified.
        """
        u_e = self.u_exact(t)
        if u_e is None:
            return None

        if atol is None or rtol is None:
            if not (u.shape == u_e.shape):
                raise ValueError('u has shape %s and u_e has %s' %
                                 (u.shape, u_e.shape))
            return np.abs(u - u_e).max()
        else:
            return np.allclose(u, u_e, rtol, atol)

    # subclasses may implement computation of derived
    # quantities, e.g., energy(u, t) etc

def convergence(u, t, u_ref=None, t_ref=None):
    """
    Given a series of solutions and corresponding t arrays in
    `u` and `t`, use the analytical solution or a reference solution
    in `u_ref` and `t_ref` to compute errors. Compute pairwise
    convergence rates for errors on consecutive time meshes.
    """
    raise NotImplementedError

class Linear1(Problem):
    """
    Linear solution of a nonlinear ODE,
    which is normally exactly reproduced by a numerical method.
    """
    short_description = "u' = a + (u - a*t - b)**c, u(0)=b"

    def __init__(self, a=1, b=0, c=2):
        self.a, self.b, self.c = a, b, c
        self.U0 = b
        self.not_suitable_solvers = [
            'AdaptiveResidual', 'Lsodi', 'Lsodis', 'Lsoibt']
        if self.a < 0:
            self.not_suitable_solvers += ['lsoda_scipy', 'odefun_sympy']

    def f(self, u, t):
        return self.a + (u - self.a*t - self.b)**self.c

    def jac(self, u, t):
        return self.c*(u - self.a*t - self.b)**(self.c-1)

    def u_exact(self, t):
        return self.a*t + self.b

    def default_parameters(self):
        return {'time_points': np.linspace(0, 2./self.a, 3),
                'max_iter': 10, 'eps_iter': 1E-6}

class Linear2(Linear1):
    """
    Linear solution of nonlinear 2x2 system,
    which is normally exactly reproduced by a numerical method.
    """
    short_description = "2x2 nonlinear system with linear solution"

    def __init__(self, a=1, b=0, c=2):
        Linear1.__init__(self, a, b, c)
        self.U0 = [self.b, self.b]

    def f(self, u, t):
        return [self.a + (u[1] - self.a*t - self.b)**self.c,
                self.a + (u[0] - self.a*t - self.b)**self.c]

    def jac(self, u, t):
        return [[0, self.c*(u[1] - self.a*t - self.b)**(self.c-1)],
                [self.c*(u[0] - self.a*t - self.b)**(self.c-1), 0]]

    def spcrad(self, u, t):
        J = self.jac(u, t)
        return np.abs(np.linalg.eigvals(J)).max()

    def u_exact(self, t):
        func = self.a*t + self.b
        return np.array([func, func]).transpose()


class Linear2t(Linear2):
    """
    Linear solution of trivial 2x2 system,
    which is normally exactly reproduced by a numerical method.
    """
    short_description = "2x2 trivial system with linear solution"

    def f(self, u, t):
        return [self.a,
                self.a]

    def jac(self, u, t):
        return [[0, 0],
                [0, 0]]

    def spcrad(self, u, t):
        return 0.0



class Exponential(Problem):
    short_description = "Exponential solution: u' = a*u + b, u(0)=A"

    def __init__(self, a=1, b=0, A=1):
        self.a, self.b, self.U0 = float(a), b, A
        if abs(a) < 1E-15:
            raise ValueError('a=%g is too small' % a)
        self.not_suitable_solvers = [
            'Lsodi', 'Lsodis', 'Lsoibt', 'MyRungeKutta',
            'MySolver', 'Lsodes', 'SymPy_odefun']

    def f(self, u, t):
        return self.a*u + self.b

    def jac(self, u, t):
        return self.a

    def f_with_args(self, u, t, a, b):
        return a*u + b

    def f_with_kwargs(self, u, t, a=1, b=0):
        return a*u + b

    def jac_with_args(self, u, t, a, b):
        return a

    def jac_with_kwargs(self, u, t, a=1, b=0):
        return a

    def u_exact(self, t):
        a, b, A = self.a, self.b, self.A
        return np.exp(a*t)*(A + b/a) -  b/a

    def default_parameters(self):
        T = 5/abs(self.a)
        d = {'time_points': np.linspace(0, T, 31),
             'atol': 1E-4, 'rtol': 1E-3}
        if self.a < 0:
            d['terminate'] = lambda u, t, step_no: u[step_no] < self.A/20.
        return d

class Logistic(Problem):
    short_description = "Logistic equation"

    def __init__(self, a, R, A):
        self.a = a
        self.R = R
        self.U0 = A

    def f(self, u, t):
        a, R = self.a, self.R  # short form
        return a*u*(1 - u/R)

    def jac(self, u, t):
        a, R = self.a, self.R  # short form
        return a*(1 - u/R) + a*u*(1 - 1./R)

    def u_exact(self, t):
        a, R, U0 = self.a, self.R, self.U0  # short form
        return R*U0*np.exp(a*t)/(R + U0*(np.exp(a*t) - 1))

class Gaussian1(Problem):
    short_description = "1 + Gaussian function as solution"

    def __init__(self, c=2.0, s=1.0):
        self.c, self.s = c, float(s)
        self.U0 = self.u_exact(0)

    def f(self, u, t):
        return -((t-self.c)/self.s**2)*(u-1)

    def jac(self, u, t):
        return -(t-self.c)/self.s**2

    def u_exact(self, t):
        return 1 + np.exp(-0.5*((t-self.c)/self.s)**2)

class Gaussian0(Problem):
    short_description = "Gaussian function as solution"

    def __init__(self, c=2.0, s=1.0):
        self.c, self.s = c, float(s)
        self.U0 = self.u_exact(0)

    def f(self, u, t):
        return -((t-self.c)/self.s**2)*u

    def jac(self, u, t):
        return -(t-self.c)/self.s**2

    def u_exact(self, t):
        return np.exp(-0.5*((t-self.c)/self.s)**2)


def default_oscillator(P, resolution_per_period=20):
    n = 4.5
    r = resolution_per_period # short form
    print 'default', P, n*P
    tp = np.linspace(0, n*P, r*n+1)
    # atol=rtol since u approx 1, atol level set at the
    # error RK4 produces with 20 steps per period, times 0.05
    error_level = 3E-3
    atol = rtol = 0.05*error_level
    # E = C*dt**q
    return dict(time_points=tp, atol=atol, rtol=rtol)

class LinearOscillator(Problem):
    short_description = "Linear oscillator: m*u'' + k*u = F(t)"

    def __init__(self, m=1, k=1, x0=1, v0=0, F=None):
        self.U0 = [x0, v0]
        self.m, self.k = float(m), k
        self.F = F

    def f(self, u, t):
        x, v = u
        k, m = self.k, self.m
        F = 0 if self.F is None else self.F(t)
        return [v, -k/m*x + F]

    def jac(self, u, t):
        x, v = u
        k, m = self.k, self.m
        return [[0, 1], [-k/m, 0]]

    def u_exact(self, t):
        if self.F is None:
            k, m = self.k, self.m
            w = np.sqrt(k/m)
            x = self.U0[0]*np.cos(w*t) + self.U0[1]*np.sin(w*t)
            v = -w*self.U0[0]*np.sin(w*t) + w*self.U0[1]*np.cos(w*t)
            return np.array([x, v]).transpose()
        else:
            return None

    def default_parameters(self, resolution_per_period=20):
        w = self.k/self.m
        P = 2*pi/np.sqrt(w)  # period
        return default_oscillator(P, resolution_per_period)

    def kinetic_energy(self, u):
        v = u[:,1]
        return 0.5*self.m*v**2

    def potential_energy(self, u):
        x = u[:,0]
        return 0.5*self.k*x**2

# offer analytical solution with damping, also with sin/cos excitation

class ComplexOscillator(Problem):
    """
    Harmonic oscillator (u'' + w*u = 0) expressed as a complex
    ODE: u' = i*w*u.
    """
    short_description = "Complex oscillator: u' = i*w*u"

    def __init__(self, w=1, U0=[1, 0]):
        self.w = w
        self.not_suitable_solvers = list_not_suitable_complex_solvers()

    def f(self, u, t):
        return 1j*self.w*u

    def jac(self, u, t):
        return 1j*self.w

    def u_exact(self, t):
        return np.exp(1j*self.w*t)

    def default_parameters(self, resolution_per_period=20):
        P = 2*pi/np.sqrt(self.w)  # period
        return default_oscillator(P, resolution_per_period)

class StiffSystem2x2(Problem):
    short_description = "Potentially stiff 2x2 system u' = 1/eps*u"
    stiff = True

    def __init__(self, eps=1E-2):
        self.eps = eps
        self.U0 = [1, 1]
        if 0.2 <= eps <= 5:
            self.stiff = False

    def f(self, u, t):
        return [-u[0], -1./self.eps*u[1]]

    def jac(self, u, t):
        return [[-1, 0], [0, -1./self.eps]]

    def u_exact(self, t):
        return np.array([np.exp(-t), np.exp(-1./self.eps*t)]).transpose()


class VanDerPolOscillator(Problem):
    """
    Classical van der Pool oscillator:

    .. math:: y'' = \mu (1 - y^2) y' - y

    with initial conditions :math:`y(0)=2, y'(0)=1`.
    The equation is rewritten as a system

    .. math::
             u_0' &= u_1, \\
             u_1' &= \mu (1-u_0^2)u_1 - u_0

    with a Jacobian

    .. math::
             \left(\begin{array}{cc}
             0 & 1\\
             -2\mu u_0 - 1 & \mu (1-u_0^2)
             \end{array}\right)
    """
    short_description = "Van der Pol oscillator"

    def __init__(self, U0=[2, 1], mu=3., f77=False):
        self.U0 = U0
        self.mu = mu

        # Compile F77 functions
        if f77:
            self.f_f77, self.jac_f77_radau5, self.jac_f77_lsode = \
                        compile_f77([self.str_f_f77(),
                                     self.str_jac_f77_radau5(),
                                     self.str_jac_f77_lsode()])

    def f(self, u, t):
        u_0, u_1 = u
        mu = self.mu
        return [u_1, mu*(1 - u_0**2)*u_1 - u_0]

    def jac(self, u, t):
        u_0, u_1 = u
        mu = self.mu
        return [[0., 1.],
                [-2*mu*u_0*u_1 - 1, mu*(1 - u_0**2)]]

    def f_with_args(self, u, t, mu):
        u_0, u_1 = u
        return [u_1, mu*(1 - u_0**2)*u_1 - u_0]

    def jac_with_args(self, u, t, mu):
        u_0, u_1 = u
        return [[0., 1.],
                [-2*mu*u_0*u_1 - 1, mu*(1 - u_0**2)]]

    def f_with_kwargs(self, u, t, mu=3):
        u_0, u_1 = u
        return [u_1, mu*(1 - u_0**2)*u_1 - u_0]

    def jac_with_kwargs(self, u, t, mu=3):
        u_0, u_1 = u
        return [[0., 1.],
                [-2*mu*u_0*u_1 - 1, mu*(1 - u_0**2)]]


    def str_f_f77(self):
        """Return f(u,t) as Fortran source code string."""
        return """
      subroutine f_f77(neq, t, u, udot)
Cf2py intent(hide) neq
Cf2py intent(out) udot
      integer neq
      double precision t, u, udot
      dimension u(neq), udot(neq)
      udot(1) = u(2)
      udot(2) = %g*(1 - u(1)**2)*u(2) - u(1)
      return
      end
""" % self.mu

    def str_jac_f77_fadau5(self):
        """Return Jacobian Fortran source code string for Radau5."""
        return """
      subroutine jac_f77_radau5(neq,t,u,dfu,ldfu,rpar,ipar)
Cf2py intent(hide) neq,rpar,ipar
Cf2py intent(in)   t,u,ldfu
Cf2py intent(out) dfu
      integer neq,ipar,ldfu
      double precision t,u,dfu,rpar
      dimension u(neq),dfu(ldfu,neq),rpar(*),ipar(*)
      dfu(1,1) = 0
      dfu(1,2) = 1
      dfu(2,1) = -2*%g*u(1)*u(2) - 1
      dfu(2,2) = %g*(1-u(1)**2)
      return
      end
""" % (self.mu, self.mu)

    def str_jac_f77_lsode_dense(self):
        """Return Fortran source for dense Jacobian matrix in LSODE format."""
        return """
      subroutine jac_f77(neq, t, u, ml, mu, pd, nrowpd)
Cf2py intent(hide) neq, ml, mu, nrowpd
Cf2py intent(out) pd
      integer neq, ml, mu, nrowpd
      double precision t, u, pd
      dimension u(neq), pd(nrowpd,neq)
      pd(1,1) = 0
      pd(1,2) = 1
      pd(2,1) = -2*%g*u(1)*u(2) - 1
      pd(2,2) = %g*(1 - u(1)**2)
      return
      end
""" % (self.mu, self.mu)

    def u_exact(self, t):
        if abs(self.mu) < 1E-14:
            x = self.U0[0]*np.cos(t) + self.U0[1]*np.sin(t)
            v = -self.U0[0]*np.sin(t) + self.U0[1]*np.cos(t)
            return np.array([x, v]).transpose()
        else:
            return None

    def default_parameters(self, resolution_per_period=20):
        # Period = 2*pi for mu=0
        P = 2*np.pi
        return default_oscillator(P, resolution_per_period)



def tester(problems, methods, time_points=None, compare_tol=1E-4,
           solver_prm={}):
    """
    `problems` is a list of Problem subclass instances made ready.
    `methods` is a list of strings holding the method names.
    """
    import odespy
    results = {}
    error_msg = {}
    for problem in problems:
        pname = problem.__class__.__name__
        print 'problem ', pname
        methods4problem = [method for method in methods
                           if method not in problem.not_suitable_solvers]
        defaults = problem.default_parameters()
        if time_points is None:
            time_points = defaults['time_points']
        results[pname] = {'t': time_points}
        error_msg[pname] = {}

        this_solver_prm = solver_prm.copy()
        names_in_problem = 'jac', 'spcrad', 'u_exact',
        for name in names_in_problem:
            # Set parameter if problem has it and user has
            # not specified it in solver_prm
            if name in problem and name not in solver_prm:
                this_solver_prm[name] = getattr(problem, name)
        names_in_defaults = 'atol', 'rtol', 'max_iter', 'eps_iter'
        for name in names_in_defaults:
            if name in defaults and name not in solver_prm:
                this_solver_prm[name] = defaults[name]

        for method in methods4problem:
            print '  testing', method,
            solver = eval('odespy.'+method)(problem.f)

            # Important to set parameters before setting initial cond.
            solver.set(**this_solver_prm)

            solver.set_initial_condition(problem.get_initial_condition())

            u, t = solver.solve(time_points, problem.terminate)

            error = problem.verify(u, t)
            if error is not None:
                results[pname][method] = (u, error)
                print error,
                if error > compare_tol:
                    print 'WARNING: tolerance %.0E exceeded' % compare_tol,
            else:
                results[pname][method] = (u,)
            print
    return results

"""

Radau5 {'rtol': 0.01, 'atol': 0.01, 'max_iter': 10, 'jac': <bound method Linear2.jac of <problems.Linear2 instance at 0x1fdd8c0>>, 'eps_iter': 1e-06, 'spcrad': <bound method Linear2.spcrad of <problems.Linear2 instance at 0x1fdd8c0>>}
1.7763568394e-15
Radau5Explicit {'rtol': 0.01, 'atol': 0.01, 'max_iter': 10, 'jac': <bound method Linear2.jac of <problems.Linear2 instance at 0x1fdd8c0>>, 'eps_iter': 1e-06, 'spcrad': <bound method Linear2.spcrad of <problems.Linear2 instance at 0x1fdd8c0>>}
1.7763568394e-15
Radau5Implicit {'rtol': 0.01, 'atol': 0.01, 'max_iter': 10, 'jac': <bound method Linear2.jac of <problems.Linear2 instance at 0x1fdd8c0>>, 'eps_iter': 1e-06, 'spcrad': <bound method Linear2.spcrad of <problems.Linear2 instance at 0x1fdd8c0>>}

    # Compare with exact solution
    # Compare solutions within tolerance
    failure = {}
    for pname in problems:
        failure[pname] = {}
        for i, mname in enumeate(methods):
            if i == 0:
                reference_solution, t = results[(mname, pname)]
                reference_name = mname
            else:
                u = results[(mname, pname)][0]
                r = np.allclose(reference_solution, u,
                                compare_tol, 0.1*compare_tol)
                if r:
                    failure[pname][mname] = False
                else:
                    failure[pname][mname] = np.abs(reference_solution - u).max()
    return failure
"""

import numpy as np
from solvers import compile_f77

class Problem:
    stiff = False
    complex_ = False

    def __init__(self):
        pass

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

    def f(self, u, t):
        """Python function for the right-hand side of the ODE u'=f(u,t)."""
        raise NotImpelementedError

    def jac(self, u, t):
        """Python function for the Jacobian of the right-hand side: df/du."""
        raise NotImpelementedError

    def constraints(self, u, t):
        """Python function for additional constraints: g(u,t)=0."""
        raise NotImpelementedError

    def define_command_line_arguments(self, parser):
        """Initialize an argparse object for reading command-line
        option-value pairs."""
        raise NotImpelementedError

    def u_exact(self, t):
        """
        Implementation of the exact solution.
        Return None if no exact solution is available.
        """
        return None

    def verify(self, u, t, atol=1e-6, rtol=1e-5):
        """
        Return True if u at time points t coincides with an exact
        solution within the prescribed tolerances.
        Return None if the solution cannot be verified.
        """
        u_e = self.u_exact(t)
        if u_e is None:
            return None

        return np.allclose(u, u_e, rtol, atol)

    # subclasses may implement computation of derived
    # quantities, e.g., energy(u, t) etc

class Logistic(Problem):
    def __init__(self, a, R, U0):
        self.a = a
        self.R = R
        self.U0 = U0

    def f(self, u, t):
        a, R = self.a, self.R  # short form
        return a*u*(1 - u/R)

    def jac(self, u, t):
        a, R = self.a, self.R  # short form
        return a*(1 - u/R) + a*u*(1 - 1./R)

    def u_exact(self, t):
        a, R, U0 = self.a, self.R, self.U0  # short form
        return R*U0*np.exp(a*t)/(R + U0*(np.exp(a*t) - 1))

class LinearOscillator(Problem):
    pass
# compute potental and kinetic energy, use their sum
# to verify
# offer analytical solution with damping, also with sin/cos excitation

class VanDerPolOscillator(Problem):
    """
    Classical van der Pool oscillator:

    ..math:: y'' = \mu (1 - y^2) y' - y

    with initial conditions :math:`y(0)=2, y'(0)=1`.
    The equation is rewritten as a system

    ..math::
             u_0' &= u_1
             u_1' &= \mu (1-u_0^2)u_1 - u_0

    with a Jacobian

    ..math::
             \left(\begin{array}{cc}
             0 & 1\\
             -2\mu u_0 - 1 & \mu (1-u_0^2)
             \end{array}\right)
    """
    def __init__(self, U0=[2, 1], mu=3., f77=False):
        self.U0 = U0
        self.mu = mu

        # Compile F77
        if f77:
            self.f_f77, self.jac_f77_radau5, self.jac_f77_lsode = \
                        compile_f77([self.str_f_f77(),
                                     self.str_jac_f77_radau5(),
                                     self.str_jac_f77_lsode])

    def f(u, t):
        u_0, u_1 = u
        mu = self.mu
        return [u_1, mu*(1 - u_0**2)*u_1 - u_0]

    def jac(u, t):
        u_0, u_1 = u
        mu = self.mu
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
        """Return f(u,t) as Fortran source code string."""
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


def tester(problems, methods, time_points, terminate=None,
           compare_tol=1E-4):
    results = {}
    for pname in problems:
        for mname in methods:
            solver = methods[mname]
            problem = problems[pname]
            solver.set_initial_condition(problem.get_initial_condition())
            u, t = solver.solve(time_points, terminate)
            results[(mname, pname)] = (u, t)
    # Compare with exact solution
    # Compare solutions within tolerance
    ok = {}
    for pname in problems:
        ok[pname] = {}
        for i, mname in enumeate(methods):
            if i == 0:
                reference_solution, t = results[(mname, pname)]
                reference_name = mname
            else:
                u = results[(mname, pname)][0]
                r = np.allclose(reference_solution, u,
                                compare_tol, 0.1*compare_tol)
                ok[pname][mname] = r

import numpy as np

class Problem:
    stiff = False
    complex_ = False

    def __init__(self):
        pass

    def get_initial_condition(self):
        """Return vector of initial conditions."""
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
    def __init__(self, mu=3.):
        self.mu = mu

    def f(u, t):
        u_0, u_1 = u
        mu = self.mu
        return [u_1, mu*(1 - u_0**2)*u_1 - u_0]

    def jac(u, t):
        u_0, u_1 = u
        mu = self.mu
        return [[0., 1.],
                [-2*mu*u_0*u_1 - 1, mu*(1 - u_0**2)]]

    def f_f77(self):
        """Return f(u,t) as Fortran source code string."""
        return """
      subroutine f_f77(neq, t, u, udot)
Cf2py intent(hide) neq
Cf2py intent(out) udot
      integer neq
      double precision t, u, udot
      dimension u(neq), udot(neq)
      udot(1) = u(2)
      udot(2) = %g*(1d0 - u(1)**2)*u(2) - u(1)
      return
      end
""" % self.mu

    def jac_f77(self):
        """Return f(u,t) as Fortran source code string."""
        return """

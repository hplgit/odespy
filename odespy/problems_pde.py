from problems import Problem, np

class Diffusion1D(Problem):
    """
    Classical 1D diffusion equation:

    .. math::

          \frac{\partial u}{\partial t} = a\frac{\partial^2 u}{\partial x^2}

    with initial condition :math:`u(x,0)=I(x)` and boundary condtions
    :math:`u(0,t)=U_L(t), u(L,t)=U_R(t)`.

    .. math::
             u_0' &= u_1
             u_1' &= \mu (1-u_0^2)u_1 - u_0

    with a Jacobian

    .. math::
             \left(\begin{array}{cc}
             0 & 1\\
             -2\mu u_0 - 1 & \mu (1-u_0^2)
             \end{array}\right)
    """
    def __init__(self, I, L, U_L, U_R, n, a, x=None,
                 jac_format='banded', f77=False,
                 I_vectorized=True):
        self.I, self.U_L, self.U_R = I, U_L, U_R
        self.L, self.n, self.a = I, n, a
        if x is None:
            self.x = np.linspace(0, L, n+1)  # spatial uniform mesh
        else:
           self.x = x

        self.u   = np.zeros_like(x)
        self.u_1 = np.zeros_like(x)

        # Set initial condition
        if I_vectorized:
            self.u_1[:] = I(x)  # I(x) is vectorized
        else:
            for i in range(len(x)):
                self.u_1[i] = I(x[i])

        # Compile F77
        if f77:
            self.f_f77, self.jac_f77_radau5, self.jac_f77_lsode = \
                        compile_f77([self.str_f_f77(),
                                     self.str_jac_f77_radau5(),
                                     self.str_jac_f77_lsode])

    def f(self, u, t):
        u_0, u_1 = u
        mu = self.mu

        return [u_1, mu*(1 - u_0**2)*u_1 - u_0]

    def jac(self, u, t):
        pass

    def jac_banded(self, u, t):
        pass
    """
  - lband : None or int
 |  - rband : None or int
 |    Jacobian band width, jac[i,j] != 0 for i-lband <= j <= i+rband.
 |    Setting these requires your jac routine to return the jacobian
 |    in packed format, jac_packed[i-j+lband, j] = jac[i,j]."""

    def jac(self, u, t):
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


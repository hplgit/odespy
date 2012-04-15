__author__ = ['Liwei Wang, Univ. of Oslo']

from odespy import Solver
import numpy as np
import sys

_parameters_Radau5 = dict(

    mas = dict(
        help='''User-supplied function to
calculate the mass matrix M in the linearly
implicit ODE system M*u' = f(u, t). This function
is used only when M is not I (the identity).
M can be supplied as full matrix or banded
matrix. The ``mas`` function returns M
and takes any number of user-defined parameters
as arguments, given through the parameter
``mas_args``. If M is sparse, a banded
matrix with dimensions ``(mlmas+mumas+1,neq)``
is returned, where ``mlmas`` and ``mumas`` are
the lower and upper half-bandwidth of the banded
matrix, and ``neq`` is the number of equations.
''',
        type=callable),

    mas_args = dict(
        help='''Extra positional arguments
to the ``mas`` function (it is called as
``mas(*mas_args)``.''',
        type=(tuple, list, np.ndarray),
        default=()),

    mas_f77 = dict(
        help='''Fortran subroutine for mas.
This subroutine has the signature
(note that rpar,ipar are not used)::

        subroutine mas_f77(neq,mas,lmas,rpar,ipar)
  Cf2py intent(hide)   rpar,ipar
  Cf2py intent(in)     neq,lmas
  Cf2py intent(out)    mas
        integer neq,lmas,ipar
        double precision mas(lmas,neq),rpar
        mas(1,1) = ...
        ...
        return
        end

''',
        type=callable),

    mlmas = dict(
        help='Lower bandwith of mass matrix M',
        type=int),

    mumas = dict(
        help='Upper bandwith of mass matrix M',
        type=int),

    jac_f77_radau5 = dict(
        help='''Fortran subroutine for jac to be
provided to Radau5. This subroutine should be
defined as (note that rpar,ipar are not used)::

        subroutine jac_f77_radau5
       &           (neq,t,u,dfu,ldfu,rpar,ipar)
  Cf2py intent(hide) neq,rpar,ipar
  Cf2py intent(in)   t,u,ldfu
  Cf2py intent(out) dfu
        integer neq,ipar,ldfu
        double precision t,u,dfu,rpar
        dimension u(neq),dfu(ldfu,neq)
        dfu(1,1) = ...
        dfu(1,2) = ...
        ...
        return
        end

''',
        type=callable),
)

import solvers
solvers._parameters.update(_parameters_Radau5)


class Radau5(Solver):
    """
    This is a wrapper for radau5.f, a file containing the
    the ``radau5`` FORTRAN subtroutine for the numerical solution
    of a stiff ODEs or differential algebraic (DAE) systems,
    involving first order equations on the form

    .. math:: Mu' = f(u,t).

    The ODE system can be linearly implicit (mass-matrix :math:`M`
    is not identity matrix) or explicit.

    The method used is an implicit Runge-Kutta method (Radau IIA)
    of order 5 with step size control. The file ``radau5.f`` can
    be downloaded from `<http://www.unige.ch/~hairer/software.html>`_.
    Details about the implementation can be found in the book
    "Solving Ordinary Differential Equations II, Stiff and
    Differential-Algebraic Problems", by E. Hairer and G. Wanner,
    2nd ed., 1996, published by Springer-Verlag (ISBN: 3-540-60452-9).
    """

    _optional_parameters = Solver._optional_parameters + \
        ['f_f77', 'jac_f77_radau5', 'atol', 'jac_banded',
         'rtol', 'nsteps', 'ml', 'mu', 'first_step', 'jac', 'safety',
         'max_step', 'nsteps', 'mas', 'mas_f77', 'mlmas', 'mumas',
         'mas_args']
    _error_messages = {
        -1: 'INPUT IS NOT CONSISTENT',
        -2: 'LARGER NMAX IS NEEDED',
        -3: 'STEP SIZE BECOMES TOO SMALL',
        -4: 'MATRIX IS REPEATEDLY SINGULAR',
        }
    _extra_args_fortran = {}

    def initialize(self):
        '''Import extension module _radau5 and check that it exists.'''
        try:
            import _radau5
            self._radau5 = _radau5
        except ImportError:
            raise ImportError('Cannot find the extension module _radau5.\nRun setup.py again and investigate why _radau5.so was not successfully built.')

    def adjust_parameters(self):
        try:
            del self._parameters['jac']['default']
            del self._parameters['nsteps']['default']
        except:
            pass

    def validate_data(self):
        # lower- & upper-bound for banded Jacobian in range [0, neq]
        for name in ('ml', 'mu', 'mlmas', 'mumas'):
            if hasattr(self, name):
                self._parameters[name]['range'] = (0, self.neq+1)

        has_banded_jac = hasattr(self, 'jac_banded')
        ml, mu = getattr(self, 'ml', None), getattr(self, 'mu', None)
        if has_banded_jac and ((ml is None) or (mu is None)):
            raise ValueError('"ml" and "mu" have to be provided when banded Jacobian matrix is involved! Your input is (%s, %s).' % (ml, mu))
        return Solver.validate_data(self)

    def func_wrappers(self):
        if hasattr(self, 'f'):
            # If f is input in form of a Python function f(u,t),
            # let self.f_f77 wrap f and have arguments t, u.
            f = self.f
            self.f_f77 = lambda t,u: np.asarray(f(u,t))
        elif hasattr(self, 'f_f77'):
            # The right-hand side "f" is input as a Fortran function
            # taking the arguments t,u.
            # Set self.f to be f_f77 wrapped to the general form f(u,t)
            # for switch_to().
            f_f77 = self.f_f77
            self.f = lambda u,t: np.asarray(f_f77(t,u))

        if hasattr(self, 'jac'):
            # If jac is input as full matrix in form of jac(u,t),
            # wrap jac to jac_f77 for Fortran code.
            jac = self.jac
            self.jac_f77_radau5 = lambda t,u: np.asarray(jac(u,t),
                                                         order='Fortran')
        elif hasattr(self, 'jac_banded'):
            # If jac is input as banded matrix in form of jac_banded(u,t,ml,mu),
            # wrap jac_banded to jac_f77_radau5 for Fortran code.
            jac_banded = self.jac_banded
            self._extra_args_fortran['jac_extra_args'] = (self.ml, self.mu)
            self.jac_f77_radau5 = lambda t,u,ml,mu: \
                np.asarray(jac_banded(u,t,ml,mu), order='Fortran')
        elif hasattr(self, 'jac_f77_radau5'):
            # If jas is input as Fortran subroutine
            jac = self.jac_f77_radau5
            if 0 < getattr(self, 'ml', -1) < self.neq:
                ljac = self.ml + getattr(self, 'mu', 0) + 1
                self.jac_banded = lambda u,t,ml,mu: jac(t,u,ljac)
            else:
                self.jac = lambda u,t: jac(t,u,self.neq)

        if hasattr(self, 'mas'):
            # If mas is input as Python function in form of mas(*mas_args),
            mas = self.mas
            self.mas_f77 = lambda : np.asarray(mas(*self.mas_args),
                                               order='Fortran')
        elif hasattr(self, 'mas_f77'):
            # If mas is input as Fortran subroutine
            mas = self.mas_f77
            if 0 < getattr(self, 'mlmas', -1) < self.neq:
                lmas = self.mlmas + getattr(self, 'mumas', 0) + 1
            else:
                lmas = self.neq
            self.mas = lambda : mas(self.neq, lmas)

    def set_dummy_functions(self):
        for name in ('jac_f77_radau5', 'f_f77'):
            if getattr(self, name, None) is None:
                setattr(self, name, lambda x, y: 0.)
        if not hasattr(self, 'mas_f77'):
            setattr(self, 'mas_f77', lambda : 0.)

    def set_tol(self):
        """
        Tolerance parameters, rtol & atol, can be provided as either two scalars
        or two vectors.
        """
        rtol_scalar = isinstance(self.rtol, (int, float))
        atol_scalar = isinstance(self.atol, (int, float))
        if (rtol_scalar and (not atol_scalar)):
            self.rtol = [self.rtol] * self.neq   # Wrap to vector
        elif (atol_scalar and (not rtol_scalar)):
            self.atol = [self.atol] * self.neq   # Wrap to vector
        # itol is a flag to indicate whether tolerance parameters are input
        # as scalars or vectors
        self.itol = int(not(atol_scalar & rtol_scalar))

    def initialize_for_solve(self):
        self.func_wrappers()
        self.set_tol()
        # Flags to indicate whether jac and mas are provided
        self.ijac = int(hasattr(self, 'jac_f77_radau5'))
        self.imas = int(hasattr(self, 'mas_f77'))
        self.set_dummy_functions()

        ljac = self.neq if ((not hasattr(self, 'ml')) or \
            (self.ml == self.neq)) else (self.ml + self.mu + 1)
        lmas = 0 if not self.imas else (self.neq \
            if ((not hasattr(self, 'mlmas')) or (self.mlmas == self.neq)) \
                 else (self.mlmas + self.mumas + 1))
        le = self.neq if ((not hasattr(self, 'ml')) or \
            (self.ml == self.neq)) else (2*self.ml + self.mu + 1)
        self.lwork = self.neq*(ljac + lmas + 3*le + 12) + 20
        self.liwork = 3*self.neq + 20

        # Arrays to specify how the problem is to be solved.
        self.iwork = np.zeros(self.liwork, int)
        self.work = np.zeros(self.lwork, float)
        self.iwork[1] = getattr(self, 'nsteps', 0)
        self.work[1] = getattr(self, 'safety')
        self.work[6] = getattr(self, 'max_step', 0.)

        Solver.initialize_for_solve(self)   # Common settings


    def advance(self):
        '''
        This function intends to one step forward with Radau5.
        '''
        neq, u, n = self.neq, self.u, self.n
        t, t_next = self.t[n], self.t[n+1]

        mas = getattr(self.mas_f77, '_cpointer', self.mas_f77)
        f = getattr(self.f_f77, '_cpointer', self.f_f77)
        jac = getattr(self.jac_f77_radau5, '_cpointer', self.jac_f77_radau5)

        h = getattr(self, 'first_step', 0.)
        ml = getattr(self, 'ml', self.neq)
        mu = getattr(self, 'mu', 0)
        mlmas = getattr(self, 'mlmas', self.neq)
        mumas = getattr(self, 'mumas', 0)

        args = (f, t, u[n].copy(), t_next, h, self.rtol, self.atol,
                self.itol, jac, self.ijac, ml, mu,
                mas, self.imas, mlmas, mumas, self.work,
                self.iwork)

        # In the last demo, do not work without printing u_n out. Why?
        # Try demo_Radau5_5_fortran.py
        """
        if not hasattr(self, 'printed'):
            print u[n].copy()
            self.printed = True"""

        u_new, t_new, iwork, idid = self._radau5.advance_radau5(
                             *args, **self._extra_args_fortran)

        if idid < 0:          # Error occurred
            print self._error_messages[idid] + str(t_new)
            sys.exit(1)   # Interrupt
	return u_new

### end of class Radau5 ###

class Radau5Explicit(Radau5):
    """
    Radau5 solver for explcit ODE problem.
    """
    _optional_parameters = Solver._optional_parameters + \
        ['f_f77', 'jac_f77_radau5', 'atol', 'jac_banded',
         'rtol', 'nsteps', 'ml', 'mu', 'first_step', 'jac', 'safety',
         'max_step', 'nsteps']

Radau5Implicit = Radau5

"""Module for wrapping rkc.f."""

from solvers import Solver, Adaptive
import numpy as np

# f_f77 and other items are defined in odepack.py and will
# be populated in solvers._parameters in any import of odepack.
# We just need to add what rkc.py and that is not defined elsewhere:

_parameters_RKC = dict(

    spcrad = dict(
        help='Python function of (u, t) returning the spectral radius '\
             'of the Jacobian in the rkc.f solver.',
        type=callable),

    spcrad_f77 = dict(
        help='''Fortran version of spcrad function.
This subroutine should be defined in form:

        double precision function spcrad_f77
       1  (neq,t,u)
  Cf2py intent(hide)  neq
        integer       neq
        double precision t,u(neq)
        spcrad_f77 =
        return
        end

''',
        type=callable),
    )

import solvers
solvers._parameters.update(_parameters_RKC)

class RKC(Adaptive):
    """
    Wrapper for rkc.f, a well-known Fortran ODE solver.

    Besides the standard attributes of class ``Solver``, class ``RKC``
    also stores a dictionary ``statistics``, which contains data and
    explanations from the execution of the ``RKC`` subroutine.

    The Fortran source code can be obtained from netlib and contains
    more details. For convenience we quote here from ``rkc.f`` the
    main description of the method:

    "ABSTRACT:  RKC integrates initial value problems for systems of first
    order ordinary differential equations.  It is based on a family of
    explicit Runge-Kutta-Chebyshev formulas of order two.  The stability
    of members of the family increases quadratically in the number of
    stages m. An estimate of the spectral radius is used at each step to
    select the smallest m resulting in a stable integration. RKC is
    appropriate for the solution to modest accuracy of mildly stiff problems
    with eigenvalues of Jacobians that are close to the negative real axis.
    For such problems it has the advantages of explicit one-step methods and
    very low storage. If it should turn out that RKC is using m far beyond
    100, the problem is not mildly stiff and alternative methods should be
    considered.  Answers can be obtained cheaply anywhere in the interval
    of integration by means of a continuous extension evaluated in the
    subroutine RKCINT.

    The initial value problems arising from semi-discretization of
    diffusion-dominated parabolic partial differential equations and of
    reaction-diffusion equations, especially in two and three spatial
    variables, exemplify the problems for which RKC was designed." (rkc.f)

    This wrapper does not call ``RKCINT`` but runs ``RKC`` between each
    time interval specified by the ``time_points`` array sent to the
    ``solve`` method.
    """
    quick_description = \
        "Explicit 2nd-order Runge-Kutta-Chebyshev method (rkc.f)"

    _optional_parameters = Adaptive._optional_parameters + \
        ['f_f77', 'spcrad', 'spcrad_f77', 'jac_constant']
    # The following step parameters are illegal for rkc.f
    _optional_parameters.remove('first_step')
    _optional_parameters.remove('min_step')
    _optional_parameters.remove('max_step')

    _idid_messages = {
        3:
        'Repeated improper error control: For some j, '
        'ATOL(j) = 0 and Y(j) = 0.',
        4:
        'Unable to achieve the desired accuracy with the precision available. '
        'A severe lack of smoothness in the solution y(t) or the function '
        'f(t,y) is likely.',
        6:
        'The method used by RKC to estimate the spectral radius of the '
        'Jacobian failed to converge.',
        0:
        'Iteration stops when function TERMINATE return with True.',
        1:'Iteration succeed.'
        }

    def adjust_parameters(self):
        self._parameters['rtol']['type'] = float
        self._parameters['rtol']['range'] = (2.22e-15, 0.1)

    def check_atol(self):
        ''' ATOL need to be supplied as scalar or vector of length NEQ. '''
        atol = self.atol
        if not isinstance(atol,float):
            if len(atol) not in (1, self.neq):
                raise ValueError,  '''
ATOL =%s should be either a scalar or a vector of length NEQ=%d.
           ''' % (str(atol), self.neq)

    def validate_data(self):
        self.check_atol()
        return Adaptive.validate_data(self)

    def initialize_for_solve(self):
        # INFO(4) is an integer array to specify how the problem
        # is to be solved
        self.info = np.zeros(4, int)

        self.info[0] = 1      # Compute solution at each time point
        if hasattr(self, 'spcrad_f77') or hasattr(self, 'spcrad'):
            self.info[1] = 1  # SPCRAD routine is supplied
        else:
            self.spcrad = lambda x,y: 0.0  # dummy function
        # Is the Jacobian constant?
        self.info[2] = self.jac_constant
        if (np.iterable(self.atol) and (len(self.atol) == self.neq)):
            self.info[3] = 1   # ATOL is a sequence of length NEQ

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
        # If spcrad is input in form of spcrad(u,t),
        # wrap spcrad to spcrad_f77 for Fortran code.
        if hasattr(self, 'spcrad'):
            # If spcrad is in form of spcrad(u,t), wrap for Fortran code.
            spcrad = self.spcrad
            self.spcrad_f77 = lambda t,u: np.asarray(spcrad(u,t))

        # We call Solver and not Adaptive below because Adaptive
        # just computes first_step, min_step and max_step, all of
        # which are non-used parameters for rkc.f
        Solver.initialize_for_solve(self)   # Common settings

    def initialize(self):
        """Import extension module _rkc and check that it exists."""
        try:
            import _rkc
            self._rkc = _rkc
        except ImportError:
            raise ImportError('Cannot find the extension module _rkc.\nRemove build directory, run setup.py again and investigate why _rkc.so was not successfully built.')


    def solve(self, time_points, terminate=None):
        # flag to indicate dummy function
        itermin = int(terminate is not None)
        # Logical value cannot be transferred with f2py.
        if terminate is None:    # Dummy function
            terminate_int = lambda u, t, step_no, nt, neq: 0
        else:
            terminate_int = lambda u, t, step_no, nt, neq: \
                                   int(terminate(u, t, step_no))

        self.t = np.asarray(time_points)

        # Setting for internal parameters, (like self.u)
        self.initialize_for_solve()

        # Validity-check for values of class attributes
        if not self.validate_data():
            raise ValueError('Invalid data in "%s":\n%s' % \
                             (self.__class__.__name__,
                              pprint.pformat(self.__dict__)))

        # Call extension module
        solve = self._rkc.solve
        # Fortran or Python function
        f = getattr(self.f_f77, '_cpointer', self.f_f77)
        spcrad = getattr(self.spcrad_f77, '_cpointer', self.spcrad_f77)

        # Start iteration
        uout, idid, nstop = solve(spcrad, f, self.u[0].copy(),
                                  self.t, self.rtol, self.atol, self.info,
                                  terminate_int, itermin)

        if idid > 1:   # abnormal status?
            raise Exception('idid=%d > 1 (abort)' % idid)
        self.idid = 'idid=%d\n%s' % (idid, RKC._idid_messages[idid])

        # Store statistics from common block
        self.statistics = {
            'nfe': (self._rkc.rkcdid.nfe,
                    'number of evaluations of f used to integrate '
                    'the initial value problem'),
            'nsteps': (self._rkc.rkcdid.nsteps,
                       'no of integration steps'),
            'naccpt': (self._rkc.rkcdid.naccpt,
                       'no of accepted steps'),
            'nrejct': (self._rkc.rkcdid.nrejct,
                       'no of rejected steps'),
            'nfesig': (self._rkc.rkcdid.nfesig,
                       'no of evaluations of f used to estimate '
                       'spectral radius'),
            'maxm': (self._rkc.rkcdid.maxm,
                     'max no of stages used')}

        self.u = np.asarray(uout[:nstop,:]).copy()
        # self.u is two-dimensional, remove 2nd dimension if scalar ODE
        if self.u.shape[1] == 1:
            self.u = self.u.reshape(self.u.size)
        return self.u, self.t



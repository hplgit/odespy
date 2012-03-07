"""Module for wrapping rkf45.py."""

from solvers import Solver, Adaptive
import numpy as np


class RKF45(Adaptive):
    """
    Wrapper for rkf45.f, a FORTRAN solver designed to solve non-stiff and
    mildly stiff differential equations by the well-known
    Runge-Kutta-Fehlberg (4,5) method.

    The FORTRAN source rkf45.f can be obtained from Netlib.
    """
    quick_description = "Adaptive Runge-Kutta-Fehlberg (4,5) method (rkf45.f)"

    _optional_parameters = Adaptive._optional_parameters + ['f_f77']
    # The following step parameters are illegal for rkf45.f
    # and therefore removed (class Adaptive adds them in the above statement)
    _optional_parameters.remove('first_step')
    _optional_parameters.remove('min_step')
    _optional_parameters.remove('max_step')

    _iflag_messages = {
        2:
        'integration reached tout. indicates successful retur '
        'and is the normal mode for continuing integration.',
        -2:
        'a single successful step in the direction of tout '
        'has been taken. normal mode for continuing '
        'integration one step at a time.',
        3:
        'integration was not completed because relative error '
        'tolerance was too small. relerr has been increased '
        'appropriately for continuing.',
        4:
        'integration was not completed because more than '
        '3000 derivative evaluations were needed. this '
        'is approximately 500 steps.',
        5:
        'integration was not completed because solution '
        'vanished making a pure relative error test'
        'impossible. must use non-zero abserr to continue. '
        'using the one-step integration mode for one step '
        'is a good way to proceed.',
        6:
        'integration was not completed because requested '
        'accuracy could not be achieved using smallest '
        'allowable stepsize. user must increase the error '
        'tolerance before continued integration can be '
        'attempted.',
        7:
        'it is likely that rkf45 is inefficient for solving '
        'this problem. too much output is restricting the '
        'natural stepsize choice. use the one-step integrator '
        'mode.',
        8:
        'invalid input parameters '
        'this indicator occurs if any of the following is '
        'satisfied - t=tout  and  iflag .ne. +1 or -1 '
        'relerr or abserr .lt. 0. iflag .eq. 0  or  .lt. -2  or  .gt. 8',
        }

    def initialize(self):
        '''Import extension module _rkf45 and check that it exists.'''
        try:
            import _rkf45
            self._rkf45 = _rkf45
        except ImportError:
            raise ImportError('Cannot find the extension module _rkf45.\nRun setup.py again and investigate why _rkf45.so was not successfully built.')

    def adjust_parameters(self):
        self._parameters['rtol']['type'] = float
        self._parameters['rtol']['extra_check'] = lambda x: x >= 0.0
        self._parameters['atol']['type'] = float
        self._parameters['atol']['extra_check'] = lambda x: x >= 0.0

    def initialize_for_solve(self):
        if hasattr(self, 'f'):
            # If f is input in form of f(u,t), wrap f to f_f77 for Fortran code.
            f = self.f
            self.f_f77 = lambda t,u: np.asarray(f(u,t))
        elif hasattr(self, 'f_f77'):
            # If f is input in form of f(t,u) (usually in Fortran),
            # wrap f_f77 to the general form f(u,t) for switch_to()
            f_f77 = self.f_f77
            self.f = lambda u,t: np.asarray(f_f77(t,u))
        Solver.initialize_for_solve(self)

    def advance(self):
        u, t, n, rtol, atol = self.u, self.t, self.n, self.rtol, self.atol

        advance = self._rkf45.advance
        f = getattr(self.f_f77, '_cpointer', self.f_f77)
        unew, iflag = advance(f, u[n], t[n], t[n+1], rtol, atol)
        if iflag > 2:
            raise Exception('rkf45.f: iflag=%d > 2 (abort)' % iflag)
        self.iflag = 'iflag=%d\n%s' % (iflag, RKF45._iflag_messages[iflag])
        return unew

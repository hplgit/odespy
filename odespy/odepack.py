__author__ = ['Liwei Wang, Univ. of Oslo']

from solvers import Solver
import numpy as np
import sys, inspect

_parameters_Odepack = dict(

    # Note that ODEPACK has its own convention for the
    # arguments to the f and jac functions

    f_f77 = dict(
        help='''Fortran subroutine for f.
This subroutine has the signature::

        subroutine f_f77(neq,t,u,udot)
  Cf2py intent(hide)   neq
  Cf2py intent(out)    udot
        integer neq
        double precision t,u(neq),udot(neq)
        udot = ...
        return
        end

''',
        type=callable),

    jac_f77 = dict(
        help='''Fortran subroutine for jac.
This subroutine has the signature::

       subroutine jac_f77
      1 (neq, t, u, ml, mu, pd, nrowpd)
 Cf2py intent(hide) neq,ml,mu,nrowpd
 Cf2py intent(out) pd
       integer neq,ml,mu,nrowpd
       double precision t,u,pd
       dimension u(neq),pd(neq,neq)
       pd = ...
       return
       end

''',
        type=callable),

    jac_banded = dict(
        help='''Function for Jacobian in banded matrix form.
Used in Lsode, Lsoda, Lsodar.
``jac_banded(u,t,ml,mu)`` returns df/du as an
array of size ``neq`` times ``ml+mu+1``.''',
        paralist_old='u,t,ml,mu',
        paralist_new='t,u,ml,mu',
        array_order='Fortran',
        name_wrapped='jac_banded_f77',
        type=callable),

    jac_banded_f77 = dict(
        help='''Fortran subroutine for jac_banded.
This subroutine has the signature::

        subroutine jac_banded_f77
       1  (neq,t,u,ml, mu,pd,nrowpd)
  Cf2py intent(hide) neq,ml,mu,nrowpd
  Cf2py intent(out) pd
        integer neq,ml,mu,nrowpd
        double precision t,u,pd
        dimension u(neq),pd(nrowpd,neq)
        pd = ...
        return
        end

''',
        type=callable),

    g = dict(
        help='''Callable object to define constraint functions.
``g(u, t)`` returns a vector of the values of the constraints
(left-hand sides in the constraint equations).
 ''',
        paralist_old='u,t',
        paralist_new='t,u',
        name_wrapped='g_f77',
        type=callable),

    g_f77 = dict(
        help='''Fortran subroutine for g.
This subroutine has the signature::

        subroutine g_f77(neq, t, u, ng, groot)
  Cf2py intent(hide) neq
  Cf2py optional, intent(hide) ng
  Cf2py intent(in) t, u
  Cf2py intent(out) groot
        integer neq, ng
        double precision t, u, groot
        dimension u(neq), groot(ng)
        groot = ...
        return
        end

''',
        type=callable),

    jac_column = dict(
        help='''\
A callable object to specify a column of the Jacobian.
``jac(u,t,j,ia,ja)`` returns the j-th column of the
Jacobian.''',
        paralist_old='u,t,j',
        paralist_new='t,u,j',
        name_wrapped='jac_column_f77',
        type=callable),

    jac_column_f77 = dict(
        help='''Fortran subroutine for jac_column.
This subroutine has the signature::

        subroutine jac_column_f77
       1  (neq, t, u, j, ia, ja, pd)
  Cf2py intent(hide) neq, ia, ja
  Cf2py intent(out) pd
        integer neq, j, ia, ja
        double precision t, u, pd
        dimension u(neq), pd(neq), ia(neq + 1), ja(*)
        pd = ...
        return
        end

''',
        type=callable),

    # parameters for linearly implicit ODE solvers: Lsodi, Lsoibt, Lsodis
    res = dict(
        help='''User-supplied function to calculate the residual
vector, defined by r =  g(u,t) - A(u,t) * s.
Used in the linearly implicit solvers: Lsodi,
Lsodis, Lsoibt. The ``res`` function has the
signature ``res(u,t,s,ires)`` and returns
the tuple ``(r,ires)``, where ``ires`` is an
int. On input, ires indicates how ODEPACK would use
the returned array "r": ``ires=1`` means the full
residual exactly, ``ires=2`` means that ``r`` is
used only to compute the Jacobian dr/du by finite
differences.
``res`` should set the flag ``ires`` if it encounters
a halt condition or illegal input. Otherwise, it
should not be reset. On output, the value 1 or -1
represents a normal return.
''',
        paralist_old='u,t,s,ires',
        paralist_new='t,u,s,ires',
        name_wrapped='res_f77',
        type=callable),

    res_f77 = dict(
        help='''Fortran subroutine for res.
This subroutine has the signature::

      subroutine res_f77(neq, t, u, s, r, ires)
 Cf2py intent(hide) neq
 Cf2py intent(out) r
 Cf2py intent(in,out) ires
      double precision t, u, s, r
      dimension u(neq, s(neq), r(neq)
      ...
      return
      end

''',
        type=callable),

    jac_lsodi = dict(
        help='''Callable object to define the full Jacobian
matrix dr/du where r = g - A*s. The
signature of this function is ``jac(u,t,s)``,
returning a matrix like the other ``jac``
functions.
''',
        paralist_old='u,t,s',
        paralist_new='t,u,s',
        name_wrapped='jac_lsodi_f77',
        type=callable),

    jac_lsodi_f77 = dict(
        help='''Fortran subroutine for jac_lsodi.
This subroutine has the signature::

       subroutine jac_lsodi_f77
      1  (neq, t, u, s, ml, mu, pd, nrowpd)
 Cf2py intent(in, hide) neq, ml, mu, nrowpd
 Cf2py intent(out) pd
       integer neq, ml, mu, nrowpd
       double precision t, u, pd, s
       dimension u(neq), s(neq), pd(nrowpd, neq)
       pd = ...
       return
       end

''',
        type=callable),

    jac_banded_lsodi = dict(
        help='''Callable object to define the banded Jacobian
matrix dr/du where r = g - A*s. The
signature is ``jac(u,t,s,ml,mu)``,
where ``ml`` and ``mu`` are the lower and
upper half-bandwidth of the banded matrix.
''',
        paralist_old='u,t,s,ml,mu',
        paralist_new='t,u,s,ml,mu',
        name_wrapped='jac_banded_lsodi_f77',
        type=callable),

    jac_banded_lsodi_f77 = dict(
        help='''Fortran subroutine for jac_banded_lsodi.
This subroutine has the signature::

       subroutine jac_banded_lsodi_f77
      1  (neq, t, u, s, ml, mu, pd, nrowpd)
 Cf2py intent(in, hide) neq, ml, mu, nrowpd
 Cf2py intent(out) pd
       integer neq, ml, mu, nrowpd
       double precision t, u, pd, s
       dimension u(neq), s(neq), pd(nrowpd, neq)
       pd = ...
       return
       end

''',
        type=callable),

    adda_lsodi = dict(
        help='''Callable object to add the matrix A = A(u,t)
to another matrix p stored in the same form as A.
The signature is ``adda(u,t,p)`` and it returns
a matrix p+A (square matrix as the Jacobian).
''',
        paralist_old='u,t,p',
        paralist_new='t,u,p',
        array_order='Fortran',
        name_wrapped='adda_lsodi_f77',
        type=callable),

    adda_lsodi_f77 = dict(
        help='''Fortran subroutine for adda_lsodi.
This subroutine has the signature::

        subroutine adda_lsodi_f77
       1  (neq, t, u, ml, mu, pd, nrowpd)
  Cf2py intent(in, hide) neq, ml, mu
  Cf2py intent(in, hide), depend(pd) nrowpd
  Cf2py intent(in, out) pd
        integer neq, ml, mu, nrowpd
        double precision t, u, pd
        dimension u(neq), pd(nrowpd, neq)
        pd = ...
        return
        end

''',
        type=callable),

    adda_banded_lsodi = dict(
        help='''Callable object to add the banded matrix
A = A(u,t) to another matrix stored P in the
same form as A. For a banded matrix, A(i,j) is
added to P(i-j+mu+1,j). The signature is
``adda(u,t,p,ml,mu)`` and it returns a banded
matrix.
''',
        paralist_old='u,t,p,ml,mu',
        paralist_new='t,u,p,ml,mu',
        name_wrapped='adda_banded_lsodi_f77',
        type=callable),

    adda_banded_lsodi_f77 = dict(
        help='''Fortran subroutine for adda_banded.
This subroutine has the signature::

       subroutine adda_banded_lsodi_f77(neq, t,
      1                  u, ml, mu, pd, nrowpd)
 Cf2py intent(in, hide) neq, ml, mu
 Cf2py intent(in, hide), depend(pd) nrowpd
 Cf2py intent(in, out) pd
       integer neq, ml, mu, nrowpd
       double precision t, u, pd
       dimension u(neq), pd(nrowpd, neq)
       pd = ...
       return
       end

''',
        type=callable),

    jac_lsodis = dict(
        help='''Callable object to supply the j-th column of
the sparse Jacobian matrix dr/du where
r = g - A*s. The signature is ``jac(u,t,s,j,ia,ja)``
and a vector is returned.
''',
        paralist_old='u,t,s,j-1,ia,ja',
        paralist_new='t,u,s,j,ia,ja',
        name_wrapped='jac_lsodis_f77',
        type=callable),

    adda_lsodis = dict(
        help='''Callable object to add j-th column of matrix
A = A(u,t) to another matrix stored in sparse
form. The signature is ``adda(u,t,j,ia,ja,p)``
and it returns a vector.
''',
        paralist_old='u,t,j-1,ia,ja,p',
        paralist_new='t,u,j,ia,ja,p',
        name_wrapped='adda_lsodis_f77',
        type=callable),
    jac_lsoibt = dict(
        help='''Callable object to supply the jth column of
the Jacobian matrix dr/du where r = g - A*s,
stored in block-tridiagonal form.
The signature is jac(u,t,s), and the return value
is a tuple (pa,pb,pc), where each of these
arrays has size mb*mb*nb, mb and nb being two
parameters that can be set (see the help line for
the parameters in the doc string of Lsoibt).
pa, pb, and pc are to be loaded with partial
derivatives (elements of the Jacobian matrix)
on output, in terms of the block-tridiagonal
structure assumed. That is, load the diagonal
blocks into pa, the superdiagonal blocks (and
block (nb,nb-2)) into pb, and the subdiagonal
blocks (and block (1,3)) into pc.
The blocks in block-row k of dr/du are to be
loaded into pa(*,*,k), pb(*,*,k), and
pc(*,*,k).
Thus the affect of this function should be::

  pa(i,j,k) = ( (i.j) element of k-th diagonal
              block of dr/du)
  pb(i,j,k) = ( (i,j) element of block (k,k+1)
              of dr/du, or block (nb,nb-2) if
              k == nb)               \
  pc(i,j,k) = ( (i,j) element of block (k,k-1)
              of dr/du, or block (1,3) if k==1).''',
        paralist_old='u,t,s',
        paralist_new='t,u,s',
        name_wrapped='jac_lsoibt_f77',
        type=callable),

    adda_lsoibt = dict(
        help='''Callable object to add matrix A = A(u,t) to
another matrix P, stored in block-tridiagonal
form. The signature is  adda(u,t,pa,pb,pc), which
should return (pa,pb,pc) as described for the
jac_lsoibt parameter.
''',
        paralist_old='u,t,pa,pb,pc',
        paralist_new='t,u,pa,pb,pc',
        name_wrapped='adda_lsoibt_f77',
        type=callable),

    seth = dict(
        help='Element threshhold for sparsity determination.',
        default=0,
        type=int),

    corrector_iter_method = dict(
        help='Corrector iteration method choice.',
        type=int),

    lrw = dict(
        help='Length of real work array.',
        type=int),

    liw = dict(
        help='Length of integer work array, similiar as <lrw>.',
        type=int),

    moss = dict(
        help=' Method to obtain sparse structure of Jacobian.',
        type=int),

    max_hnil = dict(
        help='Maximum no of warning messages to be printed.',
        type=int),

    max_ordn = dict(
        help='Maximum order in nonstiff methods.',
        type=int),

    max_ords = dict(
        help='Maximum order in stiff methods.',
        type=int),

    # ja, ia, ja, ic are used to describe sparse structure of
    # matrices in Lsodes and Lsodis.
    ja = dict(
        help='''Integer array containing the row indices where
the nonzero elements occur, in columnwise
order. Describes the sparsity matrix structure
together with ia.
In Lsodes, ia and ja describe the structure of
Jacobian matrix; while in Lsodis, ia and ja are
used to describe the structure of matrix A.''',
        type=(list, tuple, np.ndarray),
        # integer sequence
        extra_check=lambda int_seq: \
            np.array(map(lambda x: isinstance(x, int),int_seq)).all()),

    ia = dict(
        help='''Integer array with length neq+1 which contains
starting locations in ja of the descriptions
for columns 1...neq. ia(1) == 1. The last element
ia[neq+1] should equal to the total number of
nonzero locations assumed.
For each column j = 1...neq, the values of the
row index i in column j, where a nonzero element
may occur, are given by i == ja(k) where ia(j) <='
k < ia(j+1).''',
        type=(list, tuple, np.ndarray),
        # integer sequence
        extra_check=lambda int_seq: \
             np.array(map(lambda x: isinstance(x, int),int_seq)).all()),

    jc = dict(
        help='''Integer array which describes the sparsity
Jacobian structure together with ic, like ia
and ja. In Lsodis, ia and ja describe the sparse
structure of matrix A, while ic and jc describe
the sparse structure of Jacobian matrix.''',
        type=(list, tuple, np.ndarray),
        # integer sequence
        extra_check=lambda int_seq: \
            np.array(map(lambda x: isinstance(x, int),int_seq)).all()),

    ic = dict(
        help='Array which contains starting locations in jc.',
        type=(list, tuple, np.ndarray),
        # integer sequence
        extra_check=lambda int_seq: \
            np.array(map(lambda x: isinstance(x, int),int_seq)).all()),

    # mb, nb describe the block-tridiagonal form of matrix.
    mb = dict(
        help='Block size. Describe the block-tridiagonal form '\
              'of matrix A together with nb.',
        type=int,
        extra_check=lambda x: x>=1),

    nb = dict(
        help='''Number of blocks in the main diagonal.
In each of the nb block-rows of the matrix P
(each consisting of mb consecutive rows), the
nonzero elements are to lie in three
consecutive mb by mb blocks.  In block-rows 2
through nb-1, these are centered about the
main diagonal. In block rows 1 and nb, they
are the diagonal blocks and the two blocks
adjacent to the diagonal block.  (Thus block
positions (1,3) and (nb,nb-2) can be nonzero.)
Require: mb>=1, nb>=4, mb*nb==neq.''',
        type=int,
        extra_check=lambda x:x>=4),

    )

import solvers
solvers._parameters.update(_parameters_Odepack)


class Odepack(Solver):
    """
    This is a superclass for wrapping seven solvers in the Fortran
    package ODEPACK, available in the Netlib repository
    `<http://www.netlib.org/odepack>`_.

    For a quick intro to ODEPACK and its relation to other packages.
    see `<https://computation.llnl.gov/casc/odepack/odepack_home.html>`_.

    *Solvers for explicitly given systems.*
    For each of the following solvers, it is assumed that the ODEs are
    given explicitly, so that the system can be written in the form
    du/dt = f(u,t), where u is a vector of dependent variables, and t
    is a scalar.

    ===============  ==========================================================
    Name             Description
    ===============  ==========================================================
       Lsode         A wrapper of dlsode, the basic solver in ODEPACK for
                     stiff and nonstiff systems of the form u' = f.

                     In the stiff case, it treats the Jacobian matrix df/du as
                     either a dense (full) or a banded matrix, and as either
                     user-supplied or internally approximated by differences.

                     It uses Adams methods (predictor-corrector) in the
                     nonstiff case, and Backward Differentiation Formula (BDF)
                     methods (the Gear methods) in the stiff case.  The linear
                     systems that arise are solved by direct methods (LU
                     factorization/backsolve).

       Lsodes        Solves systems u' = f, and in the stiff case
                     treats Jacobian matrix in general sparse form. It can
                     determine the sparsity structure on its own, or optionally
                     accepts this information from the user.
                     It then uses parts of the Yale Sparse Matrix Package (YSMP)
                     to solve the linear systems that arise, by a sparse
                     (direct) LU factorization/backsolve method.

       Lsoda         Solves systems u' = f, with a dense or banded
                     Jacobian when the problem is stiff, but it automatically
                     selects between nonstiff (Adams) and stiff (BDF)
                     methods. It uses the nonstiff method initially, and
                     dynamically monitors data in order to decide which
                     method to use.

       Lsodar        A variant of Lsoda allowing for constraints.
                     It solves u' = with dense or banded
                     Jacobian and automatic method selection, and at the same
                     time, it solves g(u,t) = 0. This is often useful for
                     finding stopping conditions, or for finding points
                     at which a switch is to be made in the function f.
    ===============  ==========================================================

     *Solvers for linearly implicit systems.*
     The following solvers treat systems in the linearly implicit form::

            A(u,t) du/dt = g(u,t),

     where A is a square matrix, i.e., with the derivative
     du/dt implicit, but linearly so.
     These solvers allow A to be singular, in which case the system is a
     differential-algebraic equation (DAE) system.
     In that case, the user must be very careful to supply a well-posed
     problem with consistent initial conditions.

    ===============  ==========================================================
    Name             Description
    ===============  ==========================================================
       Lsodi         Solves linearly implicit systems in which the
                     matrices involved (A, dg/du, and d(A u')/du) are all
                     assumed to be either dense or banded.

       Lsodibt       Solves linearly implicit systems in which the matrices
                     involved are all assumed to be block-tridiagonal.  Linear
                     systems are solved by the LU method.

       Lsodis        Solves linearly implicit systems in which the
                     matrices involved are all assumed to be sparse.
                     Either determines the sparsity structure or accepts it from
                     the user, and uses parts of the Yale Sparse Matrix Package
                     to solve the linear systems that arise, by a direct method.
    ===============  ==========================================================

    *Note*: For large ODE systems the user is encouraged that users provide
    an f2py-compiled Fortran subroutine or a multi-line string Fortran code
    to define the ODE. This would help to improve efficiency.
    """

    # Default parameter-list for all solvers in OdePack
    # Note: f and jac are not valid in some specific solvers.

    _optional_parameters = Solver._optional_parameters + \
        ['atol', 'rtol', 'adams_or_bdf', 'order',
         'nsteps', 'first_step', 'min_step',
         'max_step', 'corrector_iter_method', 'lrw', 'liw']

    # Error messages to print out,
    # corresponding to different values of return status-flag.
    _error_messages = {\
    -1: '''\
    An excessive amount of work (more than "max_step" steps) was done on this
    call, before completing the requested task, but the integration was
    otherwise successful as far as ''',
    -2: '''\
    Too much accuracy was requested for the precision of the machine being
    used.  To continue, the tolerance parameters must be reset. This was
    detected before completing the requested task, but the integration was
    successful as far as ''',
    -3: '''\
    Illegal input was detected, before taking any integration steps.
    Current T is ''',
    -4: '''\
    There were repeated error-test failures on one attempted step, before
    completing the requested task. The problem may have a singularity, or the
    input may be inappropriate, but integration was successful as far as ''',
    -5: '''\
    There were repeated convergence-test failures on one attempted step, before
    completing the requested task. This may be caused by an inaccurate Jacobian
    matrix, if one is being used. But integration was successful as far as''',
    -6: '''\
    EWT(i) became zero for some i during the integration. Pure relative error
    control (ATOL(i)=0.0) was requested on a variable which has now vanished.
    The integration was successful as far as ''',
    -7:''}

    _extra_args_fortran = {}
    # Dictionary for extra parameters in all external functions in Fortran
    # package ODEPACK.
    # For example, _extra_args_fortran['jac_extra_args'] is defined as
    # ('self.ml','self.mu') when jac is a banded Jacobian matrix.

    _iwork_index = {4:'order', 5:'nsteps', 6:'max_hnil'}
    # Index in iwork_in to supply optional inputs.
    _rwork_index = {4:'first_step', 5:'max_step', 6:'min_step'}
    # Index in rwork_in to supply optional inputs.

    def adjust_parameters(self):
        """ Special settings for properties of input parameters."""
        # If f is input in form of f(u,t), wrap f to f_f77 for Fortran code.
        if 'f' in self._parameters:
            self._parameters['f']['paralist_old'] = 'u,t'
            self._parameters['f']['paralist_new'] = 't,u'
            self._parameters['f']['name_wrapped'] = 'f_f77'
        # If f is input in form of f(t,u),
        # wrap f_f77 to the general form f(u,t) for switch_to()
        if 'f_f77' in self._parameters:
            self._parameters['f_f77']['paralist_old'] = 't,u'
            self._parameters['f_f77']['paralist_new'] = 'u,t'
            self._parameters['f_f77']['name_wrapped'] = 'f'

        # These default values are unnecessary in ODEPACK.
        for name in ('order', 'jac', 'nsteps', 'seth'):
            if name in self._parameters:
                if 'default' in self._parameters[name]:
                    del self._parameters[name]['default']

        # If jac is input in form of jac(u,t),
        # wrap jac to jac_f77 for Fortran code.
        if 'jac' in self._parameters:
            self._parameters['jac']['paralist_old'] = 'u,t'
            self._parameters['jac']['paralist_new'] = 't,u'
            self._parameters['jac']['array_order'] = 'Fortran'
            self._parameters['jac']['name_wrapped'] = 'jac_f77'
        # If jac is input in form of jac(t,u),
        # wrap jac_f77 to the general form jac(u,t) for switch_to().
        if 'jac_f77' in self._parameters:
            self._parameters['jac_f77']['paralist_old'] = 't,u'
            self._parameters['jac_f77']['paralist_new'] = 'u,t'
            self._parameters['jac_f77']['name_wrapped'] = 'jac'
        Solver.adjust_parameters(self)
        return None

    def initialize(self):
        '''Import extension module _odesolver and check that it exists.'''
        try:
            import _odepack
            self._odepack = _odepack
        except ImportError:
            raise ImportError('Cannot find the extension module _odepack.\nRun setup.py again and investigate why _odepack.so was not successfully built.')


    def func_wrapper(self):
        '''
        This function is defined to wrap user-defined functions with new
        forms of parameter-list, or wrap the returned values as numpy arrays.

        Firstly, in odespy, all the user-supplied functions should have a
        parameter list starts with "u,t,...". But in some special subclasses,
        (like solvers in ODEPACK), all the parameter lists of user-defined
        functions start with "t,u,...". So we need this general function to
        wrap all these user-defined functions.

        Secondly, in some user-defined functions, according to the different
        start indices in Fortran and Python, we need to make special wrapping
        for these uncompability. For an example, in user-defined function
        "jac_column", column index is an internally valued parameter in
        Fortran code. In Python, it starts from 0 instead of 1 in Fortran.
        So we need to wrap the parameter list of user-defined "jac_column" from
        "u,t,j" to "t,u,j+1". That is, define the Jacobian function as
        lambda t,u,j: jac_column(u,t,j-1).

        Furthermore, the return value of user-defined functions need to be
        wrapped to Numpy arrays with great numerical features, e.g.
        vectorization and array slicing. In order to avoid unnecessary array
        copy by F2PY, it is always recommended to explicitly transform all
        Numpy arrays to Fortran order in Python code.

        This functions is not intended for simple solvers. So it is not called
        automatically in current version. But for some complicated solvers as
        ones in ODEPACK, it is very useful and convenient.

        Future developers can call this functions with appropriate locations
        and corresponding property-setting in adjust_parameters().

        '''
        import numpy as np
        parameters = self._parameters
        # Extract function parameters that are required to be wrapped
        func_list = [[name,
                      parameters[name].get('array_order', None),
                      parameters[name].get('paralist_old', None),
                      parameters[name].get('paralist_new', None),
                      parameters[name].get('name_wrapped', name)]
                     for name in parameters \
                         if name in self.__dict__ and \
                         'type' in parameters[name] and \
                         (parameters[name]['type'] is callable or \
                          parameters[name]['type'] is (callable, str)) and \
                         ('paralist_new' in parameters[name] or \
                          'array_order' in parameters[name])]
        # name in self.__dict__  --> existing attributes in current instance
        # parameters[name]['type'] is callable or (callable, str)
        #             -->     callable objects
        # 'paralist_new' in parameters[name]
        #    --> new parameter-list is defined to be wrapped
        # 'array_order' in parameters[name]
        #    --> this function return an array, and should be wrapped either in
        #    Fortran order or C (default) order.
        func_input = {}
        for name, order, arg_old, arg_new, name_new in func_list:
            # e.g. name     = 'jac'
            #      arg_old  = 'u,t'
            #      arg_new  = 't,u'
            #      order    = 'Fortran'
            #      name_new = 'jac_f77'
            #  Then:
            #  self.jac_f77 = lambda t,u: np.asarray(jac(u,t), order='Fortran')
            func_input[name] = getattr(self, name)
            wrap_string = 'lambda %s: ' % \
                (arg_new if arg_new is not None else arg_old)
            wrap_string += 'np.asarray(' if order is not None else ''
            wrap_string += 'func_input["%s"](%s)' % (name, arg_old)
            wrap_string += ', order="Fortran"' if order=='Fortran' else ''
            wrap_string += ')' if order is not None else ''
            setattr(self, name_new, eval(wrap_string, locals()))

    def check_liwlrw(self):
        '''
        If the lengths of work arrays are specified by users, check whether
        they are greater than the required lengths of Fortran solvers. '''
        for name in ('liw', 'lrw'):
            min_value = getattr(self, name+'_min')
            if not hasattr(self, name):
                setattr(self, name, min_value)
            else:
                value = getattr(self, name)
                if value < min_value:
                    print '''
          Insufficient input! "%s"=%d are reset to be the minimum size = %d '''\
                        % (name, value, min_value)
                    setattr(self, name, min_value)

    def check_tol(self):
        '''
        atol & rtol should be defined as scalars or vectors with length neq.
        '''
        for name in ('atol', 'rtol'):
            value = getattr(self, name)
            if np.iterable(value):
                if len(value)== 1:
                    setattr(self, name, value[0])
                elif len(value)!= self.neq:
                    raise AttributeError(
                    '%s has an illegal length = %d! Valid length is %d or 1.' \
                        % (name, len(value), self.neq))

    def check_pars(self):
        ''' Some pairs of parameters should be input simutaneously.'''
        for pars in (('ia', 'ja'), ('ic', 'jc'), ('ml', 'mu'), ('mb', 'nb')):
            arg_a, arg_b = pars
            if int(hasattr(self, arg_a) + hasattr(self, arg_b)) == 1:
                raise ValueError,'\
        Error! %s and %s have to be input simutaneously!' % (arg_a, arg_b)

    def check_iaja(self):
        '''
        ``ia``, ``ja``, ``ic``, ``jc`` are optional inputs to describe
        arbitrary sparse structure of matrix.

        ``ia`` and ``ja`` are used in dlsodes, dlsodis.
        ``ic`` and ``jc`` are used only in dlsodis.

        There are special requirements for their values::

            len(ia/ic) = neq + 1
            (ia/ic)[0] = 1
            (ia/ic)[-1] = 1 + len(ja/jc)

        '''
        for pars in (('ia','ja'),('ic','jc')):
            arg_a, arg_b = pars
            if hasattr(self, arg_a):
                array_a,array_b = getattr(self, arg_a), getattr(self, arg_b)
                err_messages = ['''\
        Illegal input! %s(=%s) should have a length %d (= neq+1)!\n\n'''\
                   % (arg_a, str(array_a), self.neq+1),
                                '''\
        Illegal input! %s[0](=%d) should be equal to 1!\n\n'''\
                   % (arg_a, array_a[0]),
                                '''\
        Illegal input! %s[neq+1](=%d) should be equal to 1
        plus the length of %s(=%d), which represents the total
        number of nonzero locations assumed in the matrix.\n''' \
                   % (arg_a, array_a[-1], arg_b, len(array_b)+1)]
                iaja_check = (len(array_a) != self.neq+1,
                              array_a[0] != 1,
                              array_a[-1] != len(array_b)+1)
                for error_index in range(3):
                    if iaja_check[error_index]:
                        raise ValueError, err_messages[error_index]

    def validate_data(self):
        '''
        Common validity check in Odepack.
        '''
        # lower- & upper-bound for banded Jacobian in range [0,neq]
        for name in ('ml', 'mu'):
            if hasattr(self, name):
                self._parameters[name]['range'] = (0, self.neq+1)

        if not Solver.validate_data(self):
            return False
        self.check_tol()
        self.check_pars()
        self.check_iaja()
        self.check_liwlrw()
        return True

    def set_iwork_rwork(self):
        '''
        Initialize arrays for optional inputs, and calculate the
        required length of work arrays in Fortran code.
        '''

        # calculate the required length of work arrays
        self.set_liw_min()
        self.set_lrw_min()

        # Indices of optional inputs in work arrays
        iwork_index, rwork_index = self._iwork_index, self._rwork_index

        # Initialize work arrays.
        length_iwork_in = max(iwork_index.keys()) + 1
        length_rwork_in = max(rwork_index.keys()) + 1
        self.iwork_in = np.zeros(length_iwork_in, int)
        self.rwork_in = np.zeros(length_rwork_in, float)

        # Put optional inputs into work arrays with specified indices.
        for index in iwork_index:
            self.iwork_in[index] = getattr(self, iwork_index[index], 0)
        for index in rwork_index:
            self.rwork_in[index] = getattr(self, rwork_index[index], 0.)

    def set_iopt(self):
        '''
        Initialization for "ipot", which is a flag to indicate whether optional
        parameters are specified in work arrays.
        '''
        raise NotImplementedError

    def set_ydoti(self):
        '''
        ``ydoti`` is an array used in linearly implicit solvers.
        It has to be extended if its length is smaller than neq.
        '''
        ydoti = getattr(self, 'ydoti', [])
        # flag to indicate whether ydoti is supplied
        self.ydoti_flag = int(ydoti!=[])
        # extend the length of 'ydoti' to neq
        self.ydoti = list(ydoti) + [0.]*(self.neq - len(ydoti))

    def set_liw_min(self):
        '''
        Calculate the necessary length of integer work arrays when it is not
        specified explicitly by users.
        Different solvers have different formulas.
        '''
        raise NotImplementedError

    def set_lrw_min(self):
        '''
        Calculate the necessary length of real work arrays for Fortran code.
        Different solvers have different formulas.
        '''
        raise NotImplementedError

    def set_extra_args(self):
        '''Setting for extra parameters of user-defined functions.'''
        return None

    def set_jac(self):
        '''
        Set values for Jacobian matrix. In several solvers like Lsode,
        Jacobian matrix could be supplied either in full form or in banded form.

        This function intends to tell from which kind of Jacobian matrix
        is specified.
        '''
        return None

    def set_dummy_functions(self):
        '''
        Functions have to get dummy values before they are passed to extension
        module even if they are not involved in current solver.
        '''
        for name in ('jac_f77', 'f_f77'):
            if getattr(self, name, None) is None:
                setattr(self, name, lambda x,y:0.)
        if getattr(self, 'jac_column_f77', None) is None:
            self.jac_column_f77 = lambda x,y,z: 0.
        if getattr(self, 'g_f77', None) is None:
            self.g_f77 = lambda x,y: np.array(())
        if getattr(self, 'ng', None) is None:
            self.ng = 0


    def set_iter_method(self):
        '''
        Set proper values for method-choices when it is not specified
        explicitly.
        '''
        raise NotImplementedError

    def initialize_for_solve(self):
        '''
        In the long parameter-lists for solvers in ODEPACK, quite a few
        parameters can be set automatically with proper values depending
        on values of other parameters.
        In this function, all the parameters of this kind are initialized
        with proper values.
        '''
        self.set_iter_method()     # method choice
        self.func_wrapper()        # wrap function parameters
        self.set_jac()             # tell from which kind of Jacobian matrix
                                   # are supplied:  banded, full or None
        self.set_dummy_functions() # pass dummy values to extension module
                                   # for functions that are not involved
        self.set_extra_args()      # Extra parameters for external functions
        self.mf = self.iter_method + (1 + (self.adams_or_bdf == 'bdf'))*10 + \
            getattr(self, 'moss', 0)*100
        self.set_iwork_rwork()     # work arrays initialization
        self.set_iopt()            # Flag for optional inputs
        self.set_ydoti()           # Extend length of array ydoti
        # itol is a flag to indicate whether tolerance parameters are input
        # as scalars or vectors
        self.itol = (not isinstance(self.rtol, (int, float)))*2 + \
            (not isinstance(self.atol, (int, float))) + 1
        Solver.initialize_for_solve(self)   # Common settings in super class

    def new_stepnr(self):
        '''
        When Fortran code returns a status ``istate==-1``, it indicates that
        there are excessive amount of steps detected.
        Then we could try to increase ``nsteps`` to avoid this error.
        '''
        nsteps = getattr(self, 'nsteps', 500)
        if nsteps == 2000:    # The maximum step-amount has been reached
            raise ValueError, '''
        Failed iteration although step number has been set to 2000.
        Please check your input.'''
        mx_new = min(nsteps+200, 2000)
        # maximum limitation is set to 2000
        self.nsteps = mx_new  # valid for the following steps
        print '''\
        Excessive amount of work detected!
        Input step amount "nsteps"=%d is not enough for iteration!
        nsteps has been reset to %d to avoid this error!'''\
        % (nsteps, mx_new)
        return mx_new

    def tol_multiply(self, tolsf):
        '''
        This function is used to adjust tolerance parameters for Fortran part.
        When extension module returns a status "istate" as -2 or -3, it often
        indicates that there are excessive amount of steps detected.
        Then we could try to adjust tolerance settings with suggested factor
        to avoid this error.
        '''
        print 'Tolerance is scaled by suggested factor %.2g''' % tolsf
        self.rtol *= tolsf
        self.atol *= tolsf

    def expand_rwork(self, new_lrw, expand=False):
        '''
        Length of real work array is smaller than actually required length.
	Then we could expand work array to avoid this error.
        '''
        print 'The length of real work array has been reset to %d' % new_lrw
        if expand:          # Expand real arrays for linearly implicit solvers
            self.rwork = list(self.rwork) + [0.]*(new_lrw-self.lrw)
        self.lrw = new_lrw

    def expand_iwork(self, new_liw, expand=False):
        '''
        Extension module return an actually required length for
        integer work array when it is too short.
	Then we could expand work array to required length to avoid this error.
        '''
        print 'The length of integer work array has been reset to %d' % new_liw
        if expand:          # Expand integer arrays for linearly implicit solvers
            self.iwork = list(self.iwork) + [0.]*(new_liw - self.liw)
        self.liw = new_liw

    def adjust_atol(self, u_current):
        '''
        Error tolerance ``tol(i)`` may become zero for some ``i``
        during integration, where::

          tol = rtol(i) * u(i) + atol(i)

        It indicates that pure absolute tolerance (``atol(i)=0.0``)
        was requested.  In order to avoid possible divide-by-zero
        error, we find the indices of zero items and adjust
        ``atol`` to ``1e-8``.
        '''
        tol = abs(np.asarray(u_current))*self.rtol + self.atol
        if isinstance(self.atol, float):
            # atol is set as a scalar float
            # Convert scalar "atol" to be an array
            self.atol = self.atol + np.zeros(self.neq, float)
            self.itol += 1
            for index, item in enumerate(tol):
                if item == 0. == self.atol(index):
                    # Increase absolute error
                    self.atol[index] += 1e-8

    def print_roots(self, jroot, t_current, u_current):
        '''Roots found at current T for some constraint functions. '''
        g, ng = self.g_f77, self.ng
        if hasattr(g, '_cpointer'):  # ng is required if g is in Fortran
            value = g(t_current, u_current, ng=ng)
        else:
            value = g(t_current, u_current)
        for i in range(ng):
            if jroot[i]:   # found root for i-th constraint equation
                print '''
        Root found at t = %g for %dth constraint function in g''' \
                    % (t_current, i+1)
                value_ith = value[i] if ng > 1 else value
                print 'Error in location of root is %g' % value_ith


    def solve(self, time_points, terminate=None):
        '''
        This function is involved for non-linearly implicit
        solvers in ODEPACK, i.e., Lsode, Lsoda, Lsodar, and Lsodes.
        '''

        itermin = int(terminate is not None)   # flag to indicate dummy function
        # Logical value cannot be transferred with f2py.
        if terminate is None:    # Dummy function
            terminate_int = lambda u,t,step_no,nt,neq: 0
        else:
            terminate_int = lambda u,t,step_no,nt,neq: \
                int(terminate(u, t, step_no))

        self.t = np.asarray(time_points)
        self.initialize_for_solve()
        if not self.validate_data():
            raise ValueError('Invalid data in "%s":\n%s' % \
                (self.__class__.__name__,pprint.pformat(self.__dict__)))

        # Convert class-name to name of subroutine,
        # e.g. Lsode -> dlsode or slsode
        solver_name = 'd' + self.__class__.__name__.lower()

        step_no = len(self.t)
        nstop, istate, self.finished = 1, 1, False
        u = np.asarray(self.u.copy(), order="Fortran")
        tried = 0        # Number of attempts

        while not self.finished:
            istate = 1
            # Extract _cpointers if functions are compiled with F2PY.
            f = getattr(self.f_f77, '_cpointer', self.f_f77)
            g = getattr(self.g_f77, '_cpointer', self.g_f77)
            jac = getattr(self.jac_f77, '_cpointer', self.jac_f77)
            jac_column = getattr(self.jac_column_f77, '_cpointer',
                                 self.jac_column_f77)
            # call extension module
            nstop, u, istate, rinfo, iinfo = \
                apply(self._odepack.solve,
                      (terminate_int, itermin, nstop, f, u,
                       self.t, self.itol, self.rtol, self.atol,
                       istate, self.iopt, self.rwork_in, self.lrw,
                       self.iwork_in, self.liw, jac, jac_column,
                       self.mf, g, self.ng, solver_name),
                      self._extra_args_fortran)
            tried += 1
            if nstop == step_no:    # successful
                self.finished = True
                tried = 0
            if istate == 3:         # roots founded for constraint equations
                self.print_roots(iinfo[17:], self.t[nstop-1], u[nstop-1])
                tried = 0
            elif istate == 2:
                self.finished = True
                tried = 0
            elif istate == 0:
                print "Iteration stops at step Nr.%d," % nstop
                print " when function TERMINATE return with True."
                self.finished = True
            elif istate < 0:                   # Error occurs!
                print self._error_messages[istate] + str(rinfo[1])
                if istate == -1:    # Increase maximum step-number.
                    self.iwork_in[5] = self.new_stepnr()
                    self.iopt = 1
                elif istate == -2:  # Multiply tolerance with suggested factor
                    self.tol_multiply(rinfo[3])
                elif istate == -3:
                    # Illegal input was detected,
                    # before taking any integration steps.
                    if iinfo[7] > self.lrw:   # Real work array is too short.
                        self.expand_rwork(iinfo[7])
                    elif iinfo[8] > self.liw: # Integer work array is too short.
                        self.expand_iwork(iinfo[8])
                    elif rinfo[3] > 0.:       # suggested factor for tolerance
                        self.tol_multiply(rinfo[3])
                    else:  # Other cases with istate returned as -3
                        sys.exit(1)   # Interrupt
                elif istate == -6:
                    # divide-zero-error when error(i) became 0 for some i
                    self.adjust_atol(u[nstop - 1])
                elif istate == -7:
                    if solver_name in ("dlsoda","dlsodar"):
                        if iinfo[7] > self.lrw:
                            # Real work array is too short.
                            self.expand_rwork(iinfo[7])
                        elif iinfo[8] > self.liw:
                            # Integer work array is too short.
                            self.expand_iwork(iinfo[8])
                        else:
                            sys.exit(1)
                    else:
                        sys.exit(1)
                else:
                    sys.exit(1)
                if tried > 5:     # prevent endless loop
                    raise ValueError('aborted execution of 5 tries...')
                nstart, istate = nstop, 1
        self.u, self.t = u[:nstop], self.t[:nstop]

        return self.u, self.t

    def advance(self):
        '''
        This function intends to one step forward for linearly implicit solvers
        (Lsodi, Lsodis, Lsoibt) in ODEPACK.

        For these linearly implicit solvers, if extra wrappers are
        added in Fortran code, there are often memory errors. Besides,
        sometimes there are unavoidable errors caused by bugs in the
        Ubuntu/Linux libraries as libc.  To make these solvers more
        reliable on all platforms, this function is used to call
        solvers in ODEPACK (dlsodi, dlsodis, dlsoibt) directly without
        any wrappers in Fortran. However, this would lead to
        efficiency lost with long work arrays as input parameters for
        Fortran code.  In Lsodi, Lsodis and Lsoibt, Solver.solve()
        would be applied to get the desired solution, which will
        direct to this function to step forward.
        '''
        itask = 1   # initial status
        istate = self.ydoti_flag

        solver_name = self.__class__.__name__.lower()
        # Convert class-name to name of subroutine, for example Lsodi -> lsodi

        neq, u, n = self.neq, self.u, self.n
        t, t_next = self.t[n], self.t[n+1]

        res = getattr(self.res_f77, '_cpointer', self.res_f77)
        adda = getattr(self, 'adda_%s_f77' % solver_name)
        if hasattr(adda, '_cpointer'):
            adda = adda._cpointer
        jac = getattr(self, 'jac_%s_f77' % solver_name)
        if hasattr(jac, '_cpointer'):
            jac = jac._cpointer

        tried = 0
        while tried < 5:       # prevent endless loop
            u_new, t, istate, iwork = apply(\
                eval('self._odepack.d%s' % solver_name),\
                (res, adda, jac, neq, u[n].copy(), self.ydoti, t, t_next,
                 self.itol, self.rtol, self.atol, itask, istate, self.iopt,
                 self.rwork, self.lrw, self.iwork, self.liw, self.mf),
                self._extra_args_fortran)
            tried += 1
	    # "istate" indicates the returned status
            if istate >= 1:
                # successful return status
                break
            else:       # Error occurs!
                print self._error_messages[istate] + str(self.rwork[12])
                if istate == -1:    # Increase maximum step-number.
                    self.iwork[5], self.iopt = self.new_stepnr(), 1
                elif istate == -2:  # Multiply tolerance with suggested factor
                    self.tol_multiply(self.rwork[13])
                elif istate == -3:
                    # Illegal input was detected,
                    # before taking any integration steps.
                    if iwork[16] > self.lrw:   # Real work array is too short.
                        self.expand_rwork(iwork[16], expand=True)
                    elif iwork[17] > self.liw: # Integer work array is too short
                        self.expand_iwork(iwork[17], expand=True)
                    else:  # Other abnormal cases with istate returned as -3
                        sys.exit(1)   # Interrupt
                elif istate == -6:
                    # divide-zero-error when error(i) became 0 for some i
                    self.adjust_atol(u_new)
                else:   # Unavoidable interrupts
                    sys.exit(1)  #  Interrupt
                istate = 1
	return u_new

### end of class Odepack ###

class Lsode(Odepack):
    '''

    A Python wrapper of the LSODE (Livermore Solver for Ordinary
    Differential Equations) FORTRAN subroutine.  Basic Solver in
    ODEPACK package from netib.  Solves the initial-value problem for
    stiff or nonstiff systems of first-order ODE, :math:`u' = f(u,t)`
    via Adams or BDF methods (for the nonstiff and stiff cases,
    respectively).  Can generate the Jacobian, or apply a
    user-supplied Jacobian in dense or banded format.

    '''
    quick_description = "LSODE solver for a stiff or nonstiff system"

    _optional_parameters = Odepack._optional_parameters + \
        ['jac_banded', 'jac_banded_f77', 'ml', 'mu', 'jac', 'jac_f77', 'order',
         'f_f77']

    _iwork_index = Odepack._iwork_index.copy()
    _iwork_index[0], _iwork_index[1] = 'ml', 'mu'
    _rwork_index = Odepack._rwork_index.copy()

    _extra_args_fortran = {}

    def adjust_parameters(self):
	"""Properties for new parameters in this solver."""
       # If jac_banded is input in form of jac(u,t,ml,mu),
        # wrap jac_banded to jac_banded_f77 for Fortran code
        self._parameters['jac_banded']['paralist_old'] = 'u,t,ml,mu'
        self._parameters['jac_banded']['paralist_new'] = 't,u,ml,mu'
        self._parameters['jac_banded']['array_order'] = 'Fortran'
        self._parameters['jac_banded']['name_wrapped'] = 'jac_banded_f77'
        # If jac_banded is input in form of jac(t,u,ml,mu),
        # wrap jac_banded_f77 to the general form jac_banded(u,t,ml,mu)
        # for switch_to().
        self._parameters['jac_banded_f77']['paralist_old'] = 't,u,ml,mu'
        self._parameters['jac_banded_f77']['paralist_new'] = 'u,t,ml,mu'
        self._parameters['jac_banded_f77']['name_wrapped'] = 'jac_banded'

        self._parameters['corrector_iter_method']['range'] = range(6)
        self._parameters['corrector_iter_method']['condition-list'] = \
                {'1':[('jac','jac_f77'),],\
                 '4':['ml','mu',('jac_banded','jac_banded_f77')],\
                 '5':('ml','mu')}
        self._parameters['corrector_iter_method']['help'] = """\
Corrector iteration method choice with 6 possible
values:

  0. Functional iteration without any Jacobian
     matrix involved.
  1. Chord iteration with user-supplied full
     Jacobian.
  2. Chord iteration with internally generated
     full Jacobian matrix.
  3. Chord iteration with internally generated
     diagonal Jacobian matrix.
  4. Chord iteration with user-supplied banded
     Jacobian matrix.
  5. Chord iteration with internally generated
     banded Jacobian matrix."""
        Odepack.adjust_parameters(self)

    def set_extra_args(self):
	# ml & mu are required to be extra parameters for banded Jacobian.
        if hasattr(self,'ml') and hasattr(self,'mu'):
            if self.iter_method == 4 and \
                    (not hasattr(self.jac_f77, '_cpointer')):
                self._extra_args_fortran ['jac_extra_args'] = (self.ml,self.mu)
        Odepack.set_extra_args(self)

    def set_iter_method(self):
        if not hasattr(self,'corrector_iter_method'):
            with_ml_mu = hasattr(self, 'ml') and hasattr(self, 'mu')
            with_jac_banded = \
                hasattr(self, 'jac_banded') or hasattr(self, 'jac_banded_f77')
            with_jac_full = hasattr(self, 'jac') or hasattr(self, 'jac_f77')
            if with_ml_mu:
                self.iter_method = 4 if with_jac_banded else 5
            else:
                self.iter_method = 1 if with_jac_full else 2

    def set_jac(self):
        if hasattr(self,'jac_banded_f77') and self.iter_method is 4:
            self.jac_f77 = self.jac_banded_f77

    def set_liw_min(self):
        self.liw_min = 20 if self.iter_method in (0,3,) else 20 + self.neq

    def set_lrw_min(self):
        mf_length = {10:(20,16,0),11:(22,16,1),12:(22,16,1),13:(22,17,0),\
                     14:(22,17,1),15:(22,17,1),20:(20,9,0),21:(22,9,1),\
                     22:(22,9,1),23:(22,10,0),24:(22,10,0),25:(22,10,0)}
        lrw_arg = mf_length[self.mf]
        self.lrw_min = lrw_arg[0] + self.neq*lrw_arg[1] + self.neq**2*lrw_arg[2]
        if self.mf in (24,25):  # user-supplied banded Jacobian matrics
            self.lrw_min += (2*self.ml + self.mu)*self.neq

    def set_iopt(self):
        self.iopt = int(any(self.iwork_in[4:7] > 0) or \
                            any(self.rwork_in[4:7] > 0))


### End of Lsode ###


class Lsoda(Odepack):
    '''

    A Python wrapper of the LSODA FORTRAN subroutine from ODEPACK.
    This subroutine automatically shifts between stiff (BDF) and
    nonstiff (Adams) methods, such that the user does not need to
    determine whether the problem is stiff or not (as is the case when
    using the LSODE subroutine and class :class:`Lsode`). The
    integration always starts with the nonstiff method.  Can generate
    the Jacobian, or apply a user-supplied Jacobian in dense or banded
    format.

    '''
    quick_description ="LSODA solver with stiff-nonstiff auto shift"

    _optional_parameters = Odepack._optional_parameters + \
        ['jac_banded', 'jac_banded_f77', 'ml', 'mu', 'jac', 'jac_f77',
         'max_ordn', 'max_ords', 'f_f77']


    _iwork_index = Odepack._iwork_index.copy()
    _iwork_index[0], _iwork_index[1] = 'ml', 'mu'
    _iwork_index[7], _iwork_index[8] = 'max_ordn', 'max_ords'
    _rwork_index = Odepack._rwork_index.copy()

    _error_messages = Odepack._error_messages.copy()
    _error_messages[-7] = '''
    Length of RWORK or IWORK are too small to proceed, but the integration
    was successful as far as '''

    def adjust_parameters(self):
	"""Properties for new parameters in this solver."""
        # If jac_banded is input in form of jac(u,t,ml,mu),
        # wrap jac_banded to jac_banded_f77 for Fortran code
        self._parameters['jac_banded']['paralist_old'] = 'u,t,ml,mu'
        self._parameters['jac_banded']['paralist_new'] = 't,u,ml,mu'
        self._parameters['jac_banded']['array_order'] = 'Fortran'
        self._parameters['jac_banded']['name_wrapped'] = 'jac_banded_f77'
        # If jac_banded is input in form of jac(t,u,ml,mu),
        # wrap jac_banded_f77 to the general form jac_banded(u,t,ml,mu)
        # for switch_to().
        self._parameters['jac_banded_f77']['paralist_old'] = 't,u,ml,mu'
        self._parameters['jac_banded_f77']['paralist_new'] = 'u,t,ml,mu'
        self._parameters['jac_banded_f77']['name_wrapped'] = 'jac_banded'

        self._parameters['corrector_iter_method']['range'] = [1,2,4,5]
        self._parameters['corrector_iter_method']['condition-list'] = \
                {'1':[('jac','jac_f77'),],
                 '4':[('jac_banded','jac_banded_f77'),'ml','mu'],
                 '5':('ml','mu')}
        self._parameters['corrector_iter_method']['help'] = """\
Jacobian type choice with 4 possible values:

 1. User-supplied full Jacobian matrix
 2. Internally generated full Jacobian (default)
 4. User-supplied banded Jacobian matrix
 5. Internally generated banded Jacobian matrix""",
        Odepack.adjust_parameters(self)

    def set_extra_args(self):
	# ml & mu are required to be extra parameters for banded Jacobian.
        if hasattr(self,'ml') and hasattr(self,'mu'):
            if self.iter_method == 4 and \
                    (not hasattr(self.jac_f77, '_cpointer')):
                self._extra_args_fortran ['jac_extra_args'] = (self.ml,self.mu)
        Odepack.set_extra_args(self)

    def set_iter_method(self):
        if not hasattr(self,'corrector_iter_method'):
            with_ml_mu = hasattr(self, 'ml') and hasattr(self, 'mu')
            with_jac_banded = \
                hasattr(self, 'jac_banded') or hasattr(self, 'jac_banded_f77')
            with_jac_full = hasattr(self, 'jac') or hasattr(self, 'jac_f77')
            if with_ml_mu:
                self.iter_method = 4 if with_jac_banded else 5
            else:
                self.iter_method = 1 if with_jac_full else 2

    def set_jac(self):
        if hasattr(self,'jac_banded_f77') and self.iter_method is 4:
            self.jac_f77 = self.jac_banded_f77

    def set_liw_min(self):
        self.liw_min = 20 + self.neq
        self.mf = self.iter_method

    def set_lrw_min(self):
        self.lrw_min = 20 + 16*self.neq
        if self.iter_method in (1,2):
            self.lrw_min = max(self.lrw_min, 22 + 9*self.neq + self.neq**2)
        else:
            self.lrw_min = max(self.lrw_min, \
                           22 + 10*self.neq + \
                           self.neq*(2*self.ml + self.mu))

    def set_iopt(self):
        self.iopt = int(any(self.iwork_in[5:9]>0) or any(self.rwork_in[4:7]>0))

### End of Lsoda ###

class Lsodar(Odepack):
    '''

    A Python wrapper of the LSODAR subroutine in ODEPACK.
    LSODAR is a variant of LSODE and differs from the latter in
    two ways.

    Quote from the LSODAR source code documentation:
    "(a) It switches automatically between stiff and nonstiff methods.
    This means that the user does not have to determine whether the
    problem is stiff or not, and the solver will automatically choose the
    appropriate method.  It always starts with the nonstiff method.
    (b) It finds the root of at least one of a set of constraint
    functions g(i) of the independent and dependent variables.
    It finds only those roots for which some g(i), as a function
    of t, changes sign in the interval of integration.
    It then returns the solution at the root, if that occurs
    sooner than the specified stop condition, and otherwise returns
    the solution according the specified stop condition."

    The mathematical problem reads

    ..:math::
               u' &= f(u, t),
          g(u,t)  &= 0.

    '''
    quick_description = "LSODAR method with stiff-nonstiff auto shift"

    _optional_parameters = Odepack._optional_parameters + \
        ['jac_banded', 'jac_banded_f77', 'ml', 'mu', 'jac', 'jac_f77',
         'max_ordn', 'max_ords', 'g', 'g_f77', 'ng', 'f_f77']

    _iwork_index = Odepack._iwork_index.copy()
    _iwork_index[0], _iwork_index[1] = 'ml', 'mu'
    _iwork_index[7], _iwork_index[8] = 'max_ordn', 'max_ords'
    _rwork_index = Odepack._rwork_index.copy()

    _error_messages = Odepack._error_messages.copy()
    _error_messages[-7] = '''
    Length of RWORK or IWORK are too small to proceed, but the integration
    was successful as far as '''

    def adjust_parameters(self):
	"""Properties for new parameters in this solver."""
        # If jac_banded is input in form of jac(u,t,ml,mu),
        # wrap jac_banded to jac_banded_f77 for Fortran code
        self._parameters['jac_banded']['paralist_old'] = 'u,t,ml,mu'
        self._parameters['jac_banded']['paralist_new'] = 't,u,ml,mu'
        self._parameters['jac_banded']['array_order'] = 'Fortran'
        self._parameters['jac_banded']['name_wrapped'] = 'jac_banded_f77'
        # If jac_banded is input in form of jac(t,u,ml,mu),
        # wrap jac_banded_f77 to the general form jac_banded(u,t,ml,mu)
        # for switch_to().
        self._parameters['jac_banded_f77']['paralist_old'] = 't,u,ml,mu'
        self._parameters['jac_banded_f77']['paralist_new'] = 'u,t,ml,mu'
        self._parameters['jac_banded_f77']['name_wrapped'] = 'jac_banded'

        # If g is input in form of g(u,t),
        # wrap g to g_f77 for Fortran code.
        self._parameters['g']['paralist_old'] = 'u,t'
        self._parameters['g']['paralist_new'] = 't,u'
        self._parameters['g']['name_wrapped'] = 'g_f77'
        # If g is input in form of g(t,u),
        # wrap g_f77 to the general form g(u,t) for switch_to().
        self._parameters['g_f77']['paralist_old'] = 't,u'
        self._parameters['g_f77']['paralist_new'] = 'u,t'
        self._parameters['g_f77']['name_wrapped'] = 'g'

        self._parameters['corrector_iter_method']['range'] = [1,2,4,5]
        self._parameters['corrector_iter_method']['condition-list'] = \
                {'1':[('jac','jac_f77'),],
                 '4':[('jac_banded','jac_banded_f77'),'ml','mu'],
                 '5':('ml','mu')}
        self._parameters['corrector_iter_method']['help'] = """\
Jacobian type choice with 4 possible values:

1. User-supplied full Jacobian matrix
2. Internally generated full Jacobian (default)
4. User-supplied banded Jacobian matrix
5. Internally generated banded Jacobian matrix""",
        Odepack.adjust_parameters(self)


    def set_extra_args(self):
	# ml & mu are required to be extra parameters for banded Jacobian.
        if hasattr(self,'ml') and hasattr(self,'mu'):
            if self.iter_method == 4 and \
                    (not hasattr(self.jac_f77, '_cpointer')):
                self._extra_args_fortran ['jac_extra_args'] = (self.ml,self.mu)
        Odepack.set_extra_args(self)

    def set_iter_method(self):
        if not hasattr(self,'corrector_iter_method'):
            with_ml_mu = hasattr(self, 'ml') and hasattr(self, 'mu')
            with_jac_banded = \
                hasattr(self, 'jac_banded') or hasattr(self, 'jac_banded_f77')
            with_jac_full = hasattr(self, 'jac') or hasattr(self, 'jac_f77')
            if with_ml_mu:
                self.iter_method = 4 if with_jac_banded else 5
            else:
                self.iter_method = 1 if with_jac_full else 2

    def set_jac(self):
        if hasattr(self,'jac_banded_f77') and self.iter_method is 4:
            self.jac_f77 = self.jac_banded_f77

    def set_liw_min(self):
        self.liw_min = 20 + self.neq
        self.mf = self.iter_method

    def set_lrw_min(self):
        self.lrw_min = 20 + 16*self.neq + 3*self.ng
        if self.iter_method in (1,2):
            self.lrw_min = max(self.lrw_min, 22 + 9*self.neq + self.neq**2 + \
                           3*self.ng)
        else:
            self.lrw_min = max(self.lrw_min, \
                           22 + 10*self.neq + 3*self.ng + \
                           self.neq*(2*self.ml + self.mu))

    def set_iopt(self):
        self.iopt = int(any(self.iwork_in[5:9]>0) or any(self.rwork_in[4:7]>0))

    def solve(self, time_points, terminate=None):
        # ng is numbers of constraint-functions
        if hasattr(self,'g'):
            self.ng = np.asarray(self.g(self.U0,time_points[0])).size
        elif hasattr(self,'g_f77'):
            if not hasattr(self.g_f77,'_cpointer'):
                self.ng = np.asarray(self.g(time_points[0],self.U0)).size
            elif not hasattr(self,'ng'):
                raise ValueError, '''
        Unsufficient input! ng must be specified if g is input as a
        Fortran subroutine. '''
        return Odepack.solve(self,time_points, terminate=terminate)

### End of Lsodar ###

class Lsodes(Odepack):
    """
    A Python wrapper of the LSODES subroutine from ODEPACK.
    LSODES is a variant of LSODE intended for problems where the
    Jacobian is provided via a sparse matrix data structure
    consisting of the parameters ``jac_column``, ``ia``, ``ja``,
    and optionally ``jac_column_f77``.
    """
    quick_description = "LSODES solver for sparse Jacobians"

    _optional_parameters = Odepack._optional_parameters + \
        ['order', 'moss','seth', 'jac_column', 'ia', 'ja', 'jac_column_f77',
         'f_f77']

    _extra_args_fortran = {}
    _error_messages = Odepack._error_messages.copy()
    _error_messages[-7] = '''\
    A fatal error return flag came from the sparse solver CDRV by way of DPRJS
    or DSOLSS (numerical factorization or backsolve).  This should never happen.
    The integration was successful as far as '''


    def adjust_parameters(self):
	"""Properties for new parameters in this solver."""
        # If jac_column is input in form of jac(u,t,j),
        # wrap jac_column to jac_column_f77(t,u,j-1) for Fortran code.
        self._parameters['jac_column']['paralist_old'] = 'u,t,j-1,ia,ja'
        self._parameters['jac_column']['paralist_new'] = 't,u,j,ia,ja'
        self._parameters['jac_column']['name_wrapped'] = 'jac_column_f77'
        # If jac_column is input in form of jac(t,u,j),
        # wrap it to the general form jac_column(u,t,j) for switch_to().
        self._parameters['jac_column_f77']['paralist_old'] = 't,u,j+1'
        self._parameters['jac_column_f77']['paralist_new'] = 'u,t,j'
        self._parameters['jac_column_f77']['name_wrapped'] = 'jac_column'

        self._parameters['moss']['range'] = range(4)
        self._parameters['moss']['condition-list'] = \
                        {1: [('jac_column', 'jac_column_f77'),], \
                         0: ['ia', 'ja']}

        self._parameters['moss']['help'] = """\
Method choice to obtain sparse structure with 3 possible
values:

  0. The user has supplied IA, JA.
  1  The user has supplied JAC_COLUMN and the
     sparse structure will be obtained from NEQ initial
     calls to JAC_COLUMN.
  2. The sparse structure will be obtained from
      NEQ+1 initial calls to F.""",

        self._parameters['corrector_iter_method']['condition-list'] = \
                {'1':[('jac_column','jac_column_f77'),],}
        self._parameters['corrector_iter_method']['help'] = """\
Corrector iteration method choice with 4
possible values:

 0. Functional iteration without any Jacobian
    matrix.
 1. Chord iteration with user-supplied sparse
    Jacobian.
 2. Chord iteration with internally generated
    sparse Jacobian matrix.
 3. Chord iteration with internally generated
    diagonal Jacobian matrix.""",
        Odepack.adjust_parameters(self)

    def set_iter_method(self):
        with_jac_column = hasattr(self,'jac_column') or \
            hasattr(self,'jac_column_f77')
        with_ia_ja = hasattr(self,'ia') and hasattr(self,'ja')
        if not hasattr(self,'moss'):
            if with_ia_ja:
                self.moss = 0
            elif with_jac_column:
                self.moss = 1
            else:
                self.moss = 2
        if not hasattr(self,'corrector_iter_method'):
            self.iter_method = int(with_jac_column)

    def set_iwork_rwork(self):
        '''
        Initialization of work arrays with caculated length and optional inputs.
        In ODEPACK, "iwork" & "rwork" should be initialized with the specific
        optional parameters in all the solvers.
        "liw" & "lrw" represented the length requirement of work arrays.
        Specially, in Dlsodes, ia & ja should be attached to iwork_in.
        '''

        # initialize integer work array (iwork)
        self.iwork_in = [0]*30
        if (self.moss == 0) and hasattr(self,'ia'):
            self.iwork_in += list(self.ia) + list(self.ja)
        nnz = len(getattr(self,'ja',[])) if hasattr(self,'ja') \
                      else self.neq**2/2  # default value
        self.liw_min = len(self.iwork_in)
        for index in self._iwork_index:
            self.iwork_in[index] = getattr(self,self._iwork_index[index],0)

        # calculate the minimum length of float work array (rwork)
        maxord = 5 if self.adams_or_bdf == "bdf" else 12
        self.lrw_min = 20 + self.neq*(maxord + 4)
        lrw_arg = [(0,0,0,0),(0,2,2,9),(0,2,2,10),(self.neq+2,0,0,0)]\
                  [self.iter_method]
        self.lrw_min += lrw_arg[0] + nnz*lrw_arg[1] + self.neq*lrw_arg[2] +\
                  ((nnz + lrw_arg[3]*self.neq)/2)

        # Initializereal input work arrays
        lrw_in = max(self._rwork_index.keys()) + 1
        self.rwork_in = np.zeros(lrw_in,float)
        for index in self._rwork_index:
            self.rwork_in[index] = getattr(self,self._rwork_index[index],0.)

    def set_iopt(self):
        self.iopt = int(self.iwork_in[4:7]>0 or any(self.rwork_in[4:7]>0))

### End of Lsodes ###


class Lsodi(Odepack):
    '''
    A Python wrapper of the LSODI subroutine from ODEPACK, targeting
    ODE problems of the form

    .. math::
               A(u,t) u' & = g(u,t)

    where :math:`A(u,t)` is a square matrix and :math:`g` some function.
    If :math:`A` is singular, this is a differential-algebraic system.
    The user input is in the form of a residual function
    :math:`r = g - As` for some vector :math:`s`. The residual :math:`r`
    is provided by the user function ``res`` (Python) or ``res_f77``
    (FORTRAN).
    '''
    quick_description = "LSODI solver for linearly implicit systems"

    _optional_parameters = Odepack._optional_parameters + \
        ['adda_lsodi', 'adda_banded_lsodi', 'adda_lsodi_f77',
         'adda_banded_lsodi_f77', 'jac_lsodi', 'jac_banded_lsodi',
         'jac_lsodi_f77', 'jac_banded_lsodi_f77', 'res',
         'order', 'ml', 'mu', 'ydoti', 'nsteps', 'res_f77']

    _required_parameters = []

    _iwork_index = Odepack._iwork_index.copy()
    _iwork_index[0], _iwork_index[1] = 'ml', 'mu'
    _rwork_index = Odepack._rwork_index.copy()

    _extra_args_fortran = {}
    _error_messages = Odepack._error_messages.copy()
    _error_messages[-7] = '''\
    User-supplied Subroutine RES set its error flag (IRES = 3) despite repeated
    tries by Lsodi to avoid that condition. Current T is '''
    _error_messages[-8] = '''\
    Lsodi was unable to compute the initial value of du/dt.
    Current T is '''
    _error_messages[-10] = '''\
    User-supplied Subroutine RES signalled Lsodi to halt the
    integration and return (IRES = 2). Current T is '''

    def __init__(self, **kwargs):
        # Different default setting of iteration methods
        with_jac = ('jac_lsodi' in kwargs) or ('jac_lsodi_f77' in kwargs) \
            or ('jac_banded_lsodi' in kwargs) or \
            ('jac_banded_lsodi_f77' in kwargs)
        if 'adams_or_bdf' not in kwargs:
            kwargs['adams_or_bdf'] = 'adams' if with_jac else 'bdf'
        # 'f' is not a mandatory parameter in Lsodi()
        # set f with dummy definition for general constructor
        Odepack.__init__(self, None, **kwargs)

    def adjust_parameters(self):
        self._parameters['corrector_iter_method']['range'] = [1,2,4,5]
        self._parameters['corrector_iter_method']['condition-list'] = \
                         {'1':[('jac_lsodi','jac_lsodi_f77'),
                               ('adda_lsodi','adda_lsodi_f77'),
                               ('res', 'res_f77')],
                          '2':[('adda_lsodi','adda_lsodi_f77'),
                               ('res', 'res_f77')],\
                          '4':['ml','mu', ('res', 'res_f77'),
                               ('adda_banded_lsodi','adda_banded_lsodi_f77'),
                               ('jac_banded_lsodi','jac_banded_lsodi_f77')],\
                          '5':['ml','mu', ('res', 'res_f77'),
                               ('adda_banded_lsodi','adda_banded_lsodi_f77')] }
        self._parameters['corrector_iter_method']['help'] = """\
Choice for the corrector iteration method:

  1. Chord iteration with a user-supplied full
     Jacobian matrix.
  2. Chord iteration with an internally generated
     (difference quotient) full Jacobian. This
     uses neq+1 extra calls to res per dr/du
     evaluation.(Default)
  4. Chord iteration with a user-supplied banded
     Jacobian matrix.
  5. Chord iteration with an internally generated
     banded Jacobian matrix. Using ml+mu+2
     extra calls to res per dr/du evaluation.""",
        Odepack.adjust_parameters(self)

    def set_extra_args(self):
        if self.iter_method == 4:
            self._extra_args_fortran ['jac_extra_args'] = (self.ml,self.mu)
        if self.iter_method > 3:
            self._extra_args_fortran ['adda_extra_args'] = \
                (self.ml,self.mu)
        Odepack.set_extra_args(self)

    def set_iter_method(self):
        if not hasattr(self,'corrector_iter_method'):
            with_banded_jac = hasattr(self,'jac_banded_lsodi_f77') or \
                hasattr(self,'jac_banded_lsodi')
            with_banded_adda = hasattr(self,'adda_banded_lsodi_f77') or \
                hasattr(self,'adda_banded_lsodi')
            with_ml_mu = hasattr(self,'ml') and hasattr(self,'mu')
            with_full_jac = hasattr(self,'jac_lsodi_f77') or \
                hasattr(self,'jac_lsodi')
            with_full_adda = hasattr(self,'adda_lsodi_f77') or \
                hasattr(self,'adda_lsodi')
            if with_ml_mu and with_banded_jac and with_banded_adda:
                self.iter_method = 4
            elif with_ml_mu and with_banded_adda:
                self.iter_method = 5
            elif with_full_jac and with_full_adda:
                self.iter_method = 1
            elif with_full_adda:
                self.iter_method = 2
            else:
                raise ValueError, 'adda must be supplied in Lsodi.'

    def set_jac(self):
        if self.iter_method == 4:
            self.jac_lsodi_f77 = self.jac_banded_lsodi_f77
        if self.iter_method in (4,5):
            self.adda_lsodi_f77 = self.adda_banded_lsodi_f77

    def set_liw_min(self):
        self.liw_min = 20 + self.neq

    def set_lrw_min(self):
        mf_list = [11,12,14,15,21,22,24,25]
        length_args = [(22,16,1),(22,16,1),(22,17,0),(22,17,0),\
                       (22,9,1),(22,9,1),(22,10,0),(22,10,0)]
        lrw_arg = length_args[mf_list.index(self.mf)]
        self.lrw_min = lrw_arg[0] + self.neq*lrw_arg[1] + \
                   self.neq*self.neq*lrw_arg[2] + \
                   self.neq*(2*getattr(self,'ml',0) + \
                             getattr(self,'mu',0))

    def set_iwork_rwork(self):
        # calculate the required length of work arrays
        self.set_liw_min()
        self.set_lrw_min()
        Odepack.check_liwlrw(self)

        self.iwork = np.zeros(self.liw, int)
        self.rwork = np.zeros(self.lrw, float)

        iwork_index, rwork_index = self._iwork_index, self._rwork_index
        # Put optional inputs into work arrays with specified indices.
        for index in iwork_index:
            self.iwork[index] = getattr(self,iwork_index[index], 0)
        for index in rwork_index:
            self.rwork[index] = getattr(self, rwork_index[index], 0.)

    def set_iopt(self):
        self.iopt = int(any(self.iwork[4:7]>0) or any(self.rwork[4:7]>0))

    def set_dummy_functions(self):
        if getattr(self, 'adda_lsodi_f77', None) is None:
            self.adda_lsodi_f77 = lambda x,y,z: 0.
        if getattr(self, 'jac_lsodi_f77', None) is None:
            self.jac_lsodi_f77 = lambda x,y,z: 0.

    def solve(self, time_points, terminate=None):
        # Call Solver.solve(), which will direct to Odepack.advance()
        # to step forward.
        return Solver.solve(self, time_points, terminate=terminate)


### End of Lsodi ###

class Lsodis(Odepack):
    '''
    A Python wrapper of the LSODIS subroutine from ODEPACK.
    This subroutine is a variant of LSODI, intended for stiff problems
    in which the matrix A and the Jacobian of the residual wrt the
    unknown functions have a sparse structure.
    '''
    quick_description = "LSODIS solver for linearly implicit sparse systems"

    _optional_parameters = Odepack._optional_parameters + \
        ['jac_lsodis', 'moss', 'ia', 'ja', 'ic', 'jc',
         'ydoti', 'order']

    _required_parameters = ['res', 'adda_lsodis']
    # Do not support Fortran subroutines

    _extra_args_fortran = {}
    _error_messages = Odepack._error_messages.copy()
    _error_messages[-7] = '''\
    User-supplied Subroutine RES set its error flag (IRES = 3) despite repeated
    tries by Lsodis to avoid that condition. Current T is '''
    _error_messages[-8] = '''\
    Lsodis was unable to compute the initial value of du/dt.
    Current T is '''
    _error_messages[-9] = '''\
    A fatal error return-flag came from the sparse solver CDRV.
    This should never happen. Current T is '''
    _error_messages[-10] = '''\
    User-supplied Subroutine RES signalled Lsodis to halt the
    integration and return (IRES = 2). Current T is '''

    def __init__(self,**kwargs):
        # Different default setting of iteration methods
        with_jac = 'jac_lsodis' in kwargs
        if 'adams_or_bdf' not in kwargs:
            kwargs['adams_or_bdf'] = 'adams' if with_jac else 'bdf'
        # 'f' is not a mandatory parameter in Lsodis()
        # set f with dummy definition for general constructor
        Odepack.__init__(self, None, **kwargs)

    def adjust_parameters(self):
        self._parameters['corrector_iter_method']['range'] = [1,2]
        self._parameters['corrector_iter_method']['condition-list'] = \
                         {'1':['jac_lsodis',],}
        self._parameters['corrector_iter_method']['help'] = """\
Choice for the corrector iteration method:

  1. Chord iteration with a user-supplied sparse
     Jacobian matrix.
  2. Chord iteration with an internally generated
     (difference quotient) sparse Jacobian matrix.
     This uses extra calls to "res" per dr/du
     evaluation. (Default)""",
        self._parameters['moss']['range'] = range(5)
        self._parameters['moss']['condition-list'] = \
                         {'1':['jac_lsodis',],
                          '0':['ia','ja','ic','jc'],
                          '3':['jac_lsodis','ia','ja'],
                          '4':['ia','ja']}
        self._parameters['moss']['help'] = """\
Method to obtain sparse structure of Jacobian matrix.
moss has 5 possible values:

 0. The user has supplied IA, JA, IC, and JC.
 1. The user has supplied JAC and the
    structure will be obtained from NEQ initial
    calls to JAC and NEQ initial calls to ADDA.
 2. The structure will be obtained from NEQ+1
    initial calls to RES and NEQ initial calls to ADDA.
 3. like MOSS = 1, except user has supplied IA and JA.
 4. like MOSS = 2, except user has supplied IA and JA""",
        Odepack.adjust_parameters(self)

    def set_iter_method(self):
        with_jac = hasattr(self, 'jac_lsodis')
        with_ia_ja = hasattr(self, 'ia') and  hasattr(self, 'ja')
        with_ic_jc = hasattr(self, 'ic') and hasattr(self, 'jc')
        if not hasattr(self, 'corrector_iter_method'):
            self.iter_method = 1 if with_jac else 2
        if not hasattr(self, 'moss'):
            if with_ia_ja and with_ic_jc:
                self.moss = 0
            elif with_ia_ja:
                self.moss = 3 if with_jac else 4
            else:
                self.moss = 1 if with_jac else 2

    def set_iwork_rwork(self):
        '''
        Initialization of work arrays with caculated length and optional inputs.
        In ODEPACK, "iwork" & "rwork" should be initialized with the specific
        optional parameters in all the solvers.
        "liw" & "lrw" represented the length requirement of work arrays.
        Specially, in Lsodis, (ia, ja, ic & jc) should be attached to
        iwork.
        '''
        # initialize integer work array (iwork)
        self.iwork = [0]*30
        if (self.moss in (0, 3, 4)) and hasattr(self,'ia'):
            self.iwork += list(self.ia) + list(self.ja)
            if (self.moss == 0) and hasattr(self,'ic'):
                self.iwork += list(self.ic) + list(self.jc)
        nnz = len(getattr(self,'ja',[])) if hasattr(self,'ja') \
                      else self.neq**2/2  # default value
        self.liw_min = len(self.iwork)

        # calculate the length of  float work array (rwork)
        mf_list = [11,12,21,22]
        length_args = [(20,18,9),(20,18,10),(20,11,9),(20,11,10)]
        lrw_arg = length_args[mf_list.index(self.mf % 100)]
        self.lrw_min = 2*nnz + lrw_arg[0] + self.neq*lrw_arg[1] + \
                       (nnz + lrw_arg[2]*self.neq)/2

        # Check lrw>lrw_min, liw>liw_min if lrw/liw are specified by users
        Odepack.check_liwlrw(self)
        if self.liw > self.liw_min:
            self.iwork += [0]*(self.liw - self.liw_min)

        # Read in optional inputs
        for index in self._iwork_index:
            self.iwork[index] = getattr(self,self._iwork_index[index],0)

        # Initialize real work array
        self.rwork = np.zeros(self.lrw, float)
        for index in self._rwork_index:
            self.rwork[index] = getattr(self,self._rwork_index[index],0.)

    def set_iopt(self):
        self.iopt = int(any(self.iwork[4:7])>0 or any(self.rwork[4:7])>0)

    def set_dummy_functions(self):
        for name in ('jac_lsodis_f77', 'adda_lsodis_f77'):
            if getattr(self, name, None) is None:
                setattr(self, name, lambda x,y,z,i,j,k:0.)

    def solve(self, time_points, terminate=None):
        # Call Solver.solve(), which will direct to Odepack.advance()
        # to step forward.
        return Solver.solve(self, time_points, terminate=terminate)

### End of class Lsodis ###

class Lsoibt(Odepack):
    '''
    A Python wrapper of the LSOIBT subroutine from ODEPACK.
    This subroutine is a variant of LSODI for the case where the
    matrices :math:`A`, :math:`dg/du`, and :math:`d(As)/du` are all
    block tridiagonal.
    '''
    quick_description = "LSOIBIT solver for linearly implicit block tridiag systems"

    _optional_parameters = Odepack._optional_parameters + \
        ['jac_lsoibt', 'ydoti', 'order']

    _required_parameters = ['res', 'adda_lsoibt', 'mb', 'nb']
    # Do not support Fortran subroutines

    _iwork_index = Odepack._iwork_index.copy()
    _iwork_index[0], _iwork_index[1] = 'mb', 'nb'
    _rwork_index = Odepack._rwork_index.copy()

    _extra_args_fortran = {}
    _error_messages = Odepack._error_messages.copy()
    _error_messages[-7] = '''\
    User-supplied Subroutine RES set its error flag (IRES = 3) despite repeated
    tries by Lsoibt_ODEPACK to avoid that condition. Current T is '''
    _error_messages[-8] = '''\
    Lsoibt_ODEPACK was unable to compute the initial value of du/dt.
    Current T is '''
    _error_messages[-10] = '''\
    User-supplied Subroutine RES signalled Lsoibt_ODEPACK to halt the
    integration and return (IRES = 2). Current T is '''

    def __init__(self,**kwargs):
        # Different default setting of iteration methods
        with_jac = 'jac_lsoibt' in kwargs
        if 'adams_or_bdf' not in kwargs:
            kwargs['adams_or_bdf'] = 'adams' if with_jac else 'bdf'
        # 'f' is not a mandatory parameter in Lsoibt()
        # set f with dummy definition for general constructor
        Odepack.__init__(self, None, **kwargs)

    def adjust_parameters(self):
        self._parameters['corrector_iter_method']['range'] = [1,2]
        self._parameters['corrector_iter_method']['condition-list'] = \
                         {'1':('jac_lsoibt',),}
        self._parameters['corrector_iter_method']['help'] = """\
Choice for the corrector iteration method:

  1. Chord iteration with a user-supplied block-
     tridiagonal Jacobian matrix.
  2. Chord iteration with an internally generated
     (difference quotient) block-tridiagonal
     Jacobian matrix. This uses 3*mb+1 calls to
     ``res`` per dr/du evaluation. (Default).""",
        Odepack.adjust_parameters(self)

    def set_iter_method(self):
        with_jac = hasattr(self, 'jac_lsoibt')
        if not hasattr(self, 'corrector_iter_method'):
            self.iter_method = 1 if with_jac else 2

    def set_liw_min(self):
        self.liw_min = 20 + self.neq

    def set_lrw_min(self):
        mf_list = [11,12,21,22]
        length_args = [(22,16,3),(22,16,3),(22,9,3),(22,9,3)]
        lrw_arg = length_args[mf_list.index(self.mf)]
        self.lrw_min = lrw_arg[0] + self.neq*lrw_arg[1] + \
            3*self.neq*self.mb*lrw_arg[2]

    def set_iwork_rwork(self):
        # calculate the required length of work arrays
        self.set_liw_min()
        self.set_lrw_min()
        Odepack.check_liwlrw(self)

        self.iwork = np.zeros(self.liw, int)
        self.rwork = np.zeros(self.lrw, float)

        iwork_index, rwork_index = self._iwork_index, self._rwork_index
        # Put optional inputs into work arrays with specified indices.
        for index in iwork_index:
            self.iwork[index] = getattr(self,iwork_index[index], 0)
        for index in rwork_index:
            self.rwork[index] = getattr(self, rwork_index[index], 0.)

    def set_iopt(self):
        self.iopt = int(any(self.iwork[4:7]>0) or any(self.rwork[4:7]>0))

    def set_dummy_functions(self):
        if getattr(self, 'jac_lsoibt_f77', None) is None:
            self.jac_lsoibt_f77 = lambda x,y,z:(0., 0., 0.)
        if getattr(self, 'adda_lsoibt_f77', None) is None:
            self.adda_lsoibt_f77 = lambda x,y,z,i,j:(0., 0., 0.)

    def validate_data(self):
        if not Odepack.validate_data(self):
            return False
        if self.mb*self.nb != self.neq:
                raise ValueError,'''
    The requirement for block size (mb,nb) are: mb>=1, nb>=4, mb*nb=neq=%d.
    Your block size are (%d,%d). ''' % (self.neq, mb, nb)
        return True

    def solve(self, time_points, terminate=None):
        # Call Solver.solve(), which will direct to Odepack.advance()
        # to step forward.
        return Solver.solve(self, time_points, terminate=terminate)

### end of class Lsoibt ###

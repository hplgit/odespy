# Author: Liwei Wang, Hans Petter Langtangen

'''
This module contains the base class ``Solver``
package and the implementations of many subclasses.


How to implement a new solver
=============================

Simple explicit solver
----------------------

A new non-adaptive time-stepping scheme can often be implemented by
just subclassing ``Solver`` and writing the ``advance`` method
to define the time-stepping scheme. Here is how the 4th-order
Runge-Kutta is implemented::

    class RK4(Solver):
        quick_description = "Explicit 4th-order Runge-Kutta method"

        def advance(self):
            u, f, n, t = self.u, self.f, self.n, self.t
            dt = t[n+1] - t[n]
            dt2 = dt/2.0
            K1 = dt*f(u[n], t[n])
            K2 = dt*f(u[n] + 0.5*K1, t[n] + dt2)
            K3 = dt*f(u[n] + 0.5*K2, t[n] + dt2)
            K4 = dt*f(u[n] + K3, t[n] + dt)
            u_new = u[n] + (1/6.0)*(K1 + 2*K2 + 2*K3 + K4)
            return u_new

The ``quick_description`` string is needed for the class to appear
in the package's table of contents.

Implicit solver
---------------

See examples in the ``solver.py`` file, class ``BackwardEuler``,
for instance.

Adaptive solver
---------------

See examples in the ``solver.py`` file, class ``RK43``,
for instance.

Wrapping other packages
-----------------------

This section is for developers who intend to wrap an existing
package and integrate with the present one..

If the original package is not written in Python language, developers
need to apply some specific tools (like F2PY, Cython, or SWIG) to create an
extension module to make the package accessible from Python
code.

On the other hand, if the original package is also a Python software,
developers need to install and import (in ``initialize``) the desired
package as a Python module and just write a class in the ``Solver``
hierarchy that calls the proper functions in the package. See the
classes ``SymPy_odefun`` and ``Ode_scipy`` (and its subclasses).

By an attempt to import these necessary modules (often set in method
initialize()), we can check whether the necessary dependencies are
installed properly.

Definition of parameters and their properties
---------------------------------------------

Each solver has a set of specific parameters depending on its
underlying method. For example, adaptive solvers will be more
likely to apply (many) attributes for step-size control.
When integrating a new method, first search in the ``_parameters``
dictionary in ``solver.py`` if some parameters can be reused
for the new solver. New parameters are added to the ``_parameters``
dictionary, either in the ``solver.py`` or by importing ``solver``
and updating ``solver._parameters``, see ``rkc.py`` for an example.

Each solver class lists the required and optional parameters
it needs in the class variables ``_optional_parameters`` and
``_required_parameters``. Very often these lists are inherited,
or just a few new parameters are added to the list already
defined in the superclass.

Sometimes values of parameters (or other properties) need to be
changed in a solver, e.g., because there are certain relations between
various parameters. Appropriate adjustments and checks are done in the
method ``initialize_for_solve``, which is called in the beginning of
the "solve" process (before any computations take place).  Many
classes provide examples on this, e.g., class ``RKC`` in ``rkc.py``.
This particular class shows how to generate input parameters
to the Fortran code that are not given by the user, but automatically
derived in the Python code from other data.

Another method that is called from ``solve`` is ``validate_data``. The
solver class can use this method to check for consistency of data
structures before starting the numerical computations.  As an example,
the class ``Lsodes`` in the ``odepack.py`` file imposes relations
between the input data: two input integer arrays ``ia`` and ``ja``
must be input simultaneously, ``len(ia) == neq + 1``, and ``ia[neq] =
len(ja)``. These checks are done in Python before calling the Fortran
solver.


The solve and advance methods
-----------------------------

Simple methods can just implement ``advance`` to bring the solution
one step forward in time. The existing ``solve`` method derived from
the superclass is general enough to administer the whole solution
process.

Adaptive solvers will typically think of ``advance`` as bringing the
solution to the next user-desired time level, using an unknown set of
smaller time steps whose sizes must be determined. Then there will
hence be a time loop within ``advance``, while the outer time loop in
``solve`` takes care of the stepping between the user-desired time
levels. (The simplest methods just takes one step between the
user-desired time levels.)

When wrapping solvers in external software, it is occasionally not
feasible to implement ``advance`` because one wants to rely on the
software's ability to solve the whole ODE problem. Then it is more
natural to write a new ``solve`` method (using ``Solver.solve`` as
model) and call up the solve functionality in the external
software. Class ``SymPy_odefun`` provides an example. On the other
hand, the wrappers for the ``scipy`` solvers (``vode`` for instance)
applies ``solve`` from the present package and the ``scipy`` functions
for doing one (adaptive) time step, called in the ``advance`` method.

'''

import pprint, sys, os, inspect
import numpy as np

# Collection of all possible parameters in all solvers in this package
_parameters = dict(

    f = dict(
        help='Right-hand side ``f(u,t)`` defining the ODE',
        type=callable),

    f_args = dict(
        help='Extra positional arguments to f: ``f(u, t, *f_args, **f_kwargs)``',
        type=(tuple, list, np.ndarray),
        default=()),

    f_kwargs = dict(
        help='Extra keyword arguments to f: ``f(u, t, *f_args, **f_kwargs)``',
        type=dict,
        default={}),

    complex_valued = dict(
        help='True if f is complex valued',
        default=False,
        type=bool),

    jac = dict(
        help='Jacobian of right-hand side function f (df/du)',
        default=None,
        type=callable),

    jac_args = dict(
        help='Extra positional arguments to jac: ``jac(u, t, *jac_args,'\
             '**jac_kwargs)``',
        type=(tuple,list),
        default=()),

    jac_kwargs = dict(
        help='Extra keyword arguments to jac: ``jac(u, t, *jac_args,'\
             '**jac_kwargs)``',
        type=dict,
        default={}),

    h_in_fd_jac = dict(
        help='h in finite difference approximation of the Jacobian',
        default=1E-4,
        type=float),

    verbose = dict(
        help='Integer reflecting output of intermediate quantities',
        default=0,
        type=int),

    start_method = dict(
        help='Method for the first steps in multi-step solvers',
        default='RK2',
        type=str),

    nonlinear_solver = dict(
        help='Newton or Picard nonlinear solver',
        default='Picard',
        type=str,
        range=('Newton', 'Picard')),

    eps_iter = dict(
        help='Max error measure in nonlinear solver',
        default=1E-4,
        type=float),

    max_iter = dict(
        help='Max no of iterations in nonlinear solver',
        default=25,
        type=int),

    g = dict(
        help='Constraint function of (u, t) in differential-algebraic systems.',
        type=callable),

    ng = dict(
        help='No of components in constraint function g.',
        type=int),

    theta = dict(
        help='Weight in [0,1] used for "theta-rule" finite difference approx.',
        default=0.5,
        type=(int,float),
        range=[0, 1]),

    # Parameters for adaptive methods

    atol = dict(
        help='absolute tolerance for solution',
        type=(float,list,tuple,np.ndarray),
        default=1E-8),

    rtol = dict(
        help='relative tolerance for solution',
        type=(list,tuple,np.ndarray,float),
        default=1E-6),

    min_step = dict(
        help='Minimum step size for an adaptive algorithm',
        type=float),

    max_step = dict(
        help='Maximum step size for an adaptive algorithm',
        type=float),

    first_step = dict(
        help='Suggested first time step size for an adaptive algorithm',
        type=float),

    solver = dict(
        help='Name of solver class in solvers that need an extra solver '\
             '(e.g., AdaptiveResidual)',
        default='RK4',
        type=str),

    butcher_tableau = dict(
        help='2d-array which contains the butcher table for user-supplied '\
             'Runge-Kutta method. (n,n) array for 1-level Runge-Kutta '\
             'methods.(n+1,n) array for 2-level Runge-Kutta methods.',
        type=np.ndarray),

    # vode parameters
    adams_or_bdf = dict(
        help='Method in vode or solvers in odepack: "adams" or "bdf"',
        type=str,
        default='adams',
        range=['adams', 'bdf']),

    nsteps = dict(
        help='Max no of internal solver steps per time step',
        type=int,
        default=1000),

    order = dict(
        help='Maximum order used by the integrator '\
             '(<= 12 for "adams", <= 5 for "bdf")',
        type=int,
        default=4),

    # beta, ifactor and dfactor are intended for adaptive Dormand&Prince
    # methods like dopri5 or dop853 in scipy
    beta = dict(
        help='Beta argument for stabilized step size control in '\
             'Dormand&Prince methods from scipy',
        type=float),

    ifactor = dict(
        help='Maximum factor for increasing the step size',
        type=float,
        default=2),

    dfactor = dict(
        help='Maximum factor for decreasing the step size',
        type=float,
        default=0.5),

    safety = dict(
        help='Safety factor on new step selection',
        default=0.9,
        type=float),

    # odelab parameters
    odelab_solver = dict(
        help='Name of Solver class in odelab',
        default='RungeKutta34',
        type=str),

    # Vode_PyDS parameters
    init_step = dict(
        help='Fixed step size for time mesh.',
        type=float),

    strictdt = dict(
        help='Uniform time mesh vs exact dt spacings',
        type=bool,
        default=True),

    stiff = dict(
        help='Boolean flag to indicate stiffness.',
        type=bool),

    use_special = dict(
        help='Switch for using special times',
        type=bool),

    specialtimes = dict(
        help='List of special times to use during iteration',
        type=lambda float_seq: np.asarray(map(lambda x: \
                     isinstance(x, float),float_seq)).all()),

    ode_method = dict(
        help='solver type: "adams" or "bdf"',
        alias='method',   # Different name in scipy.ode.vode
        type=str, default='adams',
        range=('adams','bdf')),

    relaxation = dict(
        help='relaxation argument (r): new_solution = r*solution + '\
             '(1-r)*old_solution',
        default=1.0, type=float),


    # parameters for Jacobian
    jac_banded = dict(
        help='Banded Jacobian (df/du) of right-hand side'\
             'function f. Arguments: u, t, ml, mu. '\
             'Used in Lsode, Lsoda, Lsodar.',
        type=callable),

    jac_column = dict(
        help='A column of the Jacobian (df/du) matrix.'\
             'Arguments: u, t, ia, ja. Used in Lsodes.',
        type=callable),

    jac_constant = dict(
        help='Flag to show whether Jacobian is constant, 0 (false) or 1 (true)',
        default=0,
        type=int),

    # parameters for linearly implicit ODE solvers: Lsodi, Lsoibt, Lsodis
    res = dict(
        help='User-supplied function to calculate the residual vector,'\
             'defined by ``r = g(t,y) - A(t,y) * s``.'\
             'Used in Lsodi, Lsoibt, Lsodis',
        type=callable),

    ydoti = dict(
        help='Real array for the initial value of dy/dt.',
        type=(list,tuple,np.ndarray),
        extra_check=lambda float_seq: np.asarray(map(lambda x: \
             isinstance(x, float),float_seq)).all(),
        default = []),

    # ja, ia, jc & ic are used to describe the sparse structure
    # of matrices
    ja = dict(
        help='Integer array containing the row indices of nonzero entries '
              'in a sparse matrix (CRS storage scheme)'.
        type=(list,tuple,np.ndarray)),

    ia = dict(
        help='Integer array containing info where the different rows '\
              'of a sparse matrix start (CRS storage scheme)'.
        type=(list,tuple,np.ndarray)),

    # ml, mu describe banded Jacobian matrix.
    ml = dict(
        help='Lower half-bandwidth of banded jacobian matrix',
        type=int),

    mu = dict(
        help='Upper half-bandwidth of banded jacobian matrix',
        type=int),

    # mb, nb describe the block-tridiagonal form of matrix.
    # Used in Lsoibt.
    mb = dict(
        help='Block size,  mb>=1, mb*nb = neq (number of equations).',
        type=int,
        extra_check=lambda x: x>=1),

    nb = dict(
        help='Number of blocks in the main diagonal. nb>=4',
        type=int,
        extra_check=lambda x:x>=4),

    # Odepack parameters
    seth = dict(
        help='Element threshhold for sparsity determination.',
        default=0,
        type=int),

    iter_method = dict(
        help='Corrector iteration method choice',
        type=int),

    lrw = dict(
        help='Length of real work array.',
        type=int),

    liw = dict(
        help='Length of integer work array, similiar as <lrw>.',
        type=int),

    method_order = dict(
        help='Method order for user-defined method if known.'\
             'A integer for 1-level methods, or a pair of   '\
             'integer for 2-levels methods.',
        type=(int,tuple,list,np.ndarray)),

    moss = dict(
        help=' Method to obtain sparse structure of Jacobian.',
        type=int),

    max_hnil = dict(
        help='Maximum no of warning messages to be printed.',
        type=int),

    max_ordn = dict(
        help='Maximum order in nonstiff methods. ',
        type=int),

    max_ords = dict(
        help='Maximum order in stiff methods. ',
        type=int),

   # Fortran versions of f, jac, g (can be used when solver is in Fortran)
    f_f77 = dict(
        help='Intended to supply a user-supplied Fortran subroutine as f.',
        type=callable),

    g_f77 = dict(
        help='Intend to supply a Fortran subroutine as g.',
        type=callable),

    jac_f77 = dict(
        help='Intend to supply a Fortran subroutine as jac.',
        type=callable),

    myadvance = dict(
        help='User supplied function to advance current solution'\
             ' one step forward. See documents of class MySolver.',
        type=callable),

    )


def _format_parameters_table(parameter_names, fixed_width=None):
    """
    Make a table of parameter names and their descriptions.
    The parameter_names list contains the names and the
    descriptions are taken from _parameters[name]['help'].
    max_name_length is the width of the first column, taken
    as the longest name in parameter_names if not specified.

    The table is formatted as a simple reST table with headings
    ("Name" and "Description") and three horizontal lines.
    """
    import textwrap
    max_line_width = 71

    if fixed_width is None:
        max_name_length = max([len(name) \
                               for name in parameter_names + ['Name']])
        c1 = max_name_length + 1     # width of column 1
        c2 = (max_line_width - c1)   # width of column 2
    else:
        c1, c2 = fixed_width

    s = ''  # string to be returned (table)
    hrule = '='*c1 + ' ' + '='*c2 + '\n'
    heading = 'Name' + ' '*(c1-3) + 'Description\n'
    s += hrule + heading + hrule

    for name in parameter_names:
        s += '%%-%ds' % (c1+1) % name
        if name in _parameters:
            text = _parameters[name]['help']
            if 'default' in _parameters[name]:
                text += ' (default: %s)' % str(_parameters[name]['default'])
            if fixed_width is not None:
                line_no = len(text)/c2 if len(text) % c2 == 0 \
                    else len(text)/c2 + 1
                # Spilt text with fixed width
                text = [text[i*49:(i+1)*49] for i in range(line_no)]
            else:
                # List of wrapped lines
                text = textwrap.wrap(text, c2, break_long_words=False)
            for i in range(1, len(text)):   # add initial space for line 2, ...
                text[i] = ' '*(c1+1) + text[i]
            text = '\n'.join(text)
            s += text
        s += '\n'

    s += hrule
    return s

def table_of_parameters(classname, fixed_width=None):
    """
    Return a table (in reST format) of the required parameters
    in a class and a table of the optional parameters.
    The returned string is typially appended to the doc string of
    a solver class so that the user can easily see which parameters
    that must and can be provided.
    """
    req_prm = getattr(classname, '_required_parameters')
    opt_prm = getattr(classname, '_optional_parameters')
    for name in opt_prm:
        if not name in _parameters:
            print 'Parameter "%s" used in class %s is not registered in _parameters.' % (name, classname.__name__)
            print 'Do that before proceeding.'
            sys.exit(1)

    s = """
Required input arguments:

""" + _format_parameters_table(req_prm, fixed_width=fixed_width) + \
"""
Optional input arguments:

""" + _format_parameters_table(opt_prm, fixed_width=fixed_width)
    # Add indent:
    indent = 4
    newlines = [' '*indent + line for line in s.splitlines()]
    s = '\n'.join(newlines)
    return s

def typeset_toc(toc):
    toc = sorted(toc)
    column1_width = max([len(classname) for classname, descr in toc])
    column2_width = max([len(descr)     for classname, descr in toc])
    hrule = '='*(column1_width + 1) + ' ' + '='*(column2_width)
    def line(name, descr):
        return '%%-%ds %%s' % (column1_width+1) % (name, descr)
    lines = [hrule, line('Classname', 'Short description'), hrule] + \
            [line(name, descr) for name, descr in toc] + [hrule]
    return '\n'.join(lines)



class Solver:
    """
    Superclass for numerical methods solving ODE problem

      u'(t) = f(u, t),  u(0) = U0

    where u and U0 are scalars (for scalar ODEs) or vectors
    (for systems of ODEs).

    Attributes stored in this class:

    =========  ========================================================
    Name       Description
    =========  ========================================================
    u          array of point values of the solution function
    t          array of time values: u[i] corresponds to t[i]
    n          the most recently computed solution is u[n+1]
    f          function wrapping the user's right-hand side f(u, t),
               used in all algorithms
    users_f    the user's original function implementing f(u, t)
    PRM        an attribute for each optional and required parameter
    =========  ========================================================

    """

    _required_parameters = ['f',]
    _optional_parameters = ['f_args', 'f_kwargs', 'complex_valued', 'verbose']

    def __init__(self, f, **kwargs):
        """
        ``f`` is the right-hand side function of the ODE u' = f(u,t).
        The legal keyword arguments (in ``kwargs``) are documented in
        the tables in the doc string of this class. The ``f`` function
        must return a ``float`` or ``complex`` object in case of a
        scalar ODE and a list or array of ``float`` or ``complex`` objects
        in case of a system of ODEs.

        This constructor makes a dictionary ``self._parameters``
        holding all the required and optional parameters for this solver
        (fetched from the global ``_parameters`` dictionary in this module).
        The method ``adjust_parameters`` (implemented in subclasses)
        is called to adjust default parameter settings if needed.
        Then all keys in ``self._parameters`` become class attributes,
        filled with default values. Thereafter, all keyword arguments
        (in ``kwargs``) with ``None`` as value are removed as keyword
        arguments. The next step is to call ``set(**kwargs)``, i.e.,
        use the keyword arguments to modify the values of the attributes
        that represent the parameters in this solver. Finally, the
        constructor calls the method ``initialize`` (to be implemeneted
        in subclasses, e.g., for importing necessary modules for the solver).

        Instead of supplying keyword arguments to this constructor, the
        user can at any time call the ``set`` method with keyword
        arguments in order to specify parameters.
        """

        # self._parameters is the union of optional and required parameters
        # for the class. self._parameters contains all the
        # legal parameters the user of the class can set.
        self._parameters = dict(
            (key, value.copy()) for key, value in _parameters.items()
            if key in self._optional_parameters or \
               key in self._required_parameters
            )

        # Compile user-supplied functions if they are supplied
        # as multi-line strings in Fortran code
        f, kwargs = self.compile_string_functions(f, **kwargs)

        # Adjust self._parameters
        self.adjust_parameters()

        # Set default values for all parameters, remove all parameters
        # with value None, and then apply set() to all the user-provided
        # parameters in kwargs

        for name in self._parameters:
            if 'default' in self._parameters[name]:
                setattr(self, name, self._parameters[name]['default'])

        nones = [name for name in kwargs.keys() if kwargs[name] is None]
        for name in nones:
            del kwargs[name]

        self.set(**kwargs)

        # Wrap user-supplied f with extra arguments
        if f is not None:
            self.users_f = f  # stored in case it is handy to have
            if not callable(f):
                raise TypeError('f is %s, not a callable function' % type(f))
            # For ODE systems, f will often return a list, but
            # arithmetic operations with f in numerical methods
            # require that f is an array. Let self.f be a function
            # that first calls f(u,t) and then ensures that the
            # result is an array (without imposing any type - if
            # U0 has integers it is detected and converted to floats
            # to ensure float results from f).
            if 'f_args' in self._optional_parameters:
                self.f = lambda u, t:  \
                    np.asarray(f(u, t, *self.f_args, **self.f_kwargs))
            else:
                self.f = lambda u, t: np.asarray(f(u,t))

        # Subclass-specific initialization
        self.initialize()


    def compile_string_functions(self, f, **kwargs):
        """
        Compile functions which are supplied as Fortran strings.
        """
        str_funcs = dict(
            (func,kwargs[func]) for func in kwargs
            if isinstance(kwargs[func], str) and \
               (func in _parameters) and \
               (_parameters[func]['type'] is callable))
        if isinstance(f, str):
            str_funcs['f'] = f

        try:
            os.remove('_callback.so')
        except:
            pass

        string_to_compile = '\n'.join(str_funcs.values())
        if string_to_compile is not '':
            try:
                # Compile these functions together into module callback
                from numpy import f2py
                f2py.compile(string_to_compile, modulename='_callback', \
                             verbose=False)
                import _callback
                for func in str_funcs.keys():
                    if func == 'f':
                        f = _callback.f
                    else:
                        kwargs[func] = getattr(_callback,func)
            except:
                raise ValueError, '''
           F2py failed to compile input string (=\n%s)
           to be callable functions (%s).''' \
                 % (string_to_compile, str_funcs.keys())
        return f, kwargs

    def adjust_parameters(self):
        '''
        This method allows subclasses to adjust (modify or add)
        entries in the self._parameters dictionary.
        The method is called from the constructor.

        Further adjustments of self._parameters can be done in
        initialize_for_solve when all data for the solver are available.
        '''
        # Define start_method and method here since the range depends
        # on return value of list_all_solvers().
        _parameters['start_method']['range'] = \
            _parameters['solver']['range'] = list_all_solvers()
        return None


    def initialize(self):
        """
        Subclass-specific initialization. Called from constructor.
        Typical use: import modules needed in methods in the class
        and provide error messages if modules are not installed.
        """
        return None


    def set(self, strict=False, **kwargs):
        """
        Assign values to one or more parameters, specified as keyword
        arguments.

        The legal parameters that can be set are contained in the dict
        self._parameters.

        If strict is true, only registered parameter names are accepted,
        otherwise unregistered parameters are ignored.

        The provided parameters (keyword arguments in kwargs) are
        first checked for legal type and legal range.

        Types and ranges of attributes are defined in self._parameters,
        which is initialized with default settings and optionally
        modified in the adjust_parameters method.
        """

        # if new string functions are supplemented
        f_dummy, kwargs = self.compile_string_functions(None, **kwargs)

        # Check for invalid names in kwargs
        kwargs_copy = kwargs.copy()
        for name in kwargs_copy:
            if name not in self._parameters:
                # invalid name
                if strict:
                    raise ValueError('set: parameter %s=%s has illegal name' % \
                                     (name, kwargs[name]))
                del kwargs[name]
            elif kwargs[name] is None:
                del kwargs[name]
                if hasattr(self, name):    # Remove this attribute
                    del self.__dict__[name]

        self.check_input_types(**kwargs)  # all values of right type?
        self.check_input_range(**kwargs)  # all values of right range?

        # Run extra check functions if specified
        self.check_extra(**kwargs)

        # Tests on right name/type/range were successful (if we come here)
        for name in kwargs:
            setattr(self, name, kwargs[name])

        # all conditional parameters are supplied?
        self.check_conditional_parameters()

    def check_input_types(self, **kwargs):
        """Check whether all existing inputs are of right specified type."""

        parameters = self._parameters
        arg_type_list = [(name,parameters[name]['type'],kwargs[name]) \
                           for name in parameters \
                           if name in kwargs and \
                              'type' in parameters[name]]
        # name in parameters -->  valid inputs in current class
        # name in kwargs       -->  existing inputs in current instance
        # 'type' in parameters[name]   --> type is specified to be checked

        for name, types, value in arg_type_list:
            #(Ex: types = (callable,int)
            if not isinstance(types, (list, tuple)):
                types = [types]  # make a type list
            checked_type = False
            for tp in types:
                if tp == callable and callable(value):
                    # value should be a callable object
                    checked_type = True
                else:
                    if isinstance(value, tp):
                        checked_type = True
            if not checked_type:
                raise TypeError('set: %s is %s, not %s' % \
                                (name, type(value), types))
        return True


    def check_input_range(self,**kwargs):
        """Check whether all existing inputs are in right specified range."""

        parameters = self._parameters
        arg_type_list = [(name,parameters[name]['range'],kwargs[name]) \
                           for name in parameters \
                           if name in kwargs and \
                              'range' in parameters[name]]
        # name in parameters -->  valid inputs in current class
        # name in kwargs       -->  existing inputs in current instance
        # 'range' in parameters[name]   --> range is specified to be checked

        for name, ranges, value in arg_type_list:
            if isinstance(value, (float, int, complex)):
               # value is a comargble number
                if len(ranges) == 2:  # ranges is an interval
                    low, high = ranges
                    if not ((low <= value <= high) or (low >= value >= high)):
                        raise ValueError('%s=%s is illegal - range=[%s, %s]' % \
                                         (name, value, low, high))
                else:    # range is a list of valid values
                    if not value in ranges:
                        raise ValueError('%s=%s is illegal - range=%s' % \
                                         (name, value, str(ranges)))
        return True

    def check_extra(self, **kwargs):
        """
        A parameter may have a keyword ``extra_check`` for user-given
        functions that performs consistency checks on the parameter.
        This method runs the user-given function(s) on the relevant
        set of parameters.
        """
        p = self._parameters
        prm_type_list = [(name, p[name]['extra_check'], kwargs[name]) \
                         for name in p \
                         if name in kwargs and \
                         'extra_check' in p[name]]
        # name in parameters -->  valid inputs in current class
        # name in kwargs       -->  existing inputs in current instance
        # 'extra_check' in parameters[name]
        #           --> extra functions is specified to check the value

        for name, check_funcs, value in prm_type_list:
            try:
                if not check_funcs(value):  # Return false
                    raise ValueError,'''
        Improper value (=%s) for parameter %s.
        Please check your input.''' % (str(value), name)
            except:     # cannot run check_function smoothly
                raise ValueError,'''
        Cannot run check function for %s=%s.
        Please check your input.''' % (name, str(value))
        return True

    def get(self, parameter_name=None, print_info=False):
        """
        Return value of specified input parameters.
        If parameter_name is None, return dict of all inputs.
        """
        if parameter_name is None:
            # Python v2.7 dict comprehension yields shorter code:
            # {name: getattr(self, name) for name in self._parameters}
            all_args = dict([(name, getattr(self, name, None)) \
                                 for name in self._parameters \
                                 if hasattr(self, name)])
            # Remove f and jac since these are wrappers of the
            # user's functions. Instead, insert an entries that
            # reflect the name of user-supplied functions
            all_args.remove('f')
            all_args['name of f'] = self.users_f.func_name
            if 'jac' in all_args:
                all_args.remove('jac')
                all_args['name of jac'] = self.users_jac.func_name

            if print_info:
                print pprint.pformat(all_args)
            return all_args

        else:
            if hasattr(self, parameter_name):
                value = getattr(self, parameter_name)
                if print_info:
                    print "%s = %s" % (parameter_name, value)
                return value
            else:
                raise AttributeError('Parameter %s is not set' % parameter_name)

    def get_parameter_info(self,print_info=False):
        '''
        Return a dictionary containing all properties of all
        legal parameters in current subclass (i.e., the parameters
        in ``self._parameters``).

        If *print_info* is *True*, the ``self._parameters`` dict
        is pretty printed, otherwise it is returned.
        '''
        if print_info:
            print 'Legal parameters for class %s are:' % self.__class__.__name__
            print pprint.pformat(self._parameters)
            return None
        else:
            return self._parameters

    def _print_method(self, with_f, default):
        """
        Return "classname(f=func, param1=..., param2=..., )".
        Skip f= if *with_f* is *False*.
        If *default* is False, skip param1= if the value equals
        the default value.
        (Helper method used by __str__ and __repr__.)
        """
        s = self.__class__.__name__
        args = []
        if with_f and hasattr(self, 'users_f'):
            if not hasattr(self.users_f, '__name__'):     # class instance?
                f_name = self.users_f.__class__.__name__
	    else:    # Ordinary functions
                f_name = self.users_f.__name__
	        if f_name == '<lambda>':   # lambda function
	  	    f_name = 'lambda u, t: ...'
            args.append('f=%s' % f_name)

        # form all parameters
        for name in self._parameters:
            if name != 'f' and hasattr(self, name):
                value = getattr(self, name)
                value_specified = True \
                    if 'default' not in self._parameters[name] \
                    else value != self._parameters[name]['default']
                if default or value_specified:
                    types = self._parameters[name]['type']
                    if types in (callable, (str, callable)):
                        value = getattr(value, '__name__', \
                                            value.__class__.__name__)
                    args.append('%s=%s' % (name, value))

        args = ', '.join(args)
        s += '(%s)' % args
        return s

    def __repr__(self):
        """Return solvername(f=..., param1=..., etc.)."""
        return self._print_method(with_f=True, default=True)

    def __str__(self):
        """
        Return solvername, plus parameters that are different from
        the default value.
        """
        return self._print_method(with_f=False, default=False)


    def set_initial_condition(self, U0):
        """
        Function set_initial_condition() is used to set initial value of
        independent variables.
        """
        # Test first if U0 is sequence (len(U0) possible),
        # and use that as indicator for system of ODEs.
        # The below code should work for U0 having
        # float,int,sympy.mpmath.mpi and other objects as elements.
        try:
            self.neq = len(U0)
            U0 = np.asarray(U0)          # (assume U0 is sequence)
        except TypeError:
            # U0 has no __len__ method, assume as a scalar
            self.neq = 1
            if isinstance(U0, int):
                U0 = float(U0)           # avoid integer division
        self.U0 = U0

    def solve(self, time_points, terminate=None):
        """
        Compute discrete solution u of the ODE problem at time points
        specified in the array time_points. An optional user-supplied
        function ``terminate(u, t, step_no)`` can be supplied to
        terminate the solution process (``terminate`` returns True
        or False) at some time earlier than ``time_points[-1]``.

        Most classes in this solver hierarchy inherit this ``solve``
        method and implement their special ``advance`` method to
        advance the solution one step.
        Some solver classes will implement their own ``solve``
        method, for instance if they wrap some underlying software
        that has a suitable ``solve`` functionality.

        The algorithm steps in this ``solve`` method goes as follows.
        The initialize_for_solve method is called to initialize
        various data needed in the solution process (self. u, for instance).
        Thereafter, ``validate_data`` is called to perform a consistency
        check on data. We are then ready for the core of the method:
        the time loop.

        Output:
           u            : array to hold solution values corresponding to points
           t            : array to hold time values.Usually same as time_points
        """
        if terminate is None:    # Default function
            terminate = lambda u, t, step_no: False

        self.t = np.asarray(time_points)
        self.n = 0  # time step counter
        self.initialize_for_solve()
        self.validate_data()

        # The time loop
        N = self.t.size - 1  # no of intervals
        for n in range(N):
            self.n = n
            self.u[n+1] = self.advance()   # new value
            if terminate(self.u, self.t, n+1):
                print self.__class__.__name__, \
                      'terminated at t=%g' % self.t[n+1]
                self.u, self.t = self.u[:n+2], self.t[:n+2]
                break  # terminate time loop over n
        return self.u, self.t


    def advance(self):
        """Advance solution one time step."""
        raise NotImplementedError


    def initialize_for_solve(self):
        """
        Setting values of internal attributes to be used in iteration.

        These internal attributes are ususally dependent on the values of
        other attributes. For example, for Rkc, self.itol should be
        initialized here as a flag to indicate whether self.atol is
        supplied as scalar or sequence.

        In subclasses, this function can be extended when required.
        """

        # Detect whether data type is in complex type or not.
        if hasattr(self, 'f'):
            value = np.array(self.f(self.U0, self.t[0]))
        else:
            value = np.asarray(self.U0)
        if value.dtype in (np.int32, np.int64, np.int):
            self.dtype = np.complex if self.complex_valued else np.float
        else:
            self.dtype = value.dtype
        self.complex_valued = (str(self.dtype)[:7] == 'complex')

        # Initialization of self.u
        N = self.t.size - 1  # no of intervals
        if self.neq == 1:  # scalar ODEs
            self.u = np.zeros(N+1, self.dtype)
        else:              # systems of ODEs
            self.u = np.zeros((N+1, self.neq), self.dtype)
        # Assume that self.t[0] corresponds to self.U0
        self.u[0] = self.U0

        return None


    def constant_time_step(self):
        """Check if self.t has a uniform partition."""
        return np.allclose(self.t,
                           np.linspace(self.t[0], self.t[-1], len(self.t)))


    def validate_data(self):
        """
        This function is used for extra checking and validating of
        attributes before the computations start.

        This version checks that the ``time_points`` is correctly
        set up. The function also check that all required parameters
        are initialized. Subclass versions may introduce additional
        tests and help.
        """
        # self.t should be a sequence of numbers
        if (not isinstance(self.t, (list, tuple, np.ndarray))) \
            or (not np.asarray(
            # all items in self.t should be numbers
            [isinstance(t, (int,float)) for t in self.t]).all()):
                raise TypeError, \
                    'solve: time_points(=%s) is not a proper '\
                    'sequence of real numbers' % str(self.t)

        # self.t should be supplied in an asscending/descending order
        t_sorted = sorted(self.t, reverse=self.t[0] > self.t[-1])
        if list(self.t) != list(t_sorted):
            raise ValueError, \
                'time_points(=%s) is not provided in an ascending/descending'\
                ' order!' % str(self.t)

        # Test whether all required parameters are provided
        for arg in self._required_parameters:
            if not hasattr(self, arg):
                raise ValueError,\
                    '"%s" has to be input as required parameter(s) for '\
                    'solver %s.' % (arg,self.__class__.__name__)
        return True

    def switch_to(self, solver_target, print_info=False, **kwargs):
        """
        Create a new solver instance which switch to another subclass with
        same values of common attributes.

        `solver_target` is either as a string (class name) or
        a class, i.e., ``RK4`` or ``'RK4'``.
        The `kwargs` arguments are optional parameters to
        reset/supplement values of argmenters in the solver we switch to.
        The instance of the target solver is returned.

        Example:

        >>> import odespy
        >>> f = lambda u,t: -u
        >>> time_points = np.linspace(0.,2.,11)
        >>> exact_u = np.exp(-time_points)
        >>> m1 = odespy.RK2(f)
        >>> m1.set_initial_condition(1.)
        >>> u1, t = m1.solve(time_points)
        >>> print 'Normarized error with RK2 is %g' % np.linalg.norm(u1 - exact_u)
        Normarized error with RK2 is 0.0077317
        >>> m2 = m1.switch_to(odespy.RKFehlberg, rtol=1e-18)
        >>> u2, t = m2.solve(time_points)
        >>> print 'Normarized error with RKFehlberg is %g' % np.linalg.norm(u2 - exact_u)
        Normarized error with RKFehlberg is 8.55517e-08
        """
        # Extract name list of all the subclasses in this module
        solver_list = list_all_solvers()
        error_message = '''
        Input error! Your input %s is not a valid solver name!
        Valid names are %s ''' % (str(solver_target), str(solver_list))

        # Check whether input solver_target is a valid subclass of Solver
        try:
            import odespy
            if type(solver_target) is str:
                # Convert string to class name
                solver_target = getattr(odespy, solver_target)
        except:
            raise ValueError, error_message

        if not solver_target.__name__ in list_all_solvers():
            raise ValueError, error_message



        # Neglect the attributes if they are illegal in target solver
        args_dict = {}
        # Extract all the existing attributes which are legal both in
        # current solver and the target solver
        common_attr = set(solver_target._optional_parameters) & \
                      set(self._optional_parameters)
        # Extract values of these common attributes
        args_dict = dict((name, getattr(self, name)) for name in common_attr \
                             if hasattr(self, name))

        # Exception: 'f' is to provided as 1st parameter to initialize
        # new instance in target solver
        if 'f' in kwargs:    # f is reset as a new parameter
            f = kwargs.pop('f')
        else:
            f = getattr(self, 'f', None)  # the wrapped general form f(u,t)
        for name in ('f_f77', 'f_args', 'f_kwargs', 'f'):
            if name in args_dict:
                del args_dict[name]

        # Union with new values in kwargs
        # Note: Old values for parameters in current solver are neglected
        #       if they are provided in kwargs.
        args_dict.update(kwargs)

        # Create new instance through normal constructor __init__().
        # It ensures all the necessary checking/setting in target solver.
        new = solver_target(f, **args_dict)

        # Set up initial value U0 if available
        if hasattr(self,'U0'):
            new.set_initial_condition(self.U0)

        # Print out information if desired
        if print_info:
            # neglected attributes in new solver
            diff_args = set(self.__dict__.keys()) - set(new.__dict__.keys()) \
                - set(('u','t','n','dtype'))
            if diff_args:
                print 'These attributes are neglected in %s: %s\n' \
                    % (solver_target.__name__, str(diff_args)[5:-2])
            print 'Switched to solver %s' % str(solver_target.__name__)

        return new

    def check_conditional_parameters(self):
        """
        This function is used to check whether conditional parameters are
        provided when specified condition fulfilled.

        This function is not intended for simple solvers.
        So it is not called automatically in current ODE.py.
        But for some complicated solvers as ones in ODEPACK, they
        are very useful and convenient.

        Future developers can apply these functions at appropriate
        locations with corresponding property-setting in
        adjust_parameters().

        For example, in Lsode_ODEPACK, when iter_method is set to 4, it
        indicates that ODEPACK would apply user-supplied banded Jacoabian
        function in corrector iteration. Then we need to confirm either
        'jac_banded' or 'jac_fortran' is supplied. Besides, 'ml' & 'mu' are
        also necessary for iteration with banded Jacobian matrix.
        Thus in order to confirm sufficient conditional inputs, we set
        parameters['iter_method']['condition_list'] =
             {...,'4': (('jac_banded','jac_fortran'),ml,mu),...}

        In this function, we would check all the legal parameters with specified
        condition-list, and make sure all the conditional parameters with
        current value is supplied.

        """
        parameters = self._parameters
        # Parameters with condition-list settings
        with_condition_args = \
            [(name, parameters[name]['condition-list'], \
              str(getattr(self,name))) \
                 for name in parameters \
                    if name in self.__dict__ and
                       'condition-list' in parameters[name]]

        # name in parameters   -->  valid inputs for current class
        # name in self.__dict__ -->  existing inputs for curremt instance
        # 'condition-list' in parameters[name] -->
        #                       'condition-list' is specified to check

        for (name, conditions, value) in with_condition_args:
            # Ex: name = 'iter_method'
            #     conditions = {'1':(('jac','jac_f77'),), '4':.., '5':..})
            #     value = '1'
            if value in conditions:
                # There is conditional requirements for current value
                condition_args = conditions[value]
                # list/tuple for conditional parameters
                for arg in condition_args:
                    if not isinstance(arg, str):
                        # arg is a list for alternative parameters
                        # e.g. ('jac', 'jac_f77')
                        # Either 'jac' or 'jac_f77' should be supplied
                        found = bool([p_name for p_name in arg \
                                          if hasattr(self,p_name)])
                        arg_print = 'One of %s' % str(arg)
                    else:      # arg is a single parameter
                        found = hasattr(self, arg)
                        arg_print = arg
                    if not found:
                        raise ValueError,'''\
        Error! Unsufficient input!
        %s must be set when %s is %s!''' % (arg_print,name,value)
        return True


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
        "u,t,j" to "t,u,j+1". That is, define the jacobian function as
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
                      parameters[name].get('returnArrayOrder', None),
                      parameters[name].get('paralist_old', None),
                      parameters[name].get('paralist_new', None),
                      parameters[name].get('name_wrapped', name)]
                     for name in parameters \
                         if name in self.__dict__ and \
                         'type' in parameters[name] and \
                         (parameters[name]['type'] is callable or \
                          parameters[name]['type'] is (callable, str)) and \
                         ('paralist_new' in parameters[name] or \
                          'returnArrayOrder' in parameters[name])]
        # name in self.__dict__  --> existing attributes in current instance
        # parameters[name]['type'] is callable or (callable, str)
        #             -->     callable objects
        # 'paralist_new' in parameters[name]
        #    --> new parameter-list is defined to be wrapped
        # 'returnArrayOrder' in parameters[name]
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


class MySolver(Solver):
    """
    Users can define a solver with supplying a function
    myadvance(), and make use of all possible parameters
    in this module::

        myadvance(MySolver_instance)  -->  return u_new

    Example::

        def myadvance_(ms):
            f, u, t, n, atol = ms.f, ms.u, ms.t, ms.n, ms.atol
            # All class attributes can be obtained
            u_new = ...
            return u_new

        def f(u,t):
            udot = ...
            return udot

        method = MySolver(f, myadvance=myadvance_)
        method.set_initial_condition(u0)
        u,t = method.solve(time_points)
    """
    _required_parameters = ['f', 'myadvance']
    _optional_parameters = _parameters.keys()
    # All arguments are valid and accessible for users.

    def advance(self):
        return self.myadvance(self)

### End of class Solver ###

class ForwardEuler(Solver):
    """
    Forward Euler scheme::

        u[n+1] = u[n] + dt*f(u[n], t[n])
    """
    quick_description = 'The simple explicit (forward) Euler scheme'

    def advance(self):
        u, f, n, t = self.u, self.f, self.n, self.t
        dt = t[n+1] - t[n]
        u_new = u[n] + dt*f(u[n], t[n])
        return u_new

Euler = ForwardEuler   # synonym


class Leapfrog(Solver):
    """
    Leapfrog scheme::

        u[n+1] = u[n-1] + dt2*f(u[n], t[n])

    with::

        dt2 = t[n+1] - t[n-1]

    Forward Euler is used for the first step.
    """
    quick_description = 'Standard explicit Leapfrog scheme'

    def advance(self):
        u, f, n, t = self.u, self.f, self.n, self.t

        if n >= 1:
            dt2 = t[n+1] - t[n-1]
            u_new = u[n-1] + dt2*f(u[n], t[n])
        else:
            dt = t[n+1] - t[n]
            u_new = u[n] + dt*f(u[n], t[n])
        return u_new



class LeapfrogFiltered(Solver):
    """
    The standard Leapfrog scheme reads::

        u[n+1] = u[k-1] + dt2*f(u[n], t[n])

    with::

        dt2 = t[n+1] - t[k-1]

    Forward Euler is used for the first step.
    Since Leapfrog gives oscillatory solutions, this class
    applies a common filtering technique::

        u[n] = u[n] + gamma*(u[n-1] - 2*u[n] + u[n+1])

    with gamma=0.6 as in the NCAR Climate Model.
    """
    quick_description = 'Filtered Leapfrog scheme'

    def advance(self):
        u, f, n, t = self.u, self.f, self.n, self.t
        gamma = 0.6  # NCAR Climate Model

        if n >= 1:
            dt2 = t[n+1] - t[n-1]
            u_new = u[n-1] + dt2*f(u[n], t[n])
            u[n] = u[n] + gamma*(u[n-1] - 2*u[n] + u_new)
        else:
            dt = t[n+1] - t[n]
            u_new = u[n] + dt*f(u[n], t[n])
        return u_new



class Heun(Solver):
    """
    Heun's method, also known as an RungeKutta2 or Trapezoidal method.
    Basically, it is a central difference method, with one
    iteration and the Forward Euler scheme as start value.
    In this sense, it is a predictor-corrector method.

    Scheme::

        u[n+1] = u[n] + 0.5*dt*(f(u[n],t[n]) + f(u[n]+dt*f(u[n],t[n]),t[n+1]))
    """
    quick_description = "Heun's explicit method (similar to RK2)"

    def advance(self):
        u, f, n, t = self.u, self.f, self.n, self.t
        dt = t[n+1] - t[n]
        u_star = u[n] + dt*f(u[n], t[n])  # Forward Euler step
        u_new = u[n] + 0.5*dt*(f(u[n], t[n]) + f(u_star, t[n+1]))
        return u_new



Trapezoidal = Heun   # alias for solver Heun


class RK2(Solver):
    """
    Standard Runge-Kutta 2nd method::

        u[n+1] = u[n] + dt*f(u[n] + 0.5*(dt*f(u[n],t[n])),t[n] + 0.5*dt)
    """
    quick_description = "Explicit 2nd-order Runge-Kutta method"

    def advance(self):
        u, f, n, t = self.u, self.f, self.n, self.t
        dt = t[n+1] - t[n]
        K1 = dt*f(u[n], t[n])
        K2 = dt*f(u[n] + 0.5*K1, t[n] + 0.5*dt)
	u_new = u[n] + K2
        return u_new


class RK4(Solver):
    """
    Standard RK4 method::

        u[n+1] = u[n] + (1/6.0)*(K1 + 2*K2 + 2*K3 + K4)

    where::

           K1 = dt*f(u[n], t[n])
           K2 = dt*f(u[n] + 0.5*K1, t[n] + 0.5*dt)
           K3 = dt*f(u[n] + 0.5*K2, t[n] + 0.5*dt)
           K4 = dt*f(u[n] + K3, t[n] + dt)
    """
    quick_description = "Explicit 4th-order Runge-Kutta method"

    def advance(self):
        u, f, n, t = self.u, self.f, self.n, self.t
        dt = t[n+1] - t[n]
        dt2 = dt/2.0
        K1 = dt*f(u[n], t[n])
        K2 = dt*f(u[n] + 0.5*K1, t[n] + dt2)
        K3 = dt*f(u[n] + 0.5*K2, t[n] + dt2)
        K4 = dt*f(u[n] + K3, t[n] + dt)
        u_new = u[n] + (1/6.0)*(K1 + 2*K2 + 2*K3 + K4)
        return u_new


class RK3(Solver):
    """
    RungeKutta3 method::

        u[n+1] = u[n] + (1/6.0)*(K1 + 4*K2 + K3)

    where::

        K1 = dt*f(u[n], t[n])
        K2 = dt*f(u[n] + 0.5*K1, t[n] + 0.5*dt)
        K3 = dt*f(u[n] - K1 + 2*K2, t[n] + dt)
    """
    quick_description = "Explicit 3rd-order Runge-Kutta method"

    def advance(self):
        u, f, n, t = self.u, self.f, self.n, self.t
        dt = t[n+1] - t[n]
        dt2 = dt/2.0
        K1 = dt*f(u[n], t[n])
        K2 = dt*f(u[n] + 0.5*K1, t[n] + dt2)
        K3 = dt*f(u[n] - K1 + 2*K2, t[n] + dt)
        u_new = u[n] + (1/6.0)*(K1 + 4*K2 + K3)
        return u_new


class AdamsBashforth2(Solver):
    """
    Second-order Adams-Bashforth method::

        u[n+1] = u[n] + dt/2.*(3*f(u[n], t[n]) - f(u[n-1], t[n-1]))

    for constant time step dt.

    RK2 is used as default solver in first step.
    """
    quick_description = "Explicit 2nd-order Adams-Bashforth method"

    _optional_parameters = Solver._optional_parameters + ['start_method',]

    def initialize_for_solve(self):
        # New solver instance for first steps
        self.starter = self.switch_to(self.start_method)
        # Create variables for holding f at previous time levels
        self.f_n_1 = None
        Solver.initialize_for_solve(self)

    def validate_data(self):
        """Check that the time steps are constant."""
        if not self.constant_time_step():
            print '%s must have constant time step' % self.__name__
            return False
        else:
            return True

    def advance(self):
        u, f, n, t = self.u, self.f, self.n, self.t

        if n >= 1:
            dt = t[n+1] - t[n]  # must be constant
            self.f_n = f(u[n], t[n])
            u_new = u[n] + dt/2.*(3*self.f_n - self.f_n_1)
            self.f_n_1 = self.f_n
        else:
            # User-specified method for the first step
            self.starter.set_initial_condition(u[n])
            time_points = [t[n], t[n+1]]
            u_starter, t_starter = self.starter.solve(time_points)
            u_new = u_starter[-1]
            self.f_n_1 = f(u[0], t[0])

        return u_new


class AdamsBashforth3(Solver):
    """
    Third-order Adams-Bashforth method::

        u[n+1] = u[n] + dt/12.*(23*f(u[n], t[n]) - 16*f(u[n-1], t[n-1])
                                + 5*f(u[n-2], t[n-2]))

    for constant time step dt.

    RK2 is used as default solver for first steps.
    """
    quick_description = "Explicit 3rd-order Adams-Bashforth method"

    _optional_parameters = Solver._optional_parameters + ['start_method',]

    def initialize_for_solve(self):
        # New solver instance for first steps
        self.starter = self.switch_to(self.start_method)
        # Create variables for f at previous time levels
        self.f_n_1 = None
        self.f_n_2 = None
        Solver.initialize_for_solve(self)

    def validate_data(self):
        if not self.constant_time_step():
            print '%s must have constant time step' % self.__name__
            return False
        else:
            return True

    def advance(self):
        u, f, n, t = self.u, self.f, self.n, self.t

        if n >= 2:
            dt = t[n+1] - t[n]  # must be constant
            self.f_n = f(u[n], t[n])
            u_new = u[n] + dt/12.*(23*self.f_n - 16*self.f_n_1 + 5*self.f_n_2)
            self.f_n_1, self.f_n_2, self.f_n = self.f_n, self.f_n_1, self.f_n_2

        else:
            # Start method
            self.starter.set_initial_condition(u[n])
            time_points = [t[n], t[n+1]]
            u_starter, t_starter = self.starter.solve(time_points)
            u_new = u_starter[-1]
            if n == 0:
                self.f_n_2 = f(u[0], t[0])
            elif n == 1:
                self.f_n_1 = f(u[1], t[1])
        return u_new


class AdamsBashMoulton2(Solver):
    """
    Two-step (3rd-order) Adams-Bashforth method::

        predictor = u[n] + dt/12.*(23.*f(u[n], t[n]) - 16*f(u[n-1], t[n-1]) +
                            5*f(u[n-2], t[n-2]))
        corrector = u[n] + dt/12.*(8.*f(u[n], t[n]) - f(u[n-1], t[n-1]) +
                            5*f(predictor, t[n+1]))

    for constant time step dt.

    RK2 is used as default solver for first steps.
    """
    quick_description = "Explicit 2nd-order Adams-Bashforth-Moulton method"

    _optional_parameters = Solver._optional_parameters + ['start_method',]

    def initialize_for_solve(self):
        # New solver instance for first steps
        self.starter = self.switch_to(self.start_method)
        # Create variables for f at previous time levels
        self.f_n_1, self.f_n_2 = None, None
        Solver.initialize_for_solve(self)

    def validate_data(self):
        if not self.constant_time_step():
            print '%s must have constant time step' % self.__name__
            return False
        else:
            return True

    def advance(self):
        u, f, n, t = self.u, self.f, self.n, self.t

        if n >= 2:
            dt = t[n+1] - t[n]  # must be constant
            self.f_n = f(u[n], t[n])
            predictor = u[n] + dt/12.*(23.*self.f_n - 16*self.f_n_1 + \
                                  5*self.f_n_2)
            u_new = u[n] + dt/12.*(8*self.f_n - self.f_n_1 + \
                                  5*f(predictor, t[n + 1]))
            self.f_n_1, self.f_n_2 = self.f_n, self.f_n_1
        else:
            # Start method
            self.starter.set_initial_condition(u[n])
            time_points = [t[n], t[n+1]]
            u_starter, t_starter = self.starter.solve(time_points)
            u_new = u_starter[-1]
            if n == 0:
                self.f_n_2 = f(u[0], t[0])
            elif n == 1:
                self.f_n_1 = f(u[1], t[1])

        return u_new


class AdamsBashforth4(Solver):
    """
    Fourth-order Adams-Bashforth method::

        u[n+1] = u[n] + dt/24.*(55.*f(u[n], t[n]) - 59*f(u[n-1], t[n-1]) +
                                37*f(u[n-2], t[n-2]) - 9*f(u[n-3], t[n-3]))

    for constant time step dt.

    RK2 is used as default solver for first steps.
    """
    quick_description = "Explicit 4th-order Adams-Bashforth method"

    _optional_parameters = Solver._optional_parameters + ['start_method',]

    def initialize_for_solve(self):
        # New solver instance for first steps
        self.starter = self.switch_to(self.start_method)
        # Create variables for f at previous time levels
        self.f_n_1, self.f_n_2, self.f_n_3 = None, None, None
        Solver.initialize_for_solve(self)

    def validate_data(self):
        if not self.constant_time_step():
            print '%s must have constant time step' % self.__name__
            return False
        else:
            return True

    def advance(self):
        u, f, n, t = self.u, self.f, self.n, self.t

        if n >= 3:
            dt = t[n+1] - t[n]  # must be constant
            self.f_n = f(u[n], t[n])
            u_new = u[n] + dt/24.*(55.*self.f_n - 59*self.f_n_1 + \
                                  37*self.f_n_2 - 9*self.f_n_3)
            self.f_n_1, self.f_n_2, self.f_n_3 = \
                self.f_n, self.f_n_1, self.f_n_2
        else:
            # Start method
            self.starter.set_initial_condition(u[n])
            time_points = [t[n], t[n+1]]
            u_starter, t_starter = self.starter.solve(time_points)
            u_new = u_starter[-1]
            if n == 0:
                self.f_n_3 = f(u[0], t[0])
            elif n == 1:
                self.f_n_2 = f(u[1], t[1])
            elif n == 2:
                self.f_n_1 = f(u[2], t[2])

        return u_new


class AdamsBashMoulton3(Solver):
    """
    Three-step (4th-order) Adams-Bashforth method::

        predictor = u[n] + dt/24.*(55.*f(u[n], t[n]) - 59*f(u[n-1], t[n-1]) +
                                   37*f(u[n-2], t[n-2]) - 9*f(u[n-3], t[n-3]))
        corrector = u[n] + dt/24.*(19.*f(u[n], t[n]) - 5*f(u[n-1], t[n-1]) +
                                   f(u[n-2], t[n-2]) + 9*f(predictor, t[n+1]))

    for constant time step dt.

    RK2 is used as default solver for first steps.
    """
    quick_description = "Explicit 3rd-order Adams-Bashforth-Moulton method"

    _optional_parameters = Solver._optional_parameters + ['start_method',]

    def initialize_for_solve(self):
        # New solver instance for first steps
        self.starter = self.switch_to(self.start_method)
        # Create variables for f at previous time levels
        self.f_n_1, self.f_n_2, self.f_n_3 = None, None, None
        Solver.initialize_for_solve(self)

    def validate_data(self):
        if not self.constant_time_step():
            print '%s must have constant time step' % self.__name__
            return False
        else:
            return True

    def advance(self):
        u, f, n, t = self.u, self.f, self.n, self.t

        if n >= 3:
            dt = t[n+1] - t[n]  # must be constant
            self.f_n = f(u[n], t[n])
            predictor = u[n] + dt/24.*(55.*self.f_n - 59*self.f_n_1 + \
                                  37*self.f_n_2 - 9*self.f_n_3)
            u_new = u[n] + dt/24.*(self.f_n_2 - 5*self.f_n_1 + 19*self.f_n + \
                                  9*f(predictor, t[n + 1]))
            self.f_n_1, self.f_n_2, self.f_n_3 = \
                self.f_n, self.f_n_1, self.f_n_2
        else:
            # Start method
            self.starter.set_initial_condition(u[n])
            time_points = [t[n], t[n+1]]
            u_starter, t_starter = self.starter.solve(time_points)
            u_new = u_starter[-1]
            if n == 0:
                self.f_n_3 = f(u[0], t[0])
            elif n == 1:
                self.f_n_2 = f(u[1], t[1])
            elif n == 2:
                self.f_n_1 = f(u[2], t[2])

        return u_new



class MidpointIter(Solver):
    """
    A midpoint/central difference method with max_iter fixed-point
    iterations to solve the nonlinear system.
    The Forward Euler scheme is recovered if max_iter=1 and f(u,t)
    is independent of t. For max_iter=2 we have the Heun/RK2 scheme.
    """
    quick_description = "Explicit 2nd-order iterated Midpoint method"

    _optional_parameters = Solver._optional_parameters + \
                           ['max_iter', 'eps_iter']

    def adjust_parameters(self):
        self._parameters['max_iter']['default'] = 3

    def advance(self):
        if not hasattr(self, 'v'):  # v is a help array needed in the method
            if self.neq == 1:
                # Scalar ODE: v can be one-dim array
                self.v = np.zeros(self.max_iter+1, self.u.dtype)
            else:
                # System of ODEs: v must be two-dim array
                self.v = np.zeros((self.max_iter+1, self.neq), self.u.dtype)

        u, f, n, t, v = \
           self.u, self.f, self.n, self.t, self.v
        dt = t[n+1] - t[n]

        v[0] = u[n]
        q = 0
        v_finished = False   # |v[q]-v[q-1]| < eps
        while not v_finished and q < self.max_iter:
            q += 1
            v[q] = u[n] + 0.5*dt*(f(v[q-1], t[n+1]) + f(u[n], t[n]))
            if abs(v[q] - v[q-1]).max() < self.eps_iter:
                v_finished = True
                self.num_iterations = q

        u_new = v[q]
        return u_new


def approximate_Jacobian(f, u0, t0, h):
    """
    Compute approximate Jacobian of fucntion f at current point (u0,t0).
    Method: forward finite difference approximation with step
    size h.
    Output: a two-dimensional array holding the Jacobian matrix.
    """
    u0 = np.asarray(u0)
    f0 = np.asarray(f(u0, t0))
    neq = len(u0)
    J = np.zeros((neq, neq), float)
    for i in range(neq):
        u_new = u0.copy()
        u_new[i] += h
        J[i] = (np.asarray(f(u_new, t0)) - f0)/h
    return J.transpose()


class SymPy_odefun(Solver):
    """
    Wrapper for the sympy.mpmath.odefun method, which applies a high-order
    Taylor series method to solve ODEs.
    """
    quick_description = "Very accurate high order Taylor method (from SymPy)"

    def initialize(self):
        try:
            import sympy
            self.sympy = sympy
        except ImportError:
            raise ImportError,'sympy is not installed - needed for sympy_odefun'

    def initialize_for_solve(self):
        # sympy.odefun requires f(t, u), not f(u, t, *args, **kwargs)
        self.f4odefun = lambda t, u: self.f(u, t, *self.f_args, **self.f_kwargs)
        Solver.initialize_for_solve(self)

    def solve(self, time_points, terminate=None):
        """
        The complete solve method must be overridded in this class
        since sympy.mpmath.odefun is such a solve method.

        The class stores an attribute ufunc (return from odefun)
        which can be used to evaluate u at any time point (ufunc(t)).
        """
        if terminate is not None:
            print 'Warning: SymPy_odefun.solve ignores the terminate function!'
        self.t = np.asarray(time_points)
        self.initialize_for_solve()

        self.sympy.mpmath.mp.dps = 15  # accuracy
        self.ufunc = self.sympy.mpmath.odefun(
            self.f4odefun, time_points[0], self.U0)

        # u and t to be returned are now computed by sampling self.ufunc
        # at the specified time points
        self.u = np.array([self.ufunc(t) for t in time_points])
        self.t = np.asarray(time_points)
        return self.u, self.t


class odelab(Adaptive):
    """
    Wrapper for the odelab package.
    """
    quick_description = "interface to all solvers in odelab"

    _required_parameters = Adaptive._required_parameters + \
        ['odelab_solver']
    _optional_parameters = Adaptive._optional_parameters + \
        ['jac', 'jac_args', 'jac_kwargs', ]

    solvers = 'ExplicitEuler ImplicitEuler RungeKutta4 ExplicitTrapezoidal RungeKutta34 SymplecticEuler Heun Kutta AdamsBashforth Butcher LDIRK343 LobattoIIIA LobattoIIIB LobattoIIIC LobattoIIICs LobattoIIID RadauIIA'.split()

    def initialize(self):
        try:
            import odelab
            self.odelab = odelab
        except ImportError:
            raise ImportError,'odelab is not installed - needed for sympy_odefun'

    def initialize_for_solve(self):
        # odelab requires f(t, u), not f(u, t, *args, **kwargs)
        # self.f4odelab(t, u) is what we pass on to odelab
        self.f4odelab = lambda t, u: self.f(u, t, *self.f_args, **self.f_kwargs)
        Solver.initialize_for_solve(self)

        if self.odelab_solver not in odelab.solvers:
            raise ValueError('requested solver %s not in %s' % \
                             (self.odelab_solver, str(odelab.solvers)))

        self.system = self.odelab.System(self.f4odelab)

        h = self.t[1] - self.t[0]
        # odelab solvers are in different modules...
        for module in [self.odelab.scheme.classic,
                       self.odelab.scheme.geometric,
                       self.odelab.scheme.rungekutta,
                       ]:
            if self.odelab_solver in dir(module):
                self.scheme = getattr(module, self.odelab_solver)(h)
                break
        self.solver = self.odelab.Solver(scheme=self.scheme,
                                         system=self.system)

        # Set initial condition
        self.solver.initialize(self.U0)
        # Solve problem in odelab
        self.solver.run(self.t[-1])
        # Retrive solution
        with self.solver.open_store() as events:
            for i in range(self.neq):
                self.u[:, i] = events[i]
        # PROBLEM: will the events correspond to self.t? Probably
        # not for adaptive methods...
        return self.u, self.t


    def solve(self, time_points, terminate=None):
        """
        The complete solve method must be overridded in this class
        since odelab. is such a solve method.
        """
        if terminate is not None:
            print 'Warning: odelab.solve ignores the terminate function!'
        self.t = np.asarray(time_points)
        self.initialize_for_solve()

        scheme = self.odelab.scheme.classic.ExplicitEuler
        self.sympy.mpmath.mp.dps = 15  # accuracy
        self.ufunc = self.sympy.mpmath.odefun(
            self.f4odefun, time_points[0], self.U0)

        # u and t to be returned are now computed by sampling self.ufunc
        # at the specified time points
        self.u = np.array([self.ufunc(t) for t in time_points])
        self.t = np.asarray(time_points)
        return self.u, self.t


class SolverImplicit(Solver):
    """
    Super class for implicit methods for ODEs.
    Existing solvers are: BackwardEuler, Backward2Step, ThetaRule
    """

    _optional_parameters = Solver._optional_parameters + \
        ['jac', 'jac_args', 'jac_kwargs', 'h_in_fd_jac',
         'nonlinear_solver', 'max_iter', 'eps_iter', 'relaxation']

    def initialize_for_solve(self):
        # Set appropriate value of nonlinear_solver if undefined
        if getattr(self, 'jac', None) is None:   # no jac provided
            if getattr(self, 'nonlinear_solver', None) is None:
                self.nonlinear_solver = 'Picard'  # default if no jac provided
            elif getattr(self, 'nonlinear_solver') == 'Newton':
                 # Approximate jacobian with finite difference approx
                self.users_jac = approximate_Jacobian
                self.jac = lambda u, t: \
                    self.users_jac(self.f, u, t, self.h_in_fd_jac)
        else:
            if getattr(self, 'nonlinear_solver', None) is None:
                self.nonlinear_solver = 'Newton'  # default if jac provided
            # Wrap user-supplied Jacobian in the way f is wrapped
            self.users_jac = self.jac
            jac = self.jac  # XXX will not just self.jac work below?
            self.jac = lambda u, t: \
                np.asarray(jac(u, t, *self.jac_args, **self.jac_kwargs))

        Solver.initialize_for_solve(self)

    def advance(self):
        n = self.n
        un, f, t_new, tn = self.u[n], self.f, self.t[n+1], self.t[n]
        dt = t_new - tn

        # General solver routine with Newton or Picard
        # Newton with Finite-Differential or exact Jac
        i, error = 1, 1E+30
        # Forward Euler step for initial guess for nonlinear solver
        u_new = un + (t_new-tn)*f(un,tn)
        # control by number of intern steps and error tolerance
        while i <= self.max_iter and error > self.eps_iter:
            if self.nonlinear_solver == 'Picard':
                u_new = self.Picard_update(u_new)
            elif self.nonlinear_solver == 'Newton':
                F, Jac = self.Newton_system(u_new)
                du = F/Jac if self.neq == 1 else np.linalg.solve(Jac, F)
                u_new = u_new - du
            elif self.nonlinear_solver == 'FixedPoint':
                u_new = un + dt*f(u_new,tn)
            error = np.abs(u_new - u_new).max()
            r = self.relaxation    # Relax factor
            u_new = r*u_new + (1-r)*un
            i += 1
        return u_new

class BackwardEuler(SolverImplicit):
    """
    Implicit Backward Euler method::

       u[n+1] = u[n] + dt*f(t[n+1], u[n+1])

    The nonlinear system is solved by Newton or Picard iteration.
    """
    quick_description = "Implicit 1st-order Backward Euler method"

    def Picard_update(self, ukp1):
        u, f, n, t = self.u, self.f, self.n, self.t
        dt = t[n+1] - t[n]
        return u[n] + dt*f(ukp1, t[n+1])

    def Newton_system(self, ukp1):
        u, f, n, t = self.u, self.f, self.n, self.t
        dt = t[n+1] - t[n]
        F = ukp1 - (u[n] + dt*f(ukp1, t[n+1]))
        J = np.eye(self.neq) - dt*self.jac(ukp1, t[n+1])
        return F, J


class Backward2Step(SolverImplicit):
    """
    Implicit Backward Euler method with 2 steps::

         u[n+1] = u[n]*4/3 - u[n-1]/3 + (t[n+1-t[n-1]])*f(t[n+1], u[n+1])/3

    The 1st-order Backward Euler method is used for the first step.
    """
    quick_description = "Implicit 2nd-order Backward Euler method"

    def Picard_update(self, ukp1):
        u, f, n, t = self.u, self.f, self.n, self.t
        if n == 0:
            # Backward Euler as starter
            dt = t[n+1] - t[n]
            return u[n] + dt*f(ukp1, t[n+1])
        else:
            dt2 = t[n+1] - t[n-1]
            return 4./3*u[n] - 1./3*u[n-1] + (1./3)*dt2*f(ukp1, t[n+1])

    def Newton_system(self, ukp1):
        u, f, n, t = self.u, self.f, self.n, self.t
        if n == 0:
            # Backward Euler as starter
            dt = t[n+1] - t[n]
            F = ukp1 - (u[n] + dt*f(ukp1, t[n+1]))
            J = np.eye(self.neq) - dt*self.jac(ukp1, t[n+1])
        else:
            dt2 = t[n+1] - t[n-1]
            F = ukp1 - (4./3*u[n] - 1./3*u[n-1] + (1./3)*dt2*f(ukp1, t[n+1]))
            J = np.eye(self.neq) - (1./3)*dt2*self.jac(ukp1, t[n+1])
        return F, J


class ThetaRule(SolverImplicit):
    """
    Theta rule method as a typical weighted method with factor theta::

       u[n+1] = u[n] + dt*(theta*f(u[n+1],t[n+1]) + (1 - theta)*f(u[n],t[n]))

    where theta is a float in [0,1].

    The nonlinear system is solved by Picard or Newton iteration.
    """
    quick_description = "Unified Forward/Backward Euler and Midpoint methods"

    _optional_parameters = SolverImplicit._optional_parameters + ['theta']

    def Picard_update(self, ukp1):
        u, f, n, t, theta = self.u, self.f, self.n, self.t, self.theta
        dt = t[n+1] - t[n]
        return u[n] + theta*dt*f(ukp1, t[n+1]) + (1-theta)*dt*f(u[n], t[n])

    def Newton_system(self, ukp1):
        u, f, n, t, theta = self.u, self.f, self.n, self.t, self.theta
        dt = t[n+1] - t[n]
        F = ukp1 - (u[n] + theta*dt*f(ukp1, t[n+1]) + (1-theta)*dt*f(u[n],t[n]))
        J = np.eye(self.neq) - theta*dt*self.jac(ukp1, t[n+1])
        return F, J


class MidpointImplicit(SolverImplicit):
    '''
    Midpoint Implicit method::

       u[n+1] = u[n] + dt*f((u[n+1] + u[n])/2., t[n] + dt/2.)

    The nonlinear system is solved by Picard or Newton iteration.
    '''
    quick_description = "Implicit 2nd-order Midpoint method"

    def Picard_update(self, ukp1):
        u, f, n, t = self.u, self.f, self.n, self.t
        dt = t[n+1] - t[n]
        return u[n] + dt*f((ukp1 + u[n])/2., t[n] + dt/2.)

    def Newton_system(self, ukp1):
        u, f, n, t = self.u, self.f, self.n, self.t
        dt = t[n+1] - t[n]
        F = ukp1 -  u[n] + dt*f((ukp1 + u[n])/2., t[n] + dt/2.)
        J = np.eye(self.neq) - dt*self.jac(ukp1, t[n+1])
        return F, J



class Adaptive(Solver):
    """Superclass for adaptive solvers."""

    _optional_parameters = Solver._optional_parameters + \
        ['rtol', 'atol', 'first_step', 'min_step', 'max_step']

    def initialize_for_solve(self):
        # Let first_step, min_step and max_ste, if not given, be
        # computed from the user's time points, available as self.t.

        if not hasattr(self, 'first_step'):
            # Use user's first step
            self.set(first_step=self.t[1] - self.t[0])
        time_steps = self.t[1:] - self.t[:-1]
        if not hasattr(self, 'min_step'):
            # Use 1/10 of the user's smallest step
            self.set(min_step=0.1*time_steps.min())
        if not hasattr(self, 'max_step'):
            # Use 10 times the user's greatest step
            self.set(max_step=10*time_steps.max())
        Solver.initialize_for_solve(self)


class AdaptiveResidual(Adaptive):
    """
    Designed for educational purposes to demonstrate a possible
    adaptive strategy.

    Currently, only scalar ODE problems can be applied for purpose of
    simplification.
    """
    quick_description = "Very simple adaptive strategy based on the residual"

    _optional_parameters = Adaptive._optional_parameters + \
                           ['solver']

    def __init__(self, f, **kwargs):
        Adaptive.__init__(self, f, **kwargs)
        if 'solver' in kwargs:
            del kwargs['solver']
        self.solver = self.switch_to(self.solver)

    def residual(self, u_n, u_next, t_n, t_next):
        dt = t_next - t_n
        t_mean = 0.5*(t_n + t_next)
        u_mean = 0.5*(u_n + u_next)
        u_diff = u_next - u_n
        # Central 2nd order difference approx to the residual
        # Valid for scalar ODE only
        return abs(u_diff/dt - self.f(u_mean, t_mean))

    def solve(self, time_points, terminate=None, print_info=False):
        self.users_time_points = np.asarray(time_points).copy()
        self.t = t = self.users_time_points
        self.initialize_for_solve()
        # Assume scalar equation...
        if not isinstance(self.U0, (float,int)):
            raise TypeError('Initial condition is not scalar - '
                            'AdaptiveResidual can only work with scalar ODEs')
        self.validate_data()
        self.u = [self.U0]
        self.t = [t[0]]
        self.solver.set_initial_condition(self.U0)
        for k in range(1,len(t)):
            R = 1E+20
            # Try to jump until next user point in time
            ntpoints = 1
            # Halve the time step until residual is small enough
            while R > self.atol:
                ntpoints *= 2
                time_points = np.linspace(t[k-1], t[k], ntpoints)
                dt = time_points[1] - time_points[0]
                if dt < self.min_step:
                    print 'AdaptiveResidual with %s solver, too small %s < %s, R=%s > %s' % (self.solver.__class__.__name__, dt, self.min_step, R, self.atol)
                    print
                    break
                if dt > self.max_step:
                    print 'AdaptiveResidual with %s solver, too large step %s > %s, R=%s > %s' % (self.solver.__class__.__name__, dt, self.max_step, R, self.atol)
                    break

                self.solver.set_initial_condition(self.u[-1])
                u_new, tnew = self.solver.solve(time_points, terminate)
                R = self.residual(u_new[-2], u_new[-1], tnew[-2], tnew[-1])
                if print_info:
                    print '\n%d time points in (t[%d], t[%d]) = (%.3g, %.3g)' \
                        % (ntpoints, k-1, k, t[k-1], t[k])
                    print 'Residual = %g, Tolerance = %g, calling %s' % \
                        (R, self.atol, self.solver.__class__.__name__)
                # reintegrate with time_step = dt/2
            self.u.extend(u_new[1:])
            self.t.extend(tnew[1:])
        return self.u, self.t


class RK34(Adaptive):
    """
    Adaptive 4th-order Runge-Kutta method.
    For each time level t[n], the method takes many
    adaptively chosen (small) time steps to hit the
    next target level t[n+1].
    All computed u and t values are available as
    self.u_adaptive and self.t_adaptive, if desired.
    """
    quick_description = "Adaptive 4th-order Runge-Kutta method"

    def initialize_for_solve(self):
        Adaptive.initialize_for_solve(self)
        self.order = 4
        self.t_all = [self.t[0]]  # computed time levels
        self.u_all = [self.u[0]]  # corresponding u values

    def advance_intermediate(self, u, t, dt):
        """
        Advance the solution an intermediate addaptive step.
        Parmeters: u at time t, dt is current step size.
        """

        f = self.f
        K1 = f(u,                   t)
        K2 = f(u + dt*K1/2.,        t + dt/2.)
        K3 = f(u + dt*K2/2,         t + dt/2)
        Z3 = f(u - dt*K1 + 2*dt*K2, t + dt)
        K4 = f(u + dt*K3,           t + dt)

        u_new = dt/6.*(K1 + 2*K2 + 2*K3 + K4)

        error = dt/6.*(2*K2 + Z3 - 2*Y3 - Y4)
        error = np.linalg.norm(error)  # scalar measure

        return u_new, error

    def adjust_timestep(self, dt, error, u_new, order):
        """
        Adjust dt, given error and order::

            tol = self.rtol*np.linalg.norm(u_new) + self.atol

            dt_new = dt*(tol/error)**(1./order)
        """
        tol = self.rtol*np.linalg.norm(u) + self.atol
        if error > 1E-14:
            dt_new = dt*(tol/error)**(1./order)
            if self.min_step <= dt_new <= self.max_step:
                dt = dt_new
        return dt


    def advance(self):
        """
        Advance from t[n] to t[n+1] in (small) adaptive steps.

        Adaptivity logic: start with dt=t[n+1]-t[n] (one big
        step), repeat computing the solution and an error
        estimate e until e < tol.
        """
        # u, t: solution u at time t in between self.t[n] and self.t[n+1]
        # dt: adaptive time step in between self.t[n] and self.t[n+1]
        n = self.n
        u, t, t_np1 = self.u[n], self.t[n], self.t[n+1]
        dt = tnp1 -t

        min_step = getattr(self, 'min_step', dt/1000.)
        max_step = getattr(self, 'max_step', dt)
        if not hasattr(self, 'h'):
            # h is current adaptive step size, use first_step if
            # h is not yet computed
            first_step = getattr(self, 'first_step', dt)
            if first_step > dt:
                first_step = dt
            self.h = first_step
        # else: start out with previous self.h step size

        if self.verbose >= 1:
            print '\nat user-specified time level %g, ',
            print 'starting with h=%g' % (t, self.h)

        while t <= t_np1:
            sufficiently_accurate = False
            while not sufficiently_accurate:
                u_new, error = self.advance_intermediate(u, t, self.h)

                if self.verbose:
                    'u=%g at t=%g is' % (u_new, t + self.h),

                # Is u_new sufficiently accurate?
                tol = self.rtol*u_norm + self.atol
                accurate = error < tol
                if accurate or self.h <= min_step or self.h >= max_step:
                    sufficiently_accurate = True
                    u = u_new
                    t = t+self.h
                    self.t_all.append(t)
                    self.u_all.append(u)

                    if self.verbose >= 1:
                        print 'accepted, ',
                else:
                    if self.verbose >= 1:
                        print 'rejected, ',

                # Adjust time step
                u_new_norm = np.linalg.norm(u_new)
                if error > 1E-14:
                    h_new = dt*(tol/error)**(1./order)
                    if min_step <= h_new <= max_step:
                        self.h = h_new

                if t+dt > t_np1:  # fit last step so we hit t_np1 exactly
                    self.h = t_np1 - t

                if self.verbose >= 1:
                    print 'new h=%g' % self.g
        return u


class RKFehlberg(Adaptive):
    """The classical adaptive Runge-Kutta-Fehlberg method of order 4-5."""
    quick_description = "Adaptive Runge-Kutta-Fehlberg (4,5) method"

    _optional_parameters = Adaptive._optional_parameters


    def initialize_for_solve(self):
        Adaptive.initialize_for_solve(self)

    def advance(self):
        # auxilatory function to pick up the middle number from 3 floats
        def middle(x, y=.1, z=4.):
           return sorted([x, y, z])[1]

        f, n, rtol, atol = self.f, self.n, self.rtol, self.atol
        u_n, t_n, t_np1 = self.u[n], self.t[n], self.t[n+1]

        dt = t_np1 - t_n

        # default setting of step size
        min_step, max_step, h = \
            getattr(self, 'min_step', dt/1000.), \
            getattr(self, 'max_step', dt),\
            getattr(self, 'first_step', dt)

        # coefficients in Butcher tableau
        c = (1/4.,
             3/8.,
             3/32.,
             9/32.,
             12/13.,
             1932/2197.,
             -7200/2197.,
             7296/2197.,
             439/216.,
             -8.,
             3680/513.,
             -845/4104.,
             1/2.,
             -8/27.,
             2.,
             -3544/2565.,
             1859/4104.,
             -11/40.,
             1/360.,
             -128/4275.,
             -2197/75240.,
             1/50.,
             2/55.,
             25/216.,
             1408/2565.,
             2197/4104.,
             -1/5.)

        # u_i and t_i are intermediate steps between t_n and t_np1
        u_i = [u_n]; t_i = [t_n]

        while abs(t - t_n) < abs(t_np1 - t_n):
            u, t = u_i[-1], t_i[-1]

            # internal steps
            k1 = h*f(u, t)
            k2 = h*f(u+k1*c[0], t+h*c[0])
            k3 = h*f(u+k1*c[2]+k2*c[3], t+h*c[1])
            k4 = h*f(u+k1*c[5]+k2*c[6]+k3*c[7], t+h*c[4])
            k5 = h*f(u+k1*c[8]+k2*c[9]+k3*c[10]+k4*c[11], t+h)
            k6 = h*f(u+k1*c[13]+k2*c[14]+k3*c[15]+k4*c[16]+k5*c[17],
                     t+h*c[12])
            u_new = u + k1*c[23] + k3*c[24] + k4*c[25] + k5*c[26]

            # local error between 2 levels
            error = np.abs(k1*c[18] + k3*c[19] + k4*c[20] + \
                          k5*c[21] + k6*c[22])
            tol = rtol*np.abs(u_new) + atol
            # Error factor = local-error/error-tolerance
            rms = error/tol
            rms_norm = np.sqrt((np.sum(rms*rms))/self.neq)

            # Close enough or step size can not be reduced any more
            if rms_norm <= 1. or h <= min_step:
                u_i.append(u_new)
                t_i.append(t+h)

            # prevent the error of dividing absolute zero
            error = np.asarray([(1e-16 if x == 0. else x) for x in error]) \
                   if self.neq > 1 else (1e-16 if error == 0. else error)

            # Factor to adjust next step size
            s = (tol/(2*error))**0.25
            # Factor should be in a reasonable range[0.1,4.0]
            s = min(map(middle, s)) if self.neq > 1 else middle(s)

            # Step size should be in range [min_step, max_step]
            h = middle(h*s, y=min_step, z=max_step)

            # h should be set to 't_np1-t_i[-1]' at the last intern step.
            h = min(h, t_np1-t_i[-1])

        return u_new

class Ode_scipy(Adaptive):
    """
    Super class wrapper for scipy.integrate.ode classes.
    Existing solvers in subclasses are: Vode, Dopri5, Dop853.
    """

    _optional_parameters = Solver._optional_parameters + \
        ['jac', 'jac_kwargs', 'jac_args', 'atol', 'rtol',
         'first_step', 'max_step', 'nsteps']

    # Common scipy.integrate.ode arguments for subclass solvers
    _arglist_scipy = ['f', 'jac', 'atol', 'rtol', 'complex_valued',]
    # Name differences between this interface and scipy arguments
    _name_differences = {'adams_or_bdf': 'method',
                         'ml': 'lband', 'mu': 'rband'}

    def initialize(self):
        try:
            import scipy.integrate.ode as si_ode
        except ImportError:
            raise ImportError('The scipy package must be installed '\
                              'in order to use class %s' % \
                              self.__class__.__name__)

    def initialize_for_solve(self):
        # scipy specifies f and jac as f(t, y, *args) and jac(t, y)
        # while the present interface assumes
        # f(u, t, *f_args, **f_kwargs), jac(u, t, *jac_args, **jac_kwargs)
        self.f4scipy = lambda t, y: self.f(y, t)
        if self.jac is not None:
            # First wrap the user's Jacobian routine as we wrap f
            # (allow jac to return list of lists)
            jac = self.jac
            self.jac = lambda u, t: \
                np.asarray(jac(u, t, *self.jac_args, **self.jac_kwargs))
            # Then switch argument sequence
            self.jac4scipy = lambda t, y: self.jac(y, t)
        else:
            self.jac4scipy = None

        import scipy.integrate.ode as si_ode
        self.integrator = si_ode(self.f4scipy, jac=self.jac4scipy)
        self.scipy_arguments = {}

        # Extract common arguments, and prepare to be transferred to scipy
        for name in self._arglist_scipy:
            value = getattr(self, name, None)
            if value is not None:
                if name in self._parameters and \
                   name in Ode_scipy._name_differences:
                    # different names in scipy
                    name = Ode_scipy._name_differences[name]
                self.scipy_arguments[name] = value

        self.integrator = self.integrator.set_integrator(
            self.__class__.__name__.lower(), **self.scipy_arguments)
        self.integrator = self.integrator.set_initial_value(self.U0, self.t[0])
        Solver.initialize_for_solve(self)

    def advance(self):
        u, f, n, t = self.u, self.f, self.n, self.t
        u_new = self.integrator.integrate(t[n+1])
        if not self.integrator.successful():
            print 'Warning: %s call to scipy.integrate.ode.method.integrate was not successful' % self.__class__.__name__
        if len(u_new) == 1:
            return u_new[0]
        else:
            return u_new

class Vode(Ode_scipy):
    '''
    Wrapper for scipy.integrate.ode.vode, which is a wrapper for vode.f,
    which intends to solve initial value problems of stiff or nonstiff
    type. The well-known vode.f solver applies backward differential
    formulae for iteration.
    '''
    quick_description = "Adams/BDF Vode adaptive method (vode.f wrapper)"

    _optional_parameters = Ode_scipy._optional_parameters + \
                           ['adams_or_bdf', 'min_step', 'order']

    # argument list to be passed to scipy.ode for 'vode' method
    _arglist_scipy = ['atol', 'rtol', 'ml', 'mu', 'adams_or_bdf',
                      'with_jacobian', 'nsteps', 'first_step',
                      'min_step', 'max_step', 'order']

    def initialize_for_solve(self):
        # internal argument to be transferred to scipy
        self.with_jacobian = getattr(self, 'jac', None) is not None
        Ode_scipy.initialize_for_solve(self)

class Dopri5(Ode_scipy):
    """
    Wrapper for scipy.integrate.ode.dopri5, which applies the
    Dormand&Prince method of order 5(4), based on the Fortran
    implementation by Hairer and Wanner.
    See http://www.unige.ch/~hairer/software.html.
    """
    quick_description = "Dormand & Prince method of order 5(4) (scipy)"

    _optional_parameters = Ode_scipy._optional_parameters + \
        ['ifactor', 'dfactor', 'beta', 'safety']

    # argument list to be passed to scipy.ode for 'dopri5' method
    _arglist_scipy = ('atol','rtol','safety','ifactor','dfactor','beta',
                      'nsteps','first_step','max_step')

class Dop853(Ode_scipy):
    """
    Wrapper for scipy.integrate.ode.dop853, which applies the
    Dormand&Prince method of order 8(5,3), based on the Fortran
    implementation by Hairer and Wanner.
    See http://www.unige.ch/~hairer/software.html.
    """
    quick_description = "Adaptive Dormand & Prince method of order 8(5,3) (scipy)"

    _optional_parameters = Ode_scipy._optional_parameters + \
        ['ifactor', 'dfactor', 'beta', 'safety']

    # argument list to be passed to scipy.ode for 'dop853' method
    _arglist_scipy = ('atol','rtol','safety','ifactor','dfactor','beta',
                      'nsteps','first_step','max_step')

def list_all_solvers():
    """Return all solver classes in this package, excluding superclasses."""
    # Important: odespy.__init__.py must import all solver classes
    # into the namespace for this function to work properly.

    superclasses = ('Solver','Adaptive', 'PyDS', 'Ode_scipy', 'Odepack',
                    'RungeKutta1level', 'RungeKutta2level', 'SolverImplicit')
    import odespy
    class_members = inspect.getmembers(odespy, inspect.isclass)
    solvers = [solver[0] for solver in class_members \
                   if solver[0] not in superclasses]
    return solvers

def list_available_solvers():
    """Return all available solver classes in this package."""
    available_solvers = []
    import odespy
    all_solvers = list_all_solvers()
    for solvername in all_solvers:
        try:      # Try to initialize solvers with f is None
            method = eval('odespy.%s' % solvername)(None)
            available_solvers.append(solvername)
        except:
            try:
                # Try the exception of linearly solvers in ODEPACK
                # f is illegal for these solvers.
                method = eval('odespy.%s' % solvername)()
                available_solvers.append(solvername)
            except:
                # Failed to initialize this solver.
                # Perhaps the required dependency is not installed.
                pass
    return available_solvers



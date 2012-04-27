Odespy (ODE Software in Python) offers a unified interface to a
a large collection of software for solving systems of ordinary
differential equations (ODEs). There is also some support for
Differential Algebraic Equations (DAEs).

Odespy features the following collection of numerical methods and
implementations:

  * Pure Python implementations of classical explicit schemes such as
    the Forward Euler method (also called Euler);
    Runge-Kutta methods of 2nd, 3rd, and 4th order; Heun's method;
    Adams-Bashforth methods of 2nd, 3rd, and 4th order;
    Adams-Bashforth-Moulton methods of 2nd and 3rd order.

  * Pure Python implementations of classical implicit schemes such as
    Backward Euler; 2-step backward scheme; the theta rule;
    the Midpoint (or Trapezoidal) method.

  * Pure Python implementations of adaptive explicit Runge-Kutta
    methods of type Runge-Kutta-Fehlberg of order (4,5), Dormand-Prince
    of order (4,5), Cash-Karp of order (4,5), Bogacki-Shampine of order (2,3).

  * Wrappers for all FORTRAN solvers in `ODEPACK <http://www.netlib.org/odepack/>`_.

  * Wrappers for the wrappers of FORTRAN solvers in `scipy <http://www.scipy.org>`_:
    ``vode`` and ``zvode`` (adaptive Adams or BDF from `vode.f <http://www.netlib.org/ode/vode.f>`_);
    ``dopri5`` (adaptive Dormand-Prince method of order (4,5));
    ``dop853`` (adaptive Dormand-Prince method of order 8(5,3));
    ``odeint`` (adaptive Adams or BDF, basically the same as ``vode``, but in the implementation ``lsoda`` from `ODEPACK <http://www.netlib.org/odepack/>`_).

  * Wrapper for the Runge-Kutta-Chebyshev formulas of order 2 as
    offered by the well-known FORTRAN code `rkc.f <http://www.netlib.org/ode/rkc.f>`_.

  * Wrapper for the Runge-Kutta-Fehlberg method of
    order (4,5) as provided by the well-known FORTRAN code `rkf45.f <http://www.netlib.org/ode/rkf45.f>`_.

  * Wrapper for the Radau5 method as provided by the well-known FORTRAN code
    `radau5.f <http://www.unige.ch/~hairer/prog/stiff/radau5.f>`_.

  * Wrapper for some solvers in the `odelab <https://github.com/olivierverdier/odelab>`_ package.

The ODE problem can always be specified in Python, but for wrappers of
FORTRAN codes one can also implement the problem in FORTRAN and avoid
callback to Python.

Here is an example on the Odespy syntax::

        def f(u, t):
            """2x2 system for a van der Pool oscillator."""
            return [u[1], 3.*(1. - u[0]*u[0])*u[1] - u[0]]

        import odespy, numpy
        solver = odespy.Vode(f, rtol=0.0, atol=1e-6,
                             adams_or_bdf='adams', order=10)
        solver.set_initial_condition([2.0, 0.0])
        t_points = numpy.linspace(0, 30, 150)
        u, t = solver.solve(t_points)

        u0 = u[:,0]
        from matplotlib.pyplot import *
        plot(t, u0)
        show()

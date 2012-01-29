===========================================================
odesolvers
===========================================================

odesolvers is an unified Python interface that collects a 
variety of ODE software and methods to solve Initial Value 
Problems(IVP). 

User interface is designed in a minimal sense to simplify 
its usage, and enough to define all kinds of ODE systems. 
Currently no additional functionalities (like plotting, 
trajectory) are included in this package. The desired solution
 is returned as a standard Numpy array. Users can trivially 
apply their favored Python libraries for data analysis.

Up till now, there are 43 numerical methods and ODE software
integrated in odesolvers. Some of them are implementation for
simple algorithms directly in Python, while the others are 
wrappers for well-known ODE solvers. 

===========================================================
Structure
==========================================================
Each solver is implemented as a separate class, and all these 
classes are collected in a class hierarchy. The base class for
 this hierarchy is named as Solver. Generic features of ODE 
solvers are defined in this base class, which ensure the 
consistency and uniform of parameters, methods, and user 
interface in the whole hierarchy.

Currently odesolvers consists of the following files and directories:

README.txt
   Introduction, structure and installation.

setup.py
   Script for building and installing..

ODE.py
   Superclass Solver, which provides all the functionalities
 that are common for all solvers. Solvers included:
   Euler, Leapfrog, LeapfrogFiltered, Heun, RK2, RK3, RK4, 
   AdamsBashforth2, AdamsBashforth3, AdamsBashforth4,
   AdamsBashMoulton2, AdamsBashMoulton3, MidpointIter, 
   SymPy_odefun, SolverImplicit, BackwardEuler, Backward2Step, 
   ThetaRule, MidpointImplicit, Ode_scipy, Vode, Dopri5, Dop853.
   Adaptive, AdaptiveResidual, RKFehlberg.

rkc_rkf45.py
   Rkc, Rkf45.

OdePack.py
   Collection of ODEPACK solvers, including Lsode, Lsodes, 
Lsoda, Lsodar, Lsodi, Lsodis, Lsoibt.

RungeKutta.py
   Collection of Runge-Kutta methods, including RungeKutta2,
   RungeKutta3, RungeKutta4, ForwardEuler, DormandPrince, 
   Fehlberg, CashKarp, BogackiShampine, MyRungeKutta.

rkc/
   Signature file and Fortran code of rkc.f, which is available
 at http://www.netlib.org/ode/rkc.f.

rkf45/
   Signature file and Fortran code of rkf45.f, which is available
 at http://www.netlib.org/ode/rkf45.f.

odepack/
   Signature files and Fortran code of ODEPACK, which are available
 at http://www.netlib.org/odepack.

tests/
   Tests and demos. 

============================================================
Installation
============================================================

odesolvers is implemented as an open-source software and an
independent package in Python. The only mandatory dependencies
are Numpy and Python, with some other optional dependencies
for those wrappers. 

Mandatory dependency:
 
1. package Numpy, available at  numpy.scipy.org

2. Package nose is used for basic tests, which requires a version later
   than 0.10.0.

Optional dependencies:

3. package Sympy, available at sympy.org
   Related solver in odesolvers: SymPy_odefun

4. package Scipy, available at www.scipy.org
   Related solvers in odesolvers: Vode, Dopri5, Dop853

Even if these optional dependencies are not installed, users can
install this package smoothly, and make use of other basic
solvers without any extra installation.

Dependency to plot curves in demos:

5. In some demos, scitools.std.plot is used to plot curves. 
Package scitools is available at http://code.google.com/p/scitools/.


==========================================================
Information
==========================================================

1. odesolvers.list_solvers()
   A function to list all possible solvers in this module.

2. pydoc odesolvers.(desired solver_name)
   Example: pydoc odesolvers.RungeKutta4
   Print out useful information for desired solver, e.g. 
   parameters, methods, usage.

3. get_parameter_info(self, printInfo=True) 
   Print out information about legal parameters in current 
   solver, e.g. types, ranges, help_info, rtc..

==========================================================
To do
==========================================================

Currently, Lsodis and Lsoibt do not support user-supplied Fortran 
subroutines.





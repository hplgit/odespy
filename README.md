### What is Odespy?

Odespy (ODE Software in Python) offers a unified interface to a
large collection of software for solving systems of ordinary
differential equations (ODEs). There is also some support for
Differential Algebraic Equations (DAEs).

You can pronounce Odespy by first saying the acronym *ODE*, for
ordinary differential equations, in English (or in another language!)
and then adding "spy" as in the English word *spy*.

### How do I install Odespy?

#### Pip

Odespy requires Python version 2.7.  The simplest procedure for
installing Odespy is to use `pip`:


```
Terminal> sudo pip install -e git+https://github.com/hplgit/odespy.git#egg=odespy
```

#### Clone/fork repo

Alternatively, you can check out this repo and run `setup.py`:


```
Terminal> git clone git@github.com:hplgit/odespy.git
Terminal> cd odespy
Terminal> sudo python setup.py install
```

#### Turning off Fortran compilation

If you face problems with compiling the Fortran parts of Odespy,
or if you do not have a Fortran compiler, you can install without
any Fortran code:


```
Terminal> sudo python setup.py install --no-fortran
```

#### Problems with compiled Fortran modules and Anaconda

Problems with the compiled Fortran libraries appear on different systems.
If you use Anaconda Python, successfully build the Fortran modules, but
get error messages like `cannot find extension module _rkf`, the problem
may be related to GFortran 1.4. You may try


```
Terminal> cd /path/to/anaconda/lib/
Terminal> mv libgfortran.so.3.0.0 libgfortran.so.3.0.0.bak
Terminal> ln -sf /usr/lib/x86_64-linux-gnu/libgfortran.so.3.0.0 \
          libgfortran.so.3.0.0
```

See also [this link](http://stackoverflow.com/questions/9628273/libgfortran-version-gfortran-1-4-not-found).

#### Problems with compiled Fortran modules on Windows

There have been various problems with compiling Odespy on Windows,
usually related to the Fortran compiler.
One recommended technique is to rely on Anaconda on Windows,
install the `ming32` compiler, and then run


```
Terminal> python setup.py install build --compiler=ming32
```

This may give problems of the type


```
File "C:\Anaconda\lib\site-packages\numpy\distutils\fcompiler\gnu.py",
line 333, in get_libraries
raise NotImplementedError("Only MS compiler supported with gfortran on win64")
NotImplementedError: Only MS compiler supported with gfortran on win64
```

A remedy is to edit the `gnu.py` file and comment out the
`NotImplementedError`:


```python
else:
    #raise NotImplementedError("Only MS compiler supported with gfortran on win64")
    pass
```

### Install Pre-compiled Version with Conda

If you are using [Anaconda](https://store.continuum.io/cshop/anaconda/) or Miniconda and a build has been published for your system, then you can install a pre-compiled version of Odespy. Install using `conda` with the following command:


```
Terminal> conda install -c https://conda.binstar.org/rothnic odespy
```

If a pre-compiled package is not available for your system, and you go
through the effort of installing the compilation tools, you can build
and upload a conda package as follows from the command line:


```shell
# Set environment variable in windows to build for python 2
# Use 34 instead of 27, to compile for python 3
set CONDA_PY=27

# Build the package
conda build odespy

# Upload the package
binstar upload <PATH to the build file>\odespy-<version>-<dependency versions>.tar.bz2
```

*Notice.* 
You may have to add/modify the conda build scripts to support your platform type. This has only been tested for Windows so far.




### Contents of Odespy

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
  * Wrappers for all FORTRAN solvers in [`ODEPACK`](http://www.netlib.org/odepack).
  * Wrappers for the wrappers of FORTRAN solvers in [`scipy`](http://www.scipy.org):
    `vode` and `zvode` (adaptive Adams or BDF from [`vode.f`](http://www.netlib.org/ode/vode.f));
    `dopri5` (adaptive Dormand-Prince method of order (4,5));
    `dop853` (adaptive Dormand-Prince method of order 8(5,3));
    `odeint` (adaptive Adams or BDF, basically the same as `vode`, but in the implementation `lsoda` from [`ODEPACK`](http://www.netlib.org/odepack/)).
  * Wrapper for the Runge-Kutta-Chebyshev formulas of order 2 as
    offered by the well-known FORTRAN code [`rkc.f`](http://www.netlib.org/ode/rkc.f).
  * Wrapper for the Runge-Kutta-Fehlberg method of
    order (4,5) as provided by the well-known FORTRAN code [`rkf45.f`](http://www.netlib.org/ode/rkf45.f).
  * Wrapper for the Radau5 method as provided by the well-known FORTRAN code
    [`radau5.f`](http://www.unige.ch/~hairer/prog/stiff/radau5.f).
    There have been some unidentified problems with running this solver (segmentation fault).
  * Wrapper for some solvers in the [`odelab`](https://github.com/olivierverdier/odelab).

The ODE problem can always be specified in Python, but for wrappers of
FORTRAN codes one can also implement the problem in FORTRAN and avoid
callback to Python.

*Warning: Potential problems with FORTRAN codes.* 
Some users have faced problems with some of the FORTRAN-based solvers
(typically segmentation fault errors), mostly `radau5`, but also `rkf45`.
It seems that Odespy's Python interface to `radau5` is broken.



### How do I use Odespy?

Here is an example on the Odespy syntax


```python
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
```

An incomplete [tutorial](http://hplgit.github.io/odespy/doc/tutorial/html/index.html) is under
development and explains much more of the syntax and provides many
examples.

### I have used Odespy in a publication - how do I cite it?

Please cite this GitHub repository:

> H. P. Langtangen and L. Wang. Odespy software package.
> URL: <https://github.com/hplgit/odespy.> 2015



BibTeX entry:


```
@misc{odespy,
  title = {{Odespy} Software Package},
  author = {H. P. Langtangen and L. Wang},
  url = {https://github.com/hplgit/odespy},
  key = {odespy},
  note = {\url{https://github.com/hplgit/odespy}},
  year = {2015},
}
```

[Publish](https://bitbucket.org/logg/publish) entry:


```
* misc
** {Odespy} Software Package
   key:       odespy
   author:    H. P. Langtangen, L. Wang
   url:       https://github.com/hplgit/odespy
   status:    published
   sortkey:   Odespy
   year:      2015
   note:      \url{https://github.com/hplgit/odespy}
   entrytype: misc
```


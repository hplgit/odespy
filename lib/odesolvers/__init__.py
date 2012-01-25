from os.path import join
from numpy.testing import rundocs, run_module_suite


# Files in this module
from ODE import *
from RungeKutta import *
from odepack import *
from rkc_rkf45 import *

# Doctests
path = __path__[0]
rundocs(join(path, 'ODE.py'))
rundocs(join(path,'RungeKutta.py'))

# Basic tests
path = join(path, 'tests')
run_module_suite(join(path, 'test_basics.py'))

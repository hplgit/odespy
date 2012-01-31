from os.path import join
from numpy.testing import rundocs, run_module_suite



from ODE import *
from RungeKutta import *
try:
    from odepack import *
except ImportError:
    print 'No odepack module'
try:
    from rkc_rkf45 import *
except ImportError:
    print 'No support for RKF45 and RKC methods in Fortran (rkc_rkf45 module)'


# Update doc strings with common info
class_, doc_str, classname = None, None, None
classes = [item[0] for item in locals().items() \
               if inspect.isclass(item[1])]
for classname in classes:
    class_ = eval(classname)
    doc_str = getattr(class_, '__doc__')
    setattr(class_, '__doc__', doc_str + doc_string_table_of_parameters(class_))
del class_, doc_str, classname  # do not pollute namespace

if __name__ == '__main__':
    # Doctests
    path = __path__[0]
    rundocs(join(path, 'ODE.py'))
    rundocs(join(path,'RungeKutta.py'))

    # Basic tests
    path = join(path, 'tests')
    run_module_suite(join(path, 'test_basics.py'))

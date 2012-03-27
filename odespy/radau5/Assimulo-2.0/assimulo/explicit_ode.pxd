
from ode cimport ODE
import numpy as N
cimport numpy as N


cdef class Explicit_ODE(ODE):

    cpdef _simulate(self, double t0, double tfinal,N.ndarray output_list,int ONE_STEP, int INTERPOLATE_OUTPUT,int TIME_EVENT, int STEP_EVENT)
    

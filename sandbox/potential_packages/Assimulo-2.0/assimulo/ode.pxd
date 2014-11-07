import numpy as N
cimport numpy as N

cdef class ODE:
    cdef public dict options, solver_options, problem_info
    cdef public dict supports, statistics
    
    cdef public list event_data
    
    cdef public object problem
    
    cdef public double t, t0
    cdef public N.ndarray y,yd, p
    cdef public N.ndarray y0, yd0, p0, sw0
    
    #cdef public list t,y,yd,p,sw_cur
    cdef public list t_sol, y_sol, yd_sol, p_sol, sw
        
    cpdef log_message(self, message, int level)
    cpdef log_event(self, double time, object event_info, int level)
    cpdef simulate(self, double tfinal, int ncp=*, object ncp_list=*)
    cpdef get_options(self)
    cpdef get_supports(self)
    cpdef get_statistics(self)
    cpdef get_event_data(self)
    cpdef print_event_data(self)
    cpdef finalize(self)
    cpdef initialize(self)
    cdef _reset_solution_variables(self)

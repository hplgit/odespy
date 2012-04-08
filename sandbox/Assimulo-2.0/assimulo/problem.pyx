#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2011 Modelon AB
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import numpy as N
cimport numpy as N

include "constants.pxi" #Includes the constants (textual include)

realtype = N.float

cdef class cProblem:
    cdef public int _sensitivity_result
    
    name = '---'
    
    cpdef initialize(self, solver):
        """
        Method for specializing initiation.
        """
        solver.log_message("No initialization defined for the problem.", LOUD)
    
    cpdef reset(self):
        """
        Resets a problem to its default values.
        """
        pass
        
    cpdef handle_event(self, object solver, event_info):
        """
        Defines how to handle a discontinuity. This functions gets called when
        a discontinuity has been found in the supplied event functions. The solver
        is the solver attribute while the event_info is a list of length 2 where
        the first element is a list containing information about state events and
        the second element is a boolean for indicating if there have been an time
        event. If there have not been a state event the first element is an empty
        list. The state event list contains a set of integers of values (-1,0,1),
        the values indicates which state event have triggered (determined from 
        state_event(...) ) and the value indicates to where the state event is 'headed'.
        """
        solver.log_message("No event handling defined.", NORMAL)
    
    cpdef finalize(self,object solver):
        """
        Method for specifying the finalization options when the simulation have
        finished.
        """
        solver.log_message("No finalization defined for the problem.", LOUD)

cdef class cImplicit_Problem(cProblem):
    
    def __init__(self, object res=None, y0=None, yd0=None,double t0=0.0, p0=None, sw0=None):
        
        if res != None:
            self.res = res
        if y0!=None:
            self.y0  = None if y0 is None else (N.array(y0,dtype=realtype) if len(N.array(y0,dtype=realtype).shape)>0 else N.array([y0],dtype=realtype))
        if yd0!=None:
            self.yd0 = None if yd0 is None else (N.array(yd0,dtype=realtype) if len(N.array(yd0,dtype=realtype).shape)>0 else N.array([yd0],dtype=realtype))
        if p0!=None:
            self.p0 = None if p0 is None else (N.array(p0,dtype=realtype) if len(N.array(p0,dtype=realtype).shape)>0 else N.array([p0],dtype=realtype))
        if sw0!=None:
            self.sw0 = None if sw0 is None else (N.array(sw0,dtype=bool) if len(N.array(sw0,dtype=bool).shape)>0 else N.array([sw0],dtype=bool))
        self.t0  = t0
    
    cpdef handle_result(self, solver, double t, N.ndarray[double, ndim=1] y, N.ndarray[double, ndim=1] yd):
        """
        Method for specifying how the result is to be handled. As default the
        data is stored in three vectors, solver.(t/y/yd).
        """
        cdef int i = 0
        
        solver.t_sol.extend([t])
        solver.y_sol.extend([y])
        solver.yd_sol.extend([yd])
        
        #Store sensitivity result (variable _sensitivity_result are set from the solver by the solver)
        if self._sensitivity_result == 1:
            for i in range(solver.problem_info["dimSens"]):
                solver.p_sol[i] += [solver.interpolate_sensitivity(t, i=i)]
        
    cpdef res_internal(self, N.ndarray[double, ndim=1] res, double t, N.ndarray[double, ndim=1] y, N.ndarray[double, ndim=1] yd):
        try:
            res[:] = self.res(t,y,yd)
        except:
            return ID_FAIL
        return ID_OK
    
cdef class cExplicit_Problem(cProblem):
    
    def __init__(self, object rhs=None, y0=None,double t0=0.0, p0=None, sw0=None):
        
        if rhs != None:
            self.rhs = rhs
        if y0!=None:
            self.y0  = None if y0 is None else (N.array(y0,dtype=realtype) if len(N.array(y0,dtype=realtype).shape)>0 else N.array([y0],dtype=realtype))
        if p0!=None:
            self.p0 = None if p0 is None else (N.array(p0,dtype=realtype) if len(N.array(p0,dtype=realtype).shape)>0 else N.array([p0],dtype=realtype))
        if sw0!=None:
            self.sw0 = None if sw0 is None else (N.array(sw0,dtype=bool) if len(N.array(sw0,dtype=bool).shape)>0 else N.array([sw0],dtype=bool))
        self.t0  = t0
    
    cpdef handle_result(self, solver, double t, N.ndarray[double, ndim=1] y):
        """
        Method for specifying how the result is to be handled. As default the
        data is stored in three vectors, solver.(t/y).
        """
        cdef int i = 0
        
        solver.t_sol.extend([t])
        solver.y_sol.extend([y])
        
        #Store sensitivity result (variable _sensitivity_result are set from the solver by the solver)
        if self._sensitivity_result == 1:
            for i in range(solver.problem_info["dimSens"]):
                solver.p_sol[i] += [solver.interpolate_sensitivity(t, i=i)]
        
    cpdef int rhs_internal(self, N.ndarray[double, ndim=1] yd, double t, N.ndarray[double, ndim=1] y):
        try:
            yd[:] = self.rhs(t,y)
        except:
            return ID_FAIL
        return ID_OK

class Implicit_Problem(cImplicit_Problem):
    """
        Problem for our implicit integrators (DAEs). A problem
        consists of the residual function and some initial conditions.
        
        Parameters ::
          
            res   
                Function that calculates the residual. Depending on
                the problem and the support of the solver, this function can
                have the following input parameters.
                
                    res(t,y,yd)      - Normal DAE
                    res(t,y,yd,sw)   - An DAE with different modes, sw is a list of
                                       switches (boolean list) which should be held
                                       constant during the integration and only be
                                       changed when an event have occured. Used together
                                       with event functions.
                    res(t,y,yd,p)    - An DAE with parameters for which sensitivities
                                       should be calculated.
                    res(t,y,yd,sw,p) - An DAE with both parameters and switches.
                    
                    Returns:
                        A numpy array of size len(y).
            y0
                Defines the starting values of y0.
            yd0
                Defines the starting values of yd0.
            t0
                Defines the starting time.
            sw0 (Depending on if the solver supports state events)
                Defines the starting values of the switches. 
                Should be a list of booleans.
            p0 (Depending on if the solver supports sensitivity calculations)
                Parameters for which sensitivites are to be calculated
                
        Parameters (optionally contained in class) ::
        
            algvar
                Defines the differential and algebraic components of the problem.
                Should be a list of integers. For more information, see the
                property algvar in IDA.
        
        Available (optional) options (depending on the solver support)::
        
            def state_events(self ,t ,y ,yd, sw)
                Defines the event (root) functions.
                
                Returns:
                    A numpy array.
                
            def time_events(self, t, y, yd, sw)
                Defines the time events. This function should return
                the next time-point for a time event. At a time-event
                the usual method handle_event is called for the specific
                handling. If there are no more time events. This function
                should return None.
                
                Returns:
                    Float
                        The time-point for the next time-event.
                    None
                        No time-event.
                
            def jac(self, c, t, y, yd, sw)
                Defines the Jacobian, which should be of the form
                J = dF/dx + c*dF/dx'.
                
                Returns:
                    A numpy array of size len(y)*len(y).
                    
            def handle_result(self, solver, t, y, yd)
                Method for specifying how the result is to be handled. 
                As default the data is stored in three vectors, solver.(t_sol/y_sol/yd_sol). 
                If the problem to be solved also involve sensitivities these results are
                stored in p_sol
                
            def handle_event(self, object solver, event_info):
                Defines how to handle a discontinuity. This functions gets called when
                a discontinuity has been found in the supplied event functions. The solver
                is the solver attribute while the event_info is a list of length 2 where
                the first element is a list containing information about state events and
                the second element is a boolean for indicating if there have been an time
                event. If there have not been a state event the first element is an empty
                list. The state event list contains a set of integers of values (-1,0,1),
                the values indicates which state event have triggered (determined from 
                state_event(...) ) and the value indicates to where the state event is 'headed'.
    """
    pass
    
class Explicit_Problem(cExplicit_Problem):
    """
        Problem for our explicit integrators (ODEs). A problem
        consists of the right-hand-side and some initial conditions.
 
        Parameters::
            
            rhs 
                Function that calculates the right-hand-side. Depending on
                the problem and the support of the solver, this function can
                have the following input parameters.
                
                    rhs(t,y)      - Normal ODE
                    rhs(t,y,sw)   - An ODE with different modes, sw is a list of
                                    switches (boolean list) which should be held
                                    constant during the integration and only be
                                    changed when an event have occured. Used together
                                    with event functions.
                    rhs(t,y,p)    - An ODE with parameters for which sensitivities
                                    should be calculated.
                    rhs(t,y,sw,p) - An ODE with both parameters and switches.
                    
                    Returns:
                        A numpy array of size len(y).
            
            y0
                Defines the starting values of y0
            t0
                Defines the starting time
            sw0 (Depending on if the solver supports state events)
                Defines the starting values of the switches. 
                Should be a list of booleans.
            p0 (Depending on if the solver supports sensitivity calculations)
                Parameters for which sensitivites are to be calculated
        
        Available (optional) options (depending on the solver support)::
        
            def state_events(self ,t ,y, sw)
                Defines the event (root) functions.
                
                Returns:
                    A numpy array.
                    
            def time_events(self, t, y, sw)
                Defines the time events. This function should return
                the next time-point for a time event. At a time-event
                the usual method handle_event is called for the specific
                handling. If there are no more time events. This function
                should return None.
                
                Returns:
                    Float
                        The time-point for the next time-event.
                    None
                        No time-event.
                
            def jac(self, t, y, sw=None)
                Defines the jacobian. J=df/dx.
                
                Returns:
                    A numpy matrix of size len(y)*len(y).
                    
            def jacv(self, t, y, fy, v)
                Defines a Jacobian Vector product. df/dx*v.
                
                Returns:
                    A numpy vector of size len(y).
            
            def handle_result(self, solver, t, y)
                Method for specifying how the result is to be handled. 
                As default the data is stored in two vectors, solver.(t_sol/y_sol). If
                the problem to be solved also involve sensitivities these results are
                stored in p_sol
                
            def handle_event(self, object solver, event_info):
                Defines how to handle a discontinuity. This functions gets called when
                a discontinuity has been found in the supplied event functions. The solver
                is the solver attribute while the event_info is a list of length 2 where
                the first element is a list containing information about state events and
                the second element is a boolean for indicating if there have been an time
                event. If there have not been a state event the first element is an empty
                list. The state event list contains a set of integers of values (-1,0,1),
                the values indicates which state event have triggered (determined from 
                state_event(...) ) and the value indicates to where the state event is 'headed'.
    """
    pass

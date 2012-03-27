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

from ode cimport ODE     
from problem import Explicit_Problem

import pylab as P
import itertools
import numpy as N
cimport numpy as N

from exception import *

include "constants.pxi" #Includes the constants (textual include)

realtype = N.float

cdef class Explicit_ODE(ODE):
    """
    Baseclass for our explicit ODE integrators.
    """
    def __init__(self, problem):
        """
        Initiates the solver.
        
            Parameters::
            
                problem     
                            - The problem to be solved. Should be an instance
                              of the 'Explicit_Problem' or 'cExplicit_Problem'
                              class.
        """
        ODE.__init__(self, problem) #Sets general attributes
        
        if isinstance(problem, Explicit_Problem):
            self.problem = problem
        else:
            raise Explicit_ODE_Exception('The problem needs to be a subclass of a Explicit_Problem.')
        
        #Check the dimension of the state event function
        if self.problem_info["state_events"]:
            self.problem_info["dimRoot"] = len(problem.state_events(self.t0,self.y0, self.sw0))
        
        self.t = self.t0
        self.y = self.y0.copy()
            
    def reset(self):
        """
        Resets the problem. If the problem is defined with a reset method, its called
        and then the method re_init. The re_init method is called with the initial
        values set in the problem, problem.t0 and problem.y0.
        
        """
        self.problem.reset()

        self.re_init(self.t0, self.y0, self.sw0 if self.problem_info["switches"] else None)
        
    def re_init(self,t0, y0, sw0=None):
        """
        Reinitiates the solver.
        
            Parameters::
                
                t0  
                    - The initial time.
                y0  
                    - The initial values for the states
                
        See information in the __init__ method.
        """
        if len(self.y) != len(y0):
            raise Explicit_ODE_Exception('y0 must be of the same length as the original problem.')
        
        #Set the new values as the current values
        self.t = float(t0)
        self.y = N.array(y0) if len(N.array(y0).shape)>0 else N.array([y0])
        
        if sw0 != None:
            self.sw = (N.array(sw0,dtype=N.bool) if len(N.array(sw0,dtype=N.bool).shape)>0 else N.array([sw0],dtype=N.bool)).tolist()

    cpdef _simulate(self, double t0, double tfinal,N.ndarray output_list,int ONE_STEP, int INTERPOLATE_OUTPUT,
                 int TIME_EVENT, int STEP_EVENT):
        """
        INTERNAL FUNCTION, FOR SIMULATION USE METHOD SIMULATE.
        
        Calls the integrator to perform the simulation over the given time-interval.
        If a second call to simulate is performed, the simulation starts from the last
        given final time.
        
            Parameters::
            
                tfinal  
                        - Final time for the simulation
                
                        - Should be a float or integer greater than the initial time.
                        
                ncp     
                        - Default '0'. Number of communication points where the 
                          solution is returned. If '0', the integrator will return 
                          at its internal steps.
                          
                        - Should be an integer.
                          
                    Example:
                    
                        __call__(10.0, 100), 10.0 is the final time and 100 is the number
                                             communication points.
        """
        cdef double t_log, tevent
        cdef int flag, output_index
        cdef dict opts
        
        y0 = self.y
        t_logg = t0

        #Log the first point
        self.problem.handle_result(self,t0,y0)

        #Reinitiate the solver
        flag_initialize = True

        #Start flag
        flag = ID_OK
        tevent = tfinal
        
        #Internal solver options
        opts = {}
        opts["initialize"] = flag_initialize
        opts["output_list"] = output_list
        opts["output_index"] = 0
        output_index = 0
        
        while (flag == ID_COMPLETE and tevent == tfinal) is False:
            
            #Time event function is specified.
            if TIME_EVENT == 1:
                tret = self.problem.time_events(self.t, self.y, self.sw)
                tevent = tfinal if tret is None else (tret if tret < tfinal else tfinal)
            else:
                tevent = tfinal
            
            if ONE_STEP == 1:
                #Run in one step mode
                [flag, t, y]         = self.step(self.t, self.y, tevent, opts)
                self.t, self.y = t, y.copy()
                
                #Store data depending on situation
                if INTERPOLATE_OUTPUT == 1:
                    try:
                        while output_list[output_index] <= t:
                            self.problem.handle_result(self, output_list[output_index], self.interpolate(output_list[output_index]))
                        
                            #Last logging point
                            t_logg = output_list[output_index]
                            
                            output_index = output_index+1
                    except IndexError:
                        pass
                else:
                    self.problem.handle_result(self,t,y)
                    
                    #Last logging point
                    t_logg = self.t
                    
                if STEP_EVENT == 1: #If the option completed step is set.
                    flag_initialize = self.problem.step_events(self)#completed_step(self)
                else:
                    flag_initialize = False
            else:
                #Run in Normal mode
                flag, tlist, ylist = self.integrate(self.t, self.y, tevent, opts)
                self.t, self.y = tlist[-1], ylist[-1].copy()
                
                #Store data
                map(self.problem.handle_result,itertools.repeat(self,len(tlist)), tlist, ylist)
                
                #Last logging point
                t_logg = self.t
                
                #Initialize flag to false
                flag_initialize = False
            
            #Event handling
            if flag == ID_EVENT or (flag == ID_COMPLETE and tevent != tfinal): #Event have been detected
                
                #Get and store event information
                event_info = [[],flag == ID_COMPLETE]
                if flag == ID_EVENT:
                    event_info[0] = self.state_event_info()
                
                #Log the information
                self.log_event(self.t, event_info, NORMAL)
                self.log_message("A discontinuity occured at t = %e."%self.t,NORMAL)
                self.log_message("Current Switches: " + str(self.sw), LOUD)
                self.log_message('Event info: ' + str(event_info), LOUD) 
                
                #Print statistics
                self.print_statistics(LOUD)
                
                try:
                    self.problem.handle_event(self, event_info) #self corresponds to the solver
                except TerminateSimulation: #Terminating the simulation after indication from handle event
                    self.log_message("Terminating simulation at t = %f after signal from handle_event."%self.t, NORMAL)
                    break
                
                flag_initialize = True

            #Update options
            opts["initialize"] = flag_initialize
            
            #Logg after the event handling if there was a communication point there.
            if flag_initialize and t_logg == self.t: 
                self.problem.handle_result(self, self.t, self.y)
            
    
    def plot(self, mask=None, **kwargs):
        """
        Plot the computed solution.
        
            Parameters::
            
                mask    
                        - Default 'None'. Used to determine which variables that is to be plotted.
                          Used as a list of integers, ones represents the variable that is to be
                          plotted and zeros that is not. 
                        
                        - Should be a list of integers.
                        
                            Example:
                                mask = [1,0] , plots the first variable.
                
                **kwargs
                        - See http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.plot
                          for information about the available options for **kwargs.
        """
        if len(self.t_sol) > 0:
            P.xlabel('time')
            P.ylabel('state')
            P.title(self.problem.name)
            
            if not mask:
                P.plot(self.t_sol, self.y_sol, **kwargs)
            else:
                if not isinstance(mask, list):
                    raise Explicit_ODE_Exception('Mask must be a list of integers')
                if not len(mask)==len(self.y_sol[-1]):
                    raise Explicit_ODE_Exception('Mask must be a list of integers of equal length as '\
                                                 'the number of variables.')
                for i in range(len(mask)):
                    if mask[i]:
                        P.plot(self.t_sol, N.array(self.y_sol)[:,i],**kwargs)
            
            
            P.show()
        else:
            self.log_message("No result for plotting found.",NORMAL)

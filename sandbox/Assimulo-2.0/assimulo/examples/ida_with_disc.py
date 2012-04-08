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
import pylab as P
import nose
from assimulo.solvers.sundials import IDA
from assimulo.problem import Implicit_Problem

"""
An example with event iteration and with three switches.

t=0     , [False, True, True]   (Start of simulation)
t=1 (1) , [False, True, False]  (Found a root at t=1)
t=1 (2) , [False, False, False] (Second iteration at t=1)
t=1 (3) , [True, False, False]  (Third iteration at t=1)
t=10    , [True, False, False]  (End of simulation)

"""

#Extend Assimulos problem definition
class Extended_Problem(Implicit_Problem):
    
    #Sets the initial conditons directly into the problem
    y0 = [0.0, -1.0, 0.0]
    yd0 = [-1.0, 0.0, 0.0]
    sw0 = [False,True,True]
    algvar = [1.0, 0.0, 0.0] #Determine which variables are differential and algebraic
    
    
    #The residual
    def res(self,t,y,yd,sw):
        """
        This is our function we are trying to simulate. During simulation
        the parameter sw should be fixed so that our function is continuous
        over the interval. The parameters sw should only be changed when the
        integrator has stopped.
        """
        res_0 = -yd[0] + (1.0 if sw[0] else -1.0)
        res_1 = -y[1] + (-1.0 if sw[1] else 3.0)
        res_2 = -y[2] + (0.0 if sw[2] else 2.0)
        
        return N.array([res_0,res_1,res_2])

    #Sets a name to our function
    name = 'Function with consistency problem'
    
    #The event function
    def state_events(self,t,y,yd,sw):
        """
        This is our function that keeps track of our events. When the sign
        of any of the events has changed, we have an event.
        """
        event_0 = y[1] - 1.0 
        event_1 = -y[2] + 1.0
        event_2 = -t + 1.0
        
        return N.array([event_0,event_1,event_2])
    
    
    #Responsible for handling the events.
    def handle_event(self, solver, event_info):
        """
        Event handling. This functions is called when Assimulo finds an event as
        specified by the event functions.
        """
        event_info = event_info[0] #We only look at the state events information.
        while True: #Event Iteration
            self.event_switch(solver, event_info) #Turns the switches
            
            b_mode = self.state_events(solver.t, solver.y, solver.yd, solver.sw)
            self.init_mode(solver) #Pass in the solver to the problem specified init_mode
            a_mode = self.state_events(solver.t, solver.y, solver.yd, solver.sw)
            
            event_info = self.check_eIter(b_mode, a_mode)
                
            if not True in event_info: #Breaks the iteration loop
                break
    
    #Helper function for handle_event
    def event_switch(self, solver, event_info):
        """
        Turns the switches.
        """
        for i in range(len(event_info)): #Loop across all event functions
            if event_info[i] != 0:
                solver.sw[i] = not solver.sw[i] #Turn the switch
        
    #Helper function for handle_event
    def check_eIter(self, before, after):
        """
        Helper function for handle_event to determine if we have event
        iteration.
        
            Input: Values of the event indicator functions (state_events)
            before and after we have changed mode of operations.
        """
        
        eIter = [False]*len(before)
        
        for i in range(len(before)):
            if (before[i] < 0.0 and after[i] > 0.0) or (before[i] > 0.0 and after[i] < 0.0):
                eIter[i] = True
                
        return eIter
    
    def init_mode(self, solver):
        """
        Initialize the DAE with the new conditions.
        """
        solver.make_consistent('IDA_YA_YDP_INIT') #Calculate new initial conditions.
                                                   #see SUNDIALS IDA documentation
                                                   #on the option 'IDA_YA_YDP_INIT'



def run_example(with_plots=True):
    
    #Create an instance of the problem
    iter_mod = Extended_Problem() #Create the problem

    iter_sim = IDA(iter_mod) #Create the solver
    
    iter_sim.verbosity = 0
    iter_sim.continuous_output = True

    #Simulate
    t, y, yd = iter_sim.simulate(10.0,1000) #Simulate 10 seconds with 1000 communications points
    
    #Basic test
    nose.tools.assert_almost_equal(y[-1][0],8.0)
    nose.tools.assert_almost_equal(y[-1][1],3.0)
    nose.tools.assert_almost_equal(y[-1][2],2.0)
    
    #Plot
    if with_plots:
        P.plot(t,y)
        P.show()
    
if __name__=="__main__":
    run_example()
    

    

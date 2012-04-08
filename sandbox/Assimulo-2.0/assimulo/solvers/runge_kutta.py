#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2010 Modelon AB
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

from assimulo.ode import *
from assimulo.explicit_ode import Explicit_ODE

from assimulo.exception import *

class RungeKutta34(Explicit_ODE):
    """
    Adaptive Runge-Kutta of order four.
    
    Obs. Step rejection not implemented.
    """
    def __init__(self, problem):
        """
        Initiates the solver.
        
            Parameters::
            
                problem     
                            - The problem to be solved. Should be an instance
                              of the 'Explicit_Problem' class.                       
        """
        Explicit_ODE.__init__(self, problem) #Calls the base class
        
        #Solver options
        self.options["atol"] = 1.0e-6
        self.options["rtol"] = 1.0e-6
        self.options["inith"] = 0.01
        self.options["maxsteps"] = 10000
        
        #Internal temporary result vector
        self.Y1 = N.array([0.0]*len(self.y0))
        self.Y2 = N.array([0.0]*len(self.y0))
        self.Y3 = N.array([0.0]*len(self.y0))
        self.Y4 = N.array([0.0]*len(self.y0))
        self.Z3 = N.array([0.0]*len(self.y0))
        
        #RHS-Function
        self.f = problem.rhs_internal
        
        #Solver support
        self.supports["one_step_mode"] = True
        
        #Internal values
        # - Statistic values
        self.statistics["nsteps"] = 0 #Number of steps
        self.statistics["nfcn"] = 0 #Number of function evaluations
    
    def initialize(self):
        #Reset statistics
        for k in self.statistics.keys():
            self.statistics[k] = 0
    
    def _set_initial_step(self, initstep):
        try:
            initstep = float(initstep)
        except (ValueError, TypeError):
            raise Explicit_ODE_Exception('The initial step must be an integer or float.')
        
        self.options["inith"] = initstep
        
    def _get_initial_step(self):
        """
        This determines the initial step-size to be used in the integration.
        
            Parameters::
            
                inith    
                            - Default '0.01'.
                            
                            - Should be float.
                            
                                Example:
                                    inith = 0.01
        """
        return self.options["inith"]
        
    inith = property(_get_initial_step,_set_initial_step)
    
    def _set_atol(self,atol):

        try:
            atol_arr = N.array(atol, dtype=float)
            if (atol_arr <= 0.0).any():
                raise Explicit_ODE_Exception('Absolute tolerance must be a positive float or a float vector.')
        except (ValueError,TypeError):
            raise Explicit_ODE_Exception('Absolute tolerance must be a positive float or a float vector.')
        if atol_arr.size == 1:
            self.options["atol"] = float(atol)
        elif atol_arr.size == len(self.y):
            self.options["atol"] = [float(x) for x in atol]
        else:
            raise Explicit_ODE_Exception('Absolute tolerance must be a float vector of same dimension as the problem or a scalar.')

    def _get_atol(self):
        """
        Sets the absolute tolerance to be used in the integration.
        
            Parameters::
            
                atol    
                            - Default 1.0e-6.
                            
                            - Should be float or an array/list of len(y)
                            
                                Example:
                                    atol=1.e5
                                    atol=[1.e-5,1.e-4]
        """
        return self.options["atol"]
    
    atol = property(_get_atol,_set_atol)
    
    def _set_rtol(self, rtol):
        try:
            rtol = float(rtol)
        except (TypeError,ValueError):
            raise Explicit_ODE_Exception('Relative tolerance must be a float.')
        if rtol <= 0.0:
            raise Explicit_ODE_Exception('Relative tolerance must be a positive (scalar) float.')
        self.options["rtol"] = rtol
            
    def _get_rtol(self):
        """
        The relative tolerance to be used in the integration.
        
            Parameters::
            
                rtol    
                            - Default 1.0e-6
                            
                            - Should be a float.
        """
        return self.options["rtol"]
    
    rtol = property(_get_rtol, _set_rtol)
    
    def _get_maxsteps(self):
        """
        The maximum number of steps allowed to be taken to reach the
        final time.
        
            Parameters::
            
                maxsteps
                            - Default 10000
                            
                            - Should be a positive integer
        """
        return self.options["maxsteps"]
    
    def _set_maxsteps(self, max_steps):
        try:
            max_steps = int(max_steps)
        except (TypeError, ValueError):
            raise Explicit_ODE_Exception("Maximum number of steps must be a positive integer.")
        self.options["maxsteps"] = max_steps
    
    maxsteps = property(_get_maxsteps, _set_maxsteps)
    
    def step(self, t, y, tf, opts):
        initialize = opts["initialize"]
        
        if initialize:
            self.solver_iterator = self._iter(t,y,tf)

        return self.solver_iterator.next()
    
    def integrate(self, t, y, tf, opts):
        """
        Integrates (t,y) values until t > tf
        """
        [flags, tlist, ylist] = zip(*list(self._iter(t, y, tf)))
        
        return flags[-1], tlist, ylist
    
    def _iter(self,t,y,tf):
        maxsteps = self.options["maxsteps"]
        h = self.options["inith"]
        h = min(h, N.abs(tf-t))
        
        for i in range(maxsteps):
            if t+h < tf:
                t, y, error = self._step(t, y, h)
                self.statistics["nsteps"] += 1
                yield ID_PY_OK, t,y
                h=self.adjust_stepsize(h,error)
                h=min(h, N.abs(tf-t))
            else:
                break
        else:
            raise Explicit_ODE_Exception('Final time not reached within maximum number of steps')
            
        t, y, error = self._step(t, y, h)
        self.statistics["nsteps"] += 1
        yield ID_PY_COMPLETE, t, y

    def adjust_stepsize(self, h, error):
        """
        Adjusts the stepsize.
        """
        fac=min((1.0/error)**(1.0/4.0),2.)
        h *= fac
        
        return h
        
    def _step(self, t, y, h):
        """
        This calculates the next step in the integration.
        """
        self.statistics["nfcn"] += 5
        
        scaling = N.array(abs(y)*self.rtol + self.atol) # to normalize the error 
        f = self.f
        
        f(self.Y1, t, y)
        f(self.Y2, t + h/2., y + h*self.Y1/2.)
        f(self.Y3, t + h/2., y + h*self.Y2/2.)
        f(self.Z3, t + h, y - h*self.Y1 + 2.0*h*self.Y2)
        f(self.Y4, t + h, y + h*self.Y3)
        
        error = N.linalg.norm(h/6*(2*self.Y2 + self.Z3 - 2.0*self.Y3 - self.Y4)/scaling) #normalized 
        
        return t+h, y + h/6.0*(self.Y1 + 2.0*self.Y2 + 2.0*self.Y3 + self.Y4), error
    
    def print_statistics(self, verbose):
        """
        Should print the statistics.
        """
        self.log_message('Final Run Statistics: %s \n' % self.problem.name,                  verbose)
        self.log_message(' Number of Steps                : %s '%(self.statistics["nsteps"]),             verbose)
        self.log_message(' Number of Function Evaluations : %s '%(self.statistics["nfcn"]),               verbose)
        self.log_message('\nSolver options:\n',                                              verbose)
        self.log_message(' Solver             : RungeKutta4',                                verbose)
        self.log_message(' Solver type        : Adaptive',                                   verbose)
        self.log_message(' Relative tolerance : ' + str(self.options["rtol"]),        verbose)
        self.log_message(' Absolute tolerance : ' + str(self.options["atol"]) + '\n', verbose)
    
    
class RungeKutta4(Explicit_ODE):
    """
    This solver solves an explicit ordinary differential equation using 
    a Runge-Kutta method of order 4.
    
    We want to approximate the solution to the ordinary differential 
    equation of the form,
    
    .. math::

        \dot{y} = f(t,y), \quad y(t_0) = y_0 .
        
    Using a Runge-Kutta method of order 4, the approximation is defined as 
    follow,
    
    .. math::
    
        y_{n+1} = y_n + \\frac{1}{6}(k_1+2k_2+2k_3+k_4)
        
    where,
    
    .. math::
    
        k_1 = hf(t_n,y_n)
        
        k_2 = hf(t_n+\\frac{1}{2}h,y_n+\\frac{1}{2}k_1)
        
        k_3 = hf(t_n+\\frac{1}{2}h,y_n+\\frac{1}{2}k_2)
        
        k_4 = hf(t_n+h,y_n+k_3)
        
    with :math:`h` being the step-size and :math:`y_n` the previous 
    solution to the equation.
    """
    def __init__(self, problem):
        Explicit_ODE.__init__(self, problem) #Calls the base class
        
        #Solver options
        self.options["h"] = 0.01
        
        #Internal temporary result vector
        self.Y1 = N.array([0.0]*len(self.y0))
        self.Y2 = N.array([0.0]*len(self.y0))
        self.Y3 = N.array([0.0]*len(self.y0))
        self.Y4 = N.array([0.0]*len(self.y0))
        
        #RHS-Function
        self.f = problem.rhs_internal
        
        #Solver support
        self.supports["one_step_mode"] = True
        
    def step(self, t, y, tf, opts):
        initialize = opts["initialize"]
        
        if initialize:
            self.solver_iterator = self._iter(t,y,tf)

        return self.solver_iterator.next()
    
    def integrate(self, t, y, tf, opts):
        """
        Integrates (t,y) values until t > tf
        """
        [flags, tlist, ylist] = zip(*list(self._iter(t, y, tf)))
        
        return flags[-1], tlist, ylist
    
    def _set_h(self,h):
        try:
            self.options["h"] = float(h)
        except:
            raise AssimuloException("Step-size must be a (scalar) float.")
    
    def _get_h(self):
        """
        Defines the step-size that is to be used by the solver.
        
            Parameters::
            
                maxh    
                        - Default '0.01'.
                          
                        - Should be a float.
                        
                            Example:
                                maxh = 0.01
                                
        """
        return self.options["h"]
        
    h=property(_get_h,_set_h)
    
    def _iter(self,t,y,tf):
        h = self.options["h"]
        h = min(h, N.abs(tf-t))
        
        while t+h < tf:
            t, y = self._step(t, y, h)
            yield ID_PY_OK, t,y
            h=min(h, N.abs(tf-t))
        else:
            t, y = self._step(t, y, h)
            yield ID_PY_COMPLETE, t, y
    
    def _step(self, t, y, h):
        """
        This calculates the next step in the integration.
        """
        f = self.f
        
        f(self.Y1, t, y)
        f(self.Y2, t + h/2., y + h*self.Y1/2.)
        f(self.Y3, t + h/2., y + h*self.Y2/2.)
        f(self.Y4, t + h, y + h*self.Y3)
        
        return t+h, y + h/6.*(self.Y1 + 2.*self.Y2 + 2.*self.Y3 + self.Y4)
        
    def print_statistics(self, verbose):
        """
        Should print the statistics.
        """
        self.log_message('Final Run Statistics: %s \n' % self.problem.name,        verbose)
        self.log_message(' Step-length          : %s '%(self.options["h"]), verbose)
        self.log_message('\nSolver options:\n',                                    verbose)
        self.log_message(' Solver            : RungeKutta4',                       verbose)
        self.log_message(' Solver type       : Fixed step\n',                      verbose)

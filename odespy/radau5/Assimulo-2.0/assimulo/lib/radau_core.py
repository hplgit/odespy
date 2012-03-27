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

import pylab as P
import numpy as N

from assimulo.ode import *

class Radau_Exception(Exception):
    pass
    

class Radau_Common(object):
    """
    The common attributes for the Radau solvers.
    """
    def _get_h(self):
        """
        Sets the stepsize.
        """
        return self.__h
    
    def _set_h(self, h):
        """
        Sets the stepsize.
        """
        self.__h = h
        
    h = property(fget=_get_h,fset=_set_h)
    
    def print_statistics(self, verbose=NORMAL):
        """
        Prints the run-time statistics for the problem.
        """
        self.log_message('Final Run Statistics: %s \n' % self.problem.name,        verbose)
        
        self.log_message(' Number of Steps                          : '+str(self.statistics["nsteps"]),          verbose)               
        self.log_message(' Number of Function Evaluations           : '+str(self.statistics["nfcn"]),         verbose)
        self.log_message(' Number of Jacobian Evaluations           : '+ str(self.statistics["njac"]),    verbose)
        self.log_message(' Number of F-Eval During Jac-Eval         : '+ str(self.statistics["njacfcn"]),  verbose)
        self.log_message(' Number of Error Test Failures            : '+ str(self.statistics["errfail"]),       verbose)
        self.log_message(' Number of Newton Iterations              : '+ str(self.statistics["nniter"]),        verbose)
        self.log_message(' Number of Newton Convergence Failures    : '+ str(self.statistics["nniterfail"]),       verbose)
        self.log_message(' Number of LU decompositions              : '+ str(self.statistics["nlu"]),       verbose)

        
        self.log_message('\nSolver options:\n',                                      verbose)
        self.log_message(' Solver                  : Radau5 ' + self._type,          verbose)
        self.log_message(' Tolerances (absolute)   : ' + str(self.options["atol"]),  verbose)
        self.log_message(' Tolerances (relative)   : ' + str(self.options["rtol"]),  verbose)
        self.log_message('',                                                         verbose)
        
    
    def plot_stepsize(self):
        """
        Plots the step-size.
        """
        P.semilogy(N.diff(self.t),drawstyle='steps-post')
        P.title(self.problem.name)
        P.ylabel('Step length')
        P.xlabel('Number of steps')
        P.show()
    
    def _set_newt(self, newt):
        """
        Maximal number of Newton iterations.
        
            Parameters::
            
                newt
                        - Default '7'.
                        
                        - Should be an integer.
                        
                            Example:
                                newt = 10
        """
        try:
            self.options["newt"] = int(newt)
        except (ValueError, TypeError):
            raise Radau_Exception('The newt must be an integer or float.')
		
    def _get_newt(self):
        """
        Maximal number of Newton iterations.
        
            Parameters::
            
                newt
                        - Default '7'.
                        
                        - Should be an integer.
                        
                            Example:
                                newt = 10
        """
        return self.options["newt"]
		
    newt = property(_get_newt,_set_newt)
    
    def _set_fnewt(self, fnewt):
        try:
            self.options["fnewt"] = float(fnewt)
        except (ValueError, TypeError):
            raise Radau_Exception('The fnewt must be an integer or float.')
        
    def _get_fnewt(self):
        """
        Stopping criterion for Newton's method, usually chosen <1.
        Smaller values of fnewt make the code slower, but safer.
        
            Parameters::
            
                fnewt
                        - Default min(0.03,rtol**0.5)
                        
                        - Should be a float.
                        
                            Example:
                                fnewt = 0.05
        """
        return self.options["fnewt"]
        
    fnewt = property(_get_fnewt,_set_fnewt)
    
    def _set_safe(self, safe):
        try:
            self.options["safe"] = float(safe)
        except (ValueError, TypeError):
            raise Radau_Exception('The safe must be an integer or float.')

    def _get_safe(self):
        """
        The safety factor in the step-size prediction.
        
            Parameters::
            
                safe
                        - Default '0.9'
                        
                        - Should be float.
                        
                            Example:
                                safe = 0.8
        """
        return self.options["safe"]
        
    safe = property(_get_safe, _set_safe)
    
    def _set_thet(self, thet):
        try:
            self.options["thet"] = float(thet)
        except (ValueError, TypeError):
            raise Radau_Exception('The thet must be an integer or float.')
        
    def _get_thet(self):
        """
        Value for determine if the Jacobian is to be recomputed or not.
        Increasing thet makes the code compute new Jacobians more seldom.
        Negative thet forces the code to compute the Jacobian after every accepted step.
        
            Parameters::
            
                thet
                        - Default '0.003'
                        
                        - Should be float.
                        
                            Example:
                                thet = 0.01
        """
        return self.options["thet"]
        
    thet = property(_get_thet, _set_thet)
    
    def _set_max_h(self,max_h):
        try:
            self.options["maxh"] = float(max_h)
        except (ValueError,TypeError):
            raise Radau_Exception('Maximal stepsize must be a (scalar) float.')
        if self.options["maxh"] < 0:
            raise Radau_Exception('Maximal stepsize must be a positiv (scalar) float.')
        
    def _get_max_h(self):
        """
        Defines the maximal step-size that is to be used by the solver.
        
            Parameters::
            
                maxh    
                        - Default final time - current time.
                          
                        - Should be a float.
                        
                            Example:
                                maxh = 0.01
                                
        """
        return self.options["maxh"]
        
    maxh=property(_get_max_h,_set_max_h)    
    
    def _set_initial_step(self, initstep):
        try:
            self.options["inith"] = float(initstep)
        except (ValueError, TypeError):
            raise Radau_Exception('The initial step must be an integer or float.')
        
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
    
    
    def _set_quot1(self, quot1):
        try:
            self.options["quot1"] = float(quot1)
        except (ValueError, TypeError):
            raise Radau_Exception('The quot1 must be an integer or float.')
    
    def _get_quot1(self):
        """
        If quot1 < current step-size / old step-size < quot2 the the step-size is not changed.
        This saves LU-decompositions and computing time for large systems.
        
            Parameters::
            
                quot1
                        - Default 1.0
                        
                        - Should be a float.
                        
                            Example:
                                quot1 = 0.9
        """
        return self.options["quot1"]
        
    quot1 = property(_get_quot1, _set_quot1)
    
    def _set_quot2(self, quot2):
        try:
            self.options["quot2"] = float(quot2)
        except (ValueError, TypeError):
            raise Radau_Exception('The quot2 must be an integer or float.')
    
    def _get_quot2(self):
        """
        If quot1 < current step-size / old step-size < quot2 the the step-size is not changed.
        This saves LU-decompositions and computing time for large systems.
        
            Parameters::
            
                quot2
                        - Default 1.2
                        
                        - Should be a float.
                        
                            Example:
                                quot2 = 1.2
        """
        return self.options["quot2"]
        
    quot2 = property(_get_quot2, _set_quot2)
    
    def _set_fac1(self, fac1):
        try:
            self.options["fac1"] = float(fac1)
        except (ValueError, TypeError):
            raise Radau_Exception('The fac1 must be an integer or float.')
            
    def _get_fac1(self):
        """
        Parameters for step-size selection. The new step-size is chosen
        subject to the restriction fac1 <= current step-size / old step-size <= fac2.
        
            Parameters::
            
                fac1
                        - Default 0.2
                        
                        - Should be a float.
                        
                            Example:
                                fac1 = 0.1
        """
        return self.options["fac1"]
        
    fac1 = property(_get_fac1, _set_fac1)
    
    def _set_fac2(self, fac2):
        try:
            self.options["fac2"] = float(fac2)
        except (ValueError, TypeError):
            raise Radau_Exception('The fac2 must be an integer or float.')
        
    def _get_fac2(self):
        """
        Parameters for step-size selection. The new step-size is chosen
        subject to the restriction fac1 <= current step-size / old step-size <= fac2.
        
            Parameters::
            
                fac2
                        - Default 8.0
                        
                        - Should be a float.
                        
                            Example:
                                fac2 = 10.0
        """
        return self.options["fac2"]
        
    fac2 = property(_get_fac2, _set_fac2)
    
    def _set_usejac(self, jac):
        self.options["usejac"] = bool(jac)
    
    def _get_usejac(self):
        """
        This sets the option to use the user defined jacobian. If a
        user provided jacobian is implemented into the problem the
        default setting is to use that jacobian. If not, an
        approximation is used.
        
            Parameters::
            
                usejac  
                        - True - use user defined jacobian
                          False - use an approximation
                    
                        - Should be a boolean.
                        
                            Example:
                                usejac = False
        """
        return self.options["usejac"]
    
    usejac = property(_get_usejac,_set_usejac)
    
    def _set_atol(self,atol):
        
        self.options["atol"] = N.array(atol,dtype=N.float) if len(N.array(atol,dtype=N.float).shape)>0 else N.array([atol],dtype=N.float)
    
        if len(self.options["atol"]) == 1:
            self.options["atol"] = self.options["atol"]*N.ones(self._leny)
        elif len(self.options["atol"]) != self._leny:
            raise Radau_Exception("atol must be of length one or same as the dimension of the problem.")

    def _get_atol(self):
        """
        Defines the absolute tolerance(s) that is to be used by the solver.
        Can be set differently for each variable.
        
            Parameters::
            
                atol    
                        - Default '1.0e-6'.
                
                        - Should be a positive float or a numpy vector
                          of floats.
                        
                            Example:
                                atol = [1.0e-4, 1.0e-6]
        """
        return self.options["atol"]
    
    atol=property(_get_atol,_set_atol)
    
    def _set_rtol(self,rtol):
        try:
            self.options["rtol"] = float(rtol)
        except (ValueError, TypeError):
            raise Radau_Exception('Relative tolerance must be a (scalar) float.')
        if self.options["rtol"] <= 0.0:
            raise Radau_Exception('Relative tolerance must be a positive (scalar) float.')
    
    def _get_rtol(self):
        """
        Defines the relative tolerance that is to be used by the solver.
        
            Parameters::
            
                rtol    
                        - Default '1.0e-6'.
                
                        - Should be a positive float.
                        
                            Example:
                                rtol = 1.0e-4
        """
        return self.options["rtol"]
        
    rtol=property(_get_rtol,_set_rtol)
    
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
            raise Radau_Exception("Maximum number of steps must be a positive integer.")
        self.options["maxsteps"] = max_steps
    
    maxsteps = property(_get_maxsteps, _set_maxsteps)

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
try:
    from lib.sundials_kinsol_core_wSLU import KINSOL_wrap, KINError
except:
    from lib.sundials_kinsol_core import KINSOL_wrap, KINError

from scipy.linalg import pinv2
from scipy.optimize import fminbound
from numpy.linalg import solve
import numpy as N
import pylab as P
import operator as O
import re
import time
from assimulo.problem_algebraic import *

class KINSOL_Exception(Exception):
    pass
    

class KINSOL:
    
    def __init__(self,problem):
        """
        Create the solver
        
        Parameters::
            problem--
                instance of ProblemAlgebraic found in problem_algebraic.py
        """
        
        self.solver = KINSOL_wrap()
        
        # extract info from problem
        self.problem = problem
        if hasattr(self.problem,'_x0'):
            try:
                _x0 = self.problem.get_x0()
            except ProblemAlg_Exception:
                raise KINSOL_Exception("Problem has not implemented method 'get_x0'")
        else:
            raise KINSOL_Exception("Problem has no instance '_x0'")
        
        # calculate dimension
        try:
            if isinstance(_x0, int) or isinstance(_x0, float):
                self.x0 = [_x0]
            else:
                self.x0 = _x0
            self.dim = len([N.array(self.x0, dtype=float)][0])
        except ValueError:
            raise KINSOL_Exception("Initial guess must be a Numpy.array with either ints or floats.")
        
        # check for functions and test them
        try:
            tmp = self.problem.f(self.x0)
            self.norm_of_res = 10000*N.linalg.norm(self.x0)
            #self.norm_of_res = 1e20
            self.func = self.problem.f
        except ProblemAlg_Exception:
            raise KINSOL_Exception("Problem has not implemented method 'f'")
        except IndexError:
            raise KINSOL_Exception("Problem has mismatching dimensions of f and initial guess")
        
        # Calculate 'stupid' scaling
        res = self.problem.f(self.x0)
        self.fscale = N.ones(self.dim)
        
        for i,val in zip(N.arange(self.dim),res):
            if val > 10:
                self.fscale[i] = 1.0/val
        
        # check for constraints and test them
        broken_constraints = []
        if hasattr(self.problem, 'get_constraints'):
            self.constraints = self.problem.get_constraints()
            
            if self.constraints != None:
                # test if constraints are of correct type
                if type(self.constraints).__name__ != 'ndarray':
                    raise KINSOL_Exception("Constraints must be of type numpy.ndarray")
                
                if len(self.constraints) != len(self.x0):
                    raise KINSOL_Exception("Constraints must have same length as x0")
                # Test if initial guess x0 is consistant with constraints
                for c,xi,i in zip(self.constraints,self.x0,N.arange(0,self.x0.__len__())):
                    if re.search('float',type(c).__name__) == None:
                        print "Type problem with: ", c, type(c).__name__
                        raise KINSOL_Exception("Constraints must contain floats.")
                    if abs(c) > 2:
                        raise KINSOL_Exception("Entries in constraint vector must be between -2 and 2, see documentation.")
                    
                    if c != 0.0:
                        if c == 1.0:
                            if xi < 0.0:
                                broken_constraints.append(i)
                        elif c == 2.0:
                            if xi <= 0.0:
                                broken_constraints.append(i)
                        elif c == -1.0:
                            if xi > 0.0:
                                broken_constraints.append(i)
                        elif c == -2.0:
                            if xi >= 0.0:
                                broken_constraints.append(i)
                        else:
                            raise KINSOL_Exception("Constraint vector contains illegal elements.")
                
        else:
            self.constraints = None
            
        if broken_constraints != []:
            print "Variables breaking initial constraint: "
            for i in broken_constraints:
                self.problem.print_var_info(i)

            raise KINSOL_Exception("Initial guess does not fulfill applied constraints.")
        
        
        if hasattr(self.problem, 'check_constraints'):
            self.check_with_model = True
        else:
            self.check_with_model = False 
            
        self._use_jac = True
        self.reg_count = 0
        self.lin_count = 0
        self.verbosity = 0
        self.max_reg = 2.0
        self.use_sparse = False
        self.reg_param = 0.0
        self.exec_time = 0
        self._use_ls = True
        self._use_fscale = False
                
    def set_jac_usage(self,use_jac):
        """
        Set whether to use the jacobian supplied by the model
        or if we are to calculate it numericaly
            
        Parameters::
        
            use_jac --
                Boolean set to True if the jacobian is to be 
                supplied by the model
                
        """
        if type(use_jac).__name__ == 'bool':
            self._use_jac = use_jac
        else:
            raise KINSOL_Exception("The variable sent to 'set_jac_usage' must be a boolean.")
        
    def use_LineSearch(self,use_ls):
        """
        Set whether the solver starts using the Linesearch Algorithm
        or not
        
        Parameters::
            
            use_ls --
                Boolean set to True if LineSearch is to
                be used.
                
        """
        if type(use_ls).__name__ == 'bool':
            self._use_ls = use_ls
        else:
            raise KINSOL_Exception("The variable sent to 'use_LineSearch' must be a boolean.")
        
    def use_fscale(self,use_fscale):
        """
        Set whether the solver starts using the Linesearch Algorithm
        or not
        
        Parameters::
            
            use_ls --
                Boolean set to True if LineSearch is to
                be used.
                
        """
        if type(use_fscale).__name__ == 'bool':
            self._use_fscale = use_fscale
        else:
            raise KINSOL_Exception("The variable sent to 'use_fscale' must be a boolean.")   
        
        
    def set_verbosity(self,verbosity):
        """
        Method used to set the verbosity of KINSOL
        
        Parameters::
            
            verbosity --
                Integer set to one of 0,1,2,3 where 0 is no info and 3 the
                most info. For more information please see the documentation for KINSOL.
        """
        type_name = type(verbosity).__name__
        if re.search('int',type_name) != None:
            
            # It is an integer, tes bounds
            if verbosity < 4 and verbosity > -1:
                self.verbosity = verbosity
            else:
                raise KINSOL_Exception("The variable sent to 'set_verbosity' must be either 0, 1, 2 or 3.")
        else:
            raise KINSOL_Exception("The variable sent to 'set_verbosity' must be an integer.")
        
    def set_sparsity(self,use_sparse):
        """
        Method used to set if the problem should be treated as sparse by
        the linear solver in KINSOL. If the problem supplied has not implemented
        a method sparse_jac an exception will be thrown
        
        Parameters::
        
            use_sparse --
                Boolean set to True if the problem is to be treated as
                sparse and False otherwise.
                
        """
        
        if hasattr(self.problem,'sparse_jac'):
            self.use_sparse = use_sparse
        else:
            raise KINSOL_Exception("The problem must have implemented a method 'sparse_jac' for sparsity to by used.")
        
    def set_reg_param(self,reg_param):
        """
        Method used to set the regularization parameter.
        If set to zero the parameter will be set by the solver.
        
        Parameters::
        
            reg_param --
                Float larger or equal to zero.
                
        """
        
        type_name = type(reg_param).__name__
        if re.search('float',type_name) != None:
            if reg_param < 0:
                raise KINSOL_Exception("Value sent to set_reg_param must be equal to or larger than zero.")
            else:
                self.reg_param = reg_param
        else:
            raise KINSOL_Exception("The variable sent to 'set_reg_param' must be a float.")
        
    def solve(self):
        """
        Function called when solving function rhs_fct
        
        """
        # check for jacobian and set it if present and to be used
        if self.use_sparse:
            if self._use_jac and hasattr(self.problem,'sparse_jac'):
                jac = self.problem.sparse_jac
            else:
                jac = None
        else:
            if self._use_jac and hasattr(self.problem,'jac'):
                jac = self.problem.jac
            else:
                jac = None
            
        # Initialize solver and solve        
        
        solved = False
        local_min = False

        res = N.zeros(self.x0.__len__())
        while (not solved) and self.reg_count < 2:
            try:
                if self._use_fscale:
                    self.solver.KINSOL_init(self.func,self.x0,self.dim,jac,self.constraints,self.use_sparse,self.verbosity,self.norm_of_res,self.reg_param,self.fscale)
                else:
                    self.solver.KINSOL_init(self.func,self.x0,self.dim,jac,self.constraints,self.use_sparse,self.verbosity,self.norm_of_res,self.reg_param,None)
                start = time.clock()
                res = self.solver.KINSOL_solve(not self._use_ls)
                stop = time.clock()
                self.exec_time += (stop - start)
                solved = True
            except KINError as error:
                if error.value == 42:
                    # Try the heuristic
                    if hasattr(self.problem, 'get_heuristic_x0'):
                        print "----------------------------------------------------"
                        print "      Solver stuck with zero step-length."
                        print "----------------------------------------------------"
                        print "The following variables have start value zero"
                        print "and min set to zero causing the zero step-lenght."
                        print "These settings are either set by default or by user."
                        print ""

                        self.x0 = self.problem.get_heuristic_x0()
                        self.reg_count += 1
                        
                        print ""
                        print "This setting (start and min to zero) can often"
                        print "cause problem when initializing the system. "
                        print ""
                        print "To avoid this the above variables have"
                        print "their start attributes reset to one."
                        print ""
                        print "Trying to solve the system again..."
                    else:
                        raise KINSOL_Exception("Regularization failed due to constraints, tried getting heuristic initial guess but failed.")
                

                elif (error.value == 2):
                    print "---------------------------------------------------------"
                    print ""
                    print " !!! WARNING !!!"
                    print ""
                    print " KINSOL has returned a result but the algorithm has converged"
                    print " to a local minima, the initial values are NOT consistant!"
                    print ""
                    print "---------------------------------------------------------"
                    solved = True
                    local_min = True
                else:
                    # Other error, send onward as exception
                    self.problem.check_constraints(res)
                    raise KINSOL_Exception(error.msg[error.value])
        
        if not solved:
            self.solver.Free_KINSOL()
            raise KINSOL_Exception("Algorithm exited solution loop without finding a solution, please contact Assimulo support.")

        if self.check_with_model:
            self.problem.check_constraints(res)
        if not local_min:
            print "Problem sent to KINSOL solved."
            
        return res
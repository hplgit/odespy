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
import scipy as S
import scipy.linalg as LIN

from assimulo.exception import *
from assimulo.ode import *

from assimulo.explicit_ode import Explicit_ODE
from assimulo.implicit_ode import Implicit_ODE
from assimulo.lib.radau_core import Radau_Common

class Radau5ODE(Radau_Common,Explicit_ODE):
    """
    Radau IIA fifth-order three-stages with step-size control and continuous output.
    Based on the FORTRAN code by E.Hairer and G.Wanner, which can be found here: 
    http://www.unige.ch/~hairer/software.html
    
    Details about the implementation (FORTRAN) can be found in the book,::
    
        Solving Ordinary Differential Equations II,
        Stiff and Differential-Algebraic Problems
        
        Authors: E. Hairer and G. Wanner
        Springer-Verlag, ISBN: 3-540-60452-9
    
    This code is aimed at providing a Python implementation of the original code.
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
        
        #Default values
        self.options["inith"] = 0.01
        self.options["newt"]     = 7 #Maximum number of newton iterations
        self.options["thet"]     = 1.e-3 #Boundary for re-calculation of jac
        self.options["fnewt"]    = 0 #Stopping critera for Newtons Method
        self.options["quot1"]    = 1.0 #Parameters for changing step-size (lower bound)
        self.options["quot2"]    = 1.2 #Parameters for changing step-size (upper bound)
        self.options["fac1"]     = 0.2 #Parameters for step-size selection (lower bound)
        self.options["fac2"]     = 8.0 #Parameters for step-size selection (upper bound)
        self.options["maxh"]     = N.inf #Maximum step-size.
        self.options["safe"]     = 0.9 #Safety factor
        self.options["atol"]     = 1.0e-6 #Absolute tolerance
        self.options["rtol"]     = 1.0e-6 #Relative tolerance
        self.options["usejac"]   = True if self.problem_info["jac_fcn"] else False
        self.options["maxsteps"] = 10000
        
        # - Statistic values
        self.statistics["nsteps"] = 0 #Number of steps
        self.statistics["nfcn"] = 0 #Number of function evaluations
        self.statistics["njac"] = 0 #Number of jacobian evaluations
        self.statistics["njacfcn"] = 0 #Number of function evaluations when evaluating the jacobian
        self.statistics["nniter"] = 0 #Number of nonlinear iterations
        self.statistics["nniterfail"] = 0 #Number of nonlinear failures
        self.statistics["errfail"] = 0 #Number of step rejections
        self.statistics["nlu"] = 0 #Number of LU decompositions
        
        #Internal values
        self._curjac = False #Current jacobian?
        self._itfail = False #Iteration failed?
        self._needjac = True #Need to update the jacobian?
        self._needLU = True #Need new LU-factorisation?
        self._first = True #First step?
        self._rejected = True #Is the last step rejected?
        self._leny = len(self.y) #Dimension of the problem
        self._oldh = 0.0 #Old stepsize
        self._olderr = 1.0 #Old error
        self._eps = N.finfo('double').eps
        self._col_poly = N.zeros(self._leny*3)
        self._type = '(explicit)'
        self._curiter = 0 #Number of current iterations
        
        #RHS-Function
        self.f = problem.rhs_internal
        
        #Internal temporary result vector
        self.Y1 = N.array([0.0]*len(self.y0))
        self.Y2 = N.array([0.0]*len(self.y0))
        self.Y3 = N.array([0.0]*len(self.y0))
        self._f0 = N.array([0.0]*len(self.y0))
        
        #Solver support
        self.supports["one_step_mode"] = True
        self.supports["interpolated_output"] = True
        
        # - Retrieve the Radau5 parameters
        self._load_parameters() #Set the Radau5 parameters
    
    def initialize(self):
        #Reset statistics
        for k in self.statistics.keys():
            self.statistics[k] = 0
    
    def step_generator(self, t, y, tf, opts):
        
        if opts["initialize"]:
            self._oldh = self.inith
            self.h = self.inith
            self._fac_con = 1.0
        
        if self.fnewt == 0:
            self.fnewt = max(10.*self._eps/self.rtol,min(0.03,self.rtol**0.5))

        self.f(self._f0,t,y)
        self.statistics["nfcn"] +=1
        self._tc = t
        self._yc = y
        
        for i in xrange(self.maxsteps):
            
            if t < tf:
                t, y = self._step(t, y)
                self._tc = t
                self._yc = y
                
                if self.h > N.abs(tf-t):
                    self.h = N.abs(tf-t)
                
                if t < tf:
                    yield ID_PY_OK, t, y
                else:
                    yield ID_PY_COMPLETE, t, y
                    break

                self._first = False 
        else:
            raise Explicit_ODE_Exception('Final time not reached within maximum number of steps')
        
        #t, y = self._step(t,y)
        #yield ID_PY_COMPLETE, t, y
    
    def step(self, t, y, tf, opts):
        if opts["initialize"]:
            self._next_step = self.step_generator(t,y,tf,opts)
        return self._next_step.next()
    
    def integrate(self, t, y, tf, opts):
        
        if opts["output_list"] != None:
            
            output_list = opts["output_list"]
            output_index = opts["output_index"]
            
            next_step = self.step_generator(t,y,tf,opts)
            
            tlist,ylist = [], []
            res = [ID_PY_OK]
            
            while res[0] != ID_PY_COMPLETE:
                res = next_step.next()
                try:
                    while output_list[output_index] <= res[1]:
                        tlist.append(output_list[output_index])
                        ylist.append(self.interpolate(output_list[output_index]))

                        output_index = output_index+1
                except IndexError:
                    pass
            return res[0], tlist, ylist
        else:
            [flags, tlist, ylist] = zip(*list(self.step_generator(t, y, tf,opts)))

            return flags[-1], tlist, ylist
        
    def _step(self, t, y):
        """
        This calculates the next step in the integration.
        """
        self._scaling = N.array(abs(y)*self.rtol + self.atol) #The scaling used.
        
        while True: #Loop for integrating one step.
            
            self.newton(t,y)
            self._err = self.estimate_error()
            
            if self._err > 1.0: #Step was rejected.
                self._rejected = True
                self.statistics["errfail"] += 1
                ho = self.h
                self.h = self.adjust_stepsize(self._err)
                
                self.log_message('Rejecting step at ' + str(t) + 'with old stepsize' + str(ho) + 'and new ' + str(self.h), SCREAM)
                
                if self._curjac or self._curiter == 1:
                    self._needjac = False
                    self._needLU = True
                else:
                    self._needjac = True
                    self._needLU = True
            else:
                self.log_message('Accepting step at ' + str(t) + 'with stepsize ' + str(self.h),SCREAM)
                
                self.statistics["nsteps"] += 1
                
                tn = t+self.h #Preform the step
                yn = y+self._Z[2*self._leny:3*self._leny]
                self.f(self._f0,tn,yn)
                self.statistics["nfcn"] += 1
                
                self._oldoldh = self._oldh #Store the old(old) step-size for use in the test below.
                self._oldh = self.h #Store the old step-size
                self._oldt = t #Store the old time-point
                self._newt = tn #Store the new time-point
                
                #Adjust the new step-size
                ht = self.adjust_stepsize(self._err, predict=True)
                self.h = min(self.h,ht) if self._rejected else ht
                
                self._rejected = False
                self._curjac = False
                
                if self._oldoldh == self.h and (self._theta <= self.thet):# or self._curiter==1):
                    self._needjac = False
                    self._needLU = False
                else:
                    if self._theta <= self.thet: #or self._curiter == 1:
                        self._needjac = False
                        self._needLU = True
                    else:
                        self._needjac = True
                        self._needLU = True
                if self.thet < 0:
                    self._needjac = True
                    self._needLU = True
                        
                self._olderr = max(self._err,1.e-2) #Store the old error
                break
                
        self._col_poly = self._collocation_pol(self._Z, self._col_poly, self._leny) #Calculate the new collocation polynomial
        
        return tn, yn #Return the step
    
    def _collocation_pol(self, Z, col_poly, leny):
        
        col_poly[2*leny:3*leny] = Z[:leny] / self.C[0,0]
        col_poly[leny:2*leny]   = ( Z[:leny] - Z[leny:2*leny] ) / (self.C[0,0]-self.C[1,0])
        col_poly[:leny]         = ( Z[leny:2*leny] -Z[2*leny:3*leny] ) / (self.C[1,0]-1.)
        col_poly[2*leny:3*leny] = ( col_poly[leny:2*leny] - col_poly[2*leny:3*leny] ) / self.C[1,0]
        col_poly[leny:2*leny]   = ( col_poly[leny:2*leny] - col_poly[:leny] ) / (self.C[0,0]-1.)
        col_poly[2*leny:3*leny] =   col_poly[leny:2*leny]-col_poly[2*leny:3*leny]
        
        return col_poly
    
    def _radau_F(self, Z, t, y):
        
        Z1 = Z[:self._leny]
        Z2 = Z[self._leny:2*self._leny]
        Z3 = Z[2*self._leny:3*self._leny]

        self.f(self.Y1,t+self.C[0]*self.h, y+Z1)
        self.f(self.Y2,t+self.C[1]*self.h, y+Z2)
        self.f(self.Y3,t+self.C[2]*self.h, y+Z3)
        
        self.statistics["nfcn"] += 3
        
        return N.hstack((N.hstack((self.Y1,self.Y2)),self.Y3))
    
    def calc_start_values(self):
        """
        Calculate newton starting values.
        """
        if self._first:
            Z = N.zeros(self._leny*3)
            W = N.zeros(self._leny*3)
        else:
            Z = self._Z
            cq = self.C*self.h/self._oldh#self._oldoldh#self._oldh
            newtval = self._col_poly
            leny = self._leny
            
            Z[:leny]        = cq[0,0]*(newtval[:leny]+(cq[0,0]-self.C[1,0]+1.)*(newtval[leny:2*leny]+(cq[0,0]-self.C[0,0]+1.)*newtval[2*leny:3*leny]))
            Z[leny:2*leny]  = cq[1,0]*(newtval[:leny]+(cq[1,0]-self.C[1,0]+1.)*(newtval[leny:2*leny]+(cq[1,0]-self.C[0,0]+1.)*newtval[2*leny:3*leny]))
            Z[2*leny:3*leny]= cq[2,0]*(newtval[:leny]+(cq[2,0]-self.C[1,0]+1.)*(newtval[leny:2*leny]+(cq[2,0]-self.C[0,0]+1.)*newtval[2*leny:3*leny]))
            
            W = N.dot(self.T2,Z)
            
        return Z, W
    
    def newton(self,t,y):
        """
        The newton iteration. 
        """
        
        for k in xrange(20):
            
            self._curiter = 0 #Reset the iteration
            self._fac_con = max(self._fac_con, self._eps)**0.8;
            self._theta = abs(self.thet);
            
            if self._needjac:
                self._jac = self.jacobian(t,y)
            
            if self._needLU:
                self.statistics["nlu"] += 1
                self._a = self._alpha/self.h
                self._b = self._beta/self.h
                self._g = self._gamma/self.h
                self._B = self._g*self.I - self._jac
                
                self._P1,self._L1,self._U1 = S.linalg.lu(self._B) #LU decomposition
                self._P2,self._L2,self._U2 = S.linalg.lu(self._a*self.I-self._jac)
                self._P3,self._L3,self._U3 = S.linalg.lu(self._b*self.I-self._jac)
                
                self._needLU = False
                
                if min(abs(N.diag(self._U1)))<self._eps:
                    raise Explicit_ODE_Exception('Error, gI-J is singular.')
                    
            Z, W = self.calc_start_values()
        
            for i in xrange(self.newt):
                self._curiter += 1 #The current iteration
                self.statistics["nniter"] += 1 #Adding one iteration
                
                #Solve the system
                Z = N.dot(self.T2,self._radau_F(Z.real,t,y))

                Z[:self._leny]              =Z[:self._leny]              -self._g*N.dot(self.I,W[:self._leny])
                Z[self._leny:2*self._leny]  =Z[self._leny:2*self._leny]  -self._a*N.dot(self.I,W[self._leny:2*self._leny])   #+self._b*N.dot(self.I,W[2*self._leny:3*self._leny])
                Z[2*self._leny:3*self._leny]=Z[2*self._leny:3*self._leny]-self._b*N.dot(self.I,W[2*self._leny:3*self._leny]) #-self._a*N.dot(self.I,W[2*self._leny:3*self._leny])
                
                Z[:self._leny]              =N.linalg.solve(self._U1,N.linalg.solve(self._L1,N.linalg.solve(self._P1,Z[:self._leny])))
                Z[self._leny:2*self._leny]  =N.linalg.solve(self._U2,N.linalg.solve(self._L2,N.linalg.solve(self._P2,Z[self._leny:2*self._leny])))
                Z[2*self._leny:3*self._leny]=N.linalg.solve(self._U3,N.linalg.solve(self._L3,N.linalg.solve(self._P3,Z[2*self._leny:3*self._leny])))
                #----
                newnrm = N.linalg.norm(Z.reshape(-1,self._leny)/self._scaling,'fro')/N.sqrt(3.*self._leny)
                      
                if i > 0:
                    thq = newnrm/oldnrm
                    if i == 1:
                        self._theta = thq
                    else:
                        self._theta = N.sqrt(thq*thqold)
                    thqold = thq
                    
                    if self._theta < 0.99: #Convergence
                        self._fac_con = self._theta/(1.-self._theta)
                        dyth = self._fac_con*newnrm*self._theta**(self.newt-(i+1)-1)/self.fnewt
                        
                        if dyth >= 1.0: #Too slow convergence
                            qnewt = max(1.e-4,min(20.,dyth))
                            self.h = 0.8*qnewt**(-1.0/(4.0+self.newt-(i+1)-1))*self.h
                            self._itfail = True
                            self._rejected = True
                            break
                    else: #Not convergence, abort
                        self._itfail = True
                        break
                
                oldnrm = max(newnrm,self._eps) #Store oldnorm
                W = W+Z #Perform the iteration

                Z = N.dot(self.T3,W) #Calculate the new Z values
                
                if self._fac_con*newnrm <= self.fnewt: #Convergence?
                    self._itfail = False;
                    break
                
            else: #Iteration failed
                self._itfail = True
                
            if not self._itfail: #Newton iteration converged
                self._Z = Z.real
                break
            else: #Iteration failed
                self.log_message('Iteration failed at time %e with step-size %e'%(t,self.h),SCREAM)
                
                self.statistics["nniterfail"] += 1
                self._rejected = True #The step is rejected
                
                if self._theta >= 0.99:
                    self.h = self.h/2.0
                if self._curjac:
                    self._needjac = False
                    self._needLU = True
                else:
                    self._needjac = True
                    self._needLU = True
        else:
            raise Explicit_ODE_Exception('Newton iteration failed at time %e with step-size %e'%(t,self.h))
        
    def adjust_stepsize(self, err, predict=False):
        
        fac = min(self.safe, self.safe*(2.*self.newt+1.)/(2.*self.newt+self._curiter))
        quot = max(1./self.fac2,min(1./self.fac1,(err**0.25)/fac))        
        hnormal = self.h/quot
        
        if predict:
            if not self._first:
                facgus = (self._hacc/self.h)*(err**2/self._olderr)**0.25/self.safe
                facgus = max(1./self.fac2,min(1./self.fac1,facgus))
                quot = max(quot,facgus)
                h = self.h/quot
            else:
                h = hnormal
            self._hacc = self.h
        else:
            h = hnormal
        
        qt = h/self.h
        
        if (qt >= self.quot1) and (qt <= self.quot2):
            h = self.h
            
        if self._first and err>=1.0:
            h = self.h/10.
        
        if h < self._eps:
            raise Explicit_ODE_Exception('Step-size to small at %e with h = %e'%(self._tc,self.h))
        
        if h > self.maxh:
            h = self.maxh
        
        return h
        
    def estimate_error(self):
        
        temp = 1./self.h*(self.E[0]*self._Z[:self._leny]+self.E[1]*self._Z[self._leny:2*self._leny]+self.E[2]*self._Z[2*self._leny:3*self._leny])

        scal = self._scaling#/self.h
        err_v = N.linalg.solve(self._U1,N.linalg.solve(self._L1,N.linalg.solve(self._P1,self._f0+temp)))
        err = N.linalg.norm(err_v/scal)
        err = max(err/N.sqrt(self._leny),1.e-10)

        if (self._rejected or self._first) and err >= 1.: #If the step was rejected, use the more expensive error estimation
            self.statistics["nfcn"] += 1
            err_new = N.array([0.0]*self._leny)
            self.f(err_new,self._tc,self._yc+err_v)
            err_v =  N.linalg.solve(self._U1,N.linalg.solve(self._L1,N.linalg.solve(self._P1,err_new+temp)))
            err = N.linalg.norm(err_v/scal)
            err = max(err/N.sqrt(self._leny),1.e-10)

        return err
    
    def jacobian(self, t, y):
        """
        Calculates the Jacobian, either by an approximation or by the user
        defined (jac specified in the problem class).
        """
        self._curjac = True #The jacobian is up to date
        self._needLU = True #A new LU-decomposition is needed
        self._needjac = False #A new jacobian is not needed
        
        if self.usejac: #Retrieve the user-defined jacobian
            cjac = self.problem.jac(t,y)
        else:           #Calculate a numeric jacobian
            delt = N.array([(self._eps*max(abs(yi),1.e-5))**0.5 for yi in y])*N.identity(self._leny) #Calculate a disturbance
            Fdelt = N.array([self.problem.rhs(t,y+e) for e in delt]) #Add the disturbance (row by row) 
            grad = ((Fdelt-self.problem.rhs(t,y)).T/delt.diagonal()).T
            cjac = N.array(grad).T

            self.statistics["njacfcn"] += 1+self._leny #Add the number of function evaluations
        
        self.statistics["njac"] += 1 #add the number of jacobian evaluation
        return cjac
    
    def interpolate(self, t, k=0):
        """
        Calculates the continuous output from Radau5.
        """
        leny = self._leny
        s = (t-self._newt)/self._oldh
        Z = self._col_poly
        
        yout = self._yc+s*(Z[:leny]+(s-self.C[1,0]+1.)*(Z[leny:2*leny]+(s-self.C[0,0]+1.)*Z[2*leny:3*leny]))
        return yout
    
    def _load_parameters(self):
        
        #Parameters
        A = N.zeros([3,3])
        A[0,0] = (88.-7.*N.sqrt(6.))/360.0
        A[0,1] = (296.-169.*N.sqrt(6.))/1800.0
        A[0,2] = (-2.0+3.0*N.sqrt(6.))/225.0
        A[1,0] = (296.0+169.0*N.sqrt(6.))/1800.0
        A[1,1] = (88.+7.*N.sqrt(6.))/360.0
        A[1,2] = (-2.-3.*N.sqrt(6.))/225.0
        A[2,0] = (16.0-N.sqrt(6.))/36.0
        A[2,1] = (16.0+N.sqrt(6.))/36.0
        A[2,2] = (1.0/9.0)
        
        C = N.zeros([3,1])
        C[0,0]=(4.0-N.sqrt(6.0))/10.0
        C[1,0]=(4.0+N.sqrt(6.0))/10.0
        C[2,0]=1.0
        
        B = N.zeros([1,3])
        B[0,0]=(16.0-N.sqrt(6.0))/36.0
        B[0,1]=(16.0+N.sqrt(6.0))/36.0
        B[0,2]=1.0/9.0
        
        E = N.zeros(3)
        E[0] = -13.0-7.*N.sqrt(6.)
        E[1] = -13.0+7.0*N.sqrt(6.)
        E[2] = -1.0
        E = 1.0/3.0*E
        
        Ainv = N.linalg.inv(A)
        [eig, T] = N.linalg.eig(Ainv)
        eig = N.array([eig[2],eig[0],eig[1]])
        J = N.diag(eig)

        self._alpha = eig[1]
        self._beta  = eig[2]
        self._gamma = eig[0].real
        
        temp0 = T[:,0].copy()
        temp1 = T[:,1].copy()
        temp2 = T[:,2].copy()
        T[:,0] = temp2
        T[:,1] = temp0
        T[:,2] = temp1
        Tinv = N.linalg.inv(T)
        
        I = N.eye(self._leny)
        I3 = N.eye(3)
        T1 = N.kron(J,I)
        T2 = N.kron(Tinv,I)
        T3 = N.kron(T,I)
        
        self.A = A
        self.B = B
        self.C = C
        self.I = I
        self.E = E
        self.T1 = T1
        self.T2 = T2
        self.T3 = T3
        self.I3 = I3
        self.EIG = eig



class Radau5DAE(Radau_Common,Implicit_ODE):
    """
    Radau IIA fifth-order three-stages with step-size control and continuous output.
    Based on the FORTRAN code by E.Hairer and G.Wanner, which can be found here: 
    http://www.unige.ch/~hairer/software.html
    
    Details about the implementation (FORTRAN) can be found in the book,::
    
        Solving Ordinary Differential Equations II,
        Stiff and Differential-Algebraic Problems
        
        Authors: E. Hairer and G. Wanner
        Springer-Verlag, ISBN: 3-540-60452-9
    
    This code is aimed at providing a Python implementation of the original code.
    """
    def __init__(self, problem):
        """
        Initiates the solver.
        
            Parameters::
            
                problem     
                            - The problem to be solved. Should be an instance
                              of the 'Implicit_Problem' class.
        """
        Implicit_ODE.__init__(self, problem) #Calls the base class
        
        #Internal values
        self._leny = len(self.y) #Dimension of the problem
        self._2leny = 2*self._leny
        
        #Default values
        self.options["inith"] = 0.01
        self.options["newt"]     = 7 #Maximum number of newton iterations
        self.options["thet"]     = 1.e-3 #Boundary for re-calculation of jac
        self.options["fnewt"]    = 0 #Stopping critera for Newtons Method
        self.options["quot1"]    = 1.0 #Parameters for changing step-size (lower bound)
        self.options["quot2"]    = 1.2 #Parameters for changing step-size (upper bound)
        self.options["fac1"]     = 0.2 #Parameters for step-size selection (lower bound)
        self.options["fac2"]     = 8.0 #Parameters for step-size selection (upper bound)
        self.options["maxh"]     = N.inf #Maximum step-size.
        self.options["safe"]     = 0.9 #Safety factor
        self.options["atol"]     = N.array([1.0e-6]*self._leny) #Absolute tolerance
        self.options["rtol"]     = 1.0e-6 #Relative tolerance
        self.options["index"]    = N.array([1]*self._leny+[2]*self._leny)
        self.options["usejac"]   = True if self.problem_info["jac_fcn"] else False
        self.options["maxsteps"] = 10000
        
        # - Statistic values
        self.statistics["nsteps"] = 0 #Number of steps
        self.statistics["nfcn"] = 0 #Number of function evaluations
        self.statistics["njac"] = 0 #Number of jacobian evaluations
        self.statistics["njacfcn"] = 0 #Number of function evaluations when evaluating the jacobian
        self.statistics["nniter"] = 0 #Number of nonlinear iterations
        self.statistics["nniterfail"] = 0 #Number of nonlinear failures
        self.statistics["errfail"] = 0 #Number of step rejections
        self.statistics["nlu"] = 0 #Number of LU decompositions
        
        #Internal values
        self._curjac = False #Current jacobian?
        self._itfail = False #Iteration failed?
        self._needjac = True #Need to update the jacobian?
        self._needLU = True #Need new LU-factorisation?
        self._first = True #First step?
        self._rejected = True #Is the last step rejected?
        self._oldh = 0.0 #Old stepsize
        self._olderr = 1.0 #Old error
        self._eps = N.finfo('double').eps
        self._col_poly = N.zeros(self._2leny*3)
        self._type = '(implicit)'
        self._curiter = 0 #Number of current iterations
        
        #RES-Function
        self.f = problem.res_internal
        self.RES =  N.array([0.0]*len(self.y0))
        
        #Internal temporary result vector
        self.Y1 = N.array([0.0]*len(self.y0))
        self.Y2 = N.array([0.0]*len(self.y0))
        self.Y3 = N.array([0.0]*len(self.y0))
        self._f0 = N.array([0.0]*len(self.y0))
        
        
        #Solver support
        self.supports["one_step_mode"] = True
        self.supports["interpolated_output"] = True
        
        # - Retrieve the Radau5 parameters
        self._load_parameters() #Set the Radau5 parameters
    
    def _set_index(self, index):
        """
        Sets the index of the variables in the problem which in turn
        determine the error estimations.
        
            Parameters::
            
                    index - A list of integers, indicating the index
                            (1,2,3) of the variable.
                            
                            Example:
                                Radau5.index = [2,1]
                            
        """
        if len(index) == self._2leny:
            ind = N.array(index)
        elif len(index) == self._leny:
            ind = N.array(index+(N.array(index)+1).tolist())
        else:
            raise Implicit_ODE_Exception('Wrong number of variables in the index vector.')
        self.options["index"] = ind
            
    def _get_index(self):
        """
        Sets the index of the variables in the problem which in turn
        determine the error estimations.
        
            Parameters::
            
                    index - A list of integers, indicating the index
                            (1,2,3) of the variable.
                            
                            Example:
                                Radau5.index = [2,1]
                            
        """
        return self.options["index"]
        
    index = property(_get_index,_set_index)
    
    def initialize(self):
        #Reset statistics
        for k in self.statistics.keys():
            self.statistics[k] = 0
    
    def step_generator(self, t, y, yd, tf, opts):
        
        if opts["initialize"]:
            self._oldh = self.inith
            self.h = self.inith
            self._fac_con = 1.0
        
        if self.fnewt == 0:
            self.fnewt = max(10.*self._eps/self.rtol,min(0.03,self.rtol**0.5))
            
        self._f0 = self._ode_f(t,N.append(y,yd))
        self.statistics["nfcn"] +=1
        self._tc = t
        self._yc = y
        self._ydc = yd
        
        for i in xrange(self.maxsteps):
            
            if t < tf:
                t, y, yd = self._step(t, y, yd)
                self._tc = t
                self._yc = y
                self._ydc = yd
                
                if self.h > N.abs(tf-t):
                    self.h = N.abs(tf-t)
                
                if t < tf:
                    yield ID_PY_OK, t,y,yd
                else:
                    yield ID_PY_COMPLETE, t, y, yd
                    break

                self._first = False 
        else:
            raise Implicit_ODE_Exception('Final time not reached within maximum number of steps')
    
    def step(self, t, y, yd, tf, opts):
        
        if opts["initialize"]:
            self._next_step = self.step_generator(t,y,yd,tf,opts)
        return self._next_step.next()
    
    def integrate(self, t, y, yd, tf, opts):
        
        if opts["output_list"] != None:
            
            output_list = opts["output_list"]
            output_index = opts["output_index"]
            
            next_step = self.step_generator(t,y,yd,tf,opts)
            
            tlist,ylist,ydlist = [], [], []
            res = [ID_PY_OK]
            
            while res[0] != ID_PY_COMPLETE:
                res = next_step.next()
                try:
                    while output_list[output_index] <= res[1]:
                        tlist.append(output_list[output_index])
                        ylist.append(self.interpolate(output_list[output_index]))
                        ydlist.append(self.interpolate(output_list[output_index],k=1))

                        output_index = output_index+1
                except IndexError:
                    pass
            return res[0], tlist, ylist, ydlist
        else:
            [flags, tlist, ylist, ydlist] = zip(*list(self.step_generator(t, y, yd, tf,opts)))
            
            return flags[-1], tlist, ylist, ydlist
    
    def _ode_f(self, t, y):
        
        #self.res_fcn(t,y[:self._leny],y[self._leny:])
        #return N.hstack((y[self._leny:],self.res_fcn(t,y[:self._leny],y[self._leny:])))
        
        self.f(self.RES,t,y[:self._leny],y[self._leny:])
        return N.hstack((y[self._leny:],self.RES))
    
    def _radau_F(self, Z, t, y, yd):
        
        Z1 = Z[:self._2leny]
        Z2 = Z[self._2leny:2*self._2leny]
        Z3 = Z[2*self._2leny:3*self._2leny]
        
        q = N.append(y,yd)
        
        sol1 = self._ode_f(t+self.C[0]*self.h, q+Z1)
        sol2 = self._ode_f(t+self.C[1]*self.h, q+Z2)
        sol3 = self._ode_f(t+self.C[2]*self.h, q+Z3)
        
        self.statistics["nfcn"] += 3
        
        return N.hstack((N.hstack((sol1,sol2)),sol3))
    
    def _step(self, t, y, yd):
        """
        This calculates the next step in the integration.
        """
        self._scaling = N.array(abs(N.append(y,yd))*self.rtol + self.atol.tolist()*2) #The scaling used.
        
        while True: #Loop for integrating one step.
            
            self.newton(t,y,yd)
            self._err = self.estimate_error()
            
            if self._err > 1.0: #Step was rejected.
                self._rejected = True
                self.statistics["errfail"] += 1
                ho = self.h
                self.h = self.adjust_stepsize(self._err)
                
                self.log_message('Rejecting step at ' + str(t) + 'with old stepsize' + str(ho) + 'and new ' +
                                   str(self.h) + '. Error: ' + str(self._err),SCREAM)
                
                if self._curjac or self._curiter == 1:
                    self._needjac = False
                    self._needLU = True
                else:
                    self._needjac = True
                    self._needLU = True
            else:
                
                self.log_message("Accepting step at " + str(t) + ' with stepsize ' + str(self.h) + '. Error: ' + str(self._err),SCREAM)
                self.statistics["nsteps"] += 1
                
                tn = t+self.h #Preform the step
                yn = y+self._Z[2*self._2leny:3*self._2leny][:self._leny]
                ydn = yd+self._Z[2*self._2leny:3*self._2leny][self._leny:]
                self._f0 = self._ode_f(t,N.append(yn,ydn))
                self.statistics["nfcn"] += 1
                
                self._oldoldh = self._oldh #Store the old(old) step-size for use in the test below.
                self._oldh = self.h #Store the old step-size
                self._oldt = t #Store the old time-point
                self._newt = tn #Store the new time-point
                
                #Adjust the new step-size
                ht = self.adjust_stepsize(self._err, predict=True)
                self.h = min(self.h,ht) if self._rejected else ht
                
                self._rejected = False
                self._curjac = False
                
                if self._oldoldh == self.h and (self._theta <= self.thet or self._curiter==1):
                    self._needjac = False
                    self._needLU = False
                else:
                    if self._theta <= self.thet or self._curiter == 1:
                        self._needjac = False
                        self._needLU = True
                    else:
                        self._needjac = True
                        self._needLU = True
                if self.thet < 0:
                    self._needjac = True
                    self._needLU = True
                        
                self._olderr = max(self._err,1.e-2) #Store the old error
                break
                
        self._col_poly = self._collocation_pol(self._Z, self._col_poly, self._2leny) #Calculate the new collocation polynomial
        
        return tn, yn, ydn #Return the step
    
    def newton(self,t,y,yd):
        """
        The newton iteration. 
        """
        
        for k in xrange(20):
            
            self._curiter = 0 #Reset the iteration
            self._fac_con = max(self._fac_con, self._eps)**0.8;
            self._theta = abs(self.thet);
            
            if self._needjac:
                self._jac = self.jacobian(t,y,yd)
            
            if self._needLU:
                self.statistics["nlu"] += 1
                self._a = self._alpha/self.h
                self._b = self._beta/self.h
                self._g = self._gamma/self.h
                self._B = self._g*self.M - self._jac
                
                self._P1,self._L1,self._U1 = S.linalg.lu(self._B) #LU decomposition
                self._P2,self._L2,self._U2 = S.linalg.lu(self._a*self.M-self._jac)
                self._P3,self._L3,self._U3 = S.linalg.lu(self._b*self.M-self._jac)
                
                self._needLU = False
                
                if min(abs(N.diag(self._U1)))<self._eps:
                    raise Implicit_ODE_Exception('Error, gM-J is singular at ',self._tc)
                    
            Z, W = self.calc_start_values()

            for i in xrange(self.newt):
                self._curiter += 1 #The current iteration
                self.statistics["nniter"] += 1 #Adding one iteration

                #Solve the system
                Z = N.dot(self.T2,self._radau_F(Z.real,t,y,yd))

                Z[:self._2leny]               =Z[:self._2leny]               -self._g*N.dot(self.M,W[:self._2leny])
                Z[self._2leny:2*self._2leny]  =Z[self._2leny:2*self._2leny]  -self._a*N.dot(self.M,W[self._2leny:2*self._2leny])   #+self._b*N.dot(self.I,W[2*self._leny:3*self._leny])
                Z[2*self._2leny:3*self._2leny]=Z[2*self._2leny:3*self._2leny]-self._b*N.dot(self.M,W[2*self._2leny:3*self._2leny]) #-self._a*N.dot(self.I,W[2*self._leny:3*self._leny])
                
                Z[:self._2leny]               =N.linalg.solve(self._U1,N.linalg.solve(self._L1,N.linalg.solve(self._P1,Z[:self._2leny])))
                Z[self._2leny:2*self._2leny]  =N.linalg.solve(self._U2,N.linalg.solve(self._L2,N.linalg.solve(self._P2,Z[self._2leny:2*self._2leny])))
                Z[2*self._2leny:3*self._2leny]=N.linalg.solve(self._U3,N.linalg.solve(self._L3,N.linalg.solve(self._P3,Z[2*self._2leny:3*self._2leny])))
                #----
                
                self._scaling = self._scaling/self.h**(self.index-1)#hfac
                
                newnrm = N.linalg.norm(Z.reshape(-1,self._2leny)/self._scaling,'fro')/N.sqrt(3.*self._2leny)
                
                if i > 0:
                    thq = newnrm/oldnrm
                    if i == 1:
                        self._theta = thq
                    else:
                        self._theta = N.sqrt(thq*thqold)
                    thqold = thq
                    
                    if self._theta < 0.99: #Convergence
                        self._fac_con = self._theta/(1.-self._theta)
                        dyth = self._fac_con*newnrm*self._theta**(self.newt-(i+1)-1)/self.fnewt
                        
                        if dyth >= 1.0: #Too slow convergence
                            qnewt = max(1.e-4,min(20.,dyth))
                            self._hhfac = 0.8*qnewt**(-1.0/(4.0+self.newt-(i+1)-1))
                            self.h = self._hhfac*self.h
                            self._itfail = True
                            self._rejected = True
                            break
                    else: #Not convergence, abort
                        self._itfail = True
                        break
                
                oldnrm = max(newnrm,self._eps) #Store oldnorm
                W = W+Z #Perform the iteration
                
                Z = N.dot(self.T3,W) #Calculate the new Z values
                
                if self._fac_con*newnrm <= self.fnewt: #Convergence?
                    self._itfail = False;
                    break
                
            else: #Iteration failed
                self._itfail = True
                
            if not self._itfail: #Newton iteration converged
                self._Z = Z.real
                break
            else: #Iteration failed
                self.log_message("Iteration failed at time %e with step-size %e"%(t,self.h),SCREAM)
                self.statistics["nniterfail"] += 1
                self._rejected = True #The step is rejected
                
                if self._theta >= 0.99:
                    self._hhfac = 0.5
                    self.h = self.h*self._hhfac
                if self._curjac:
                    self._needjac = False
                    self._needLU = True
                else:
                    self._needjac = True
                    self._needLU = True
        else:
            raise Implicit_ODE_Exception('Newton iteration failed at time %e with step-size %e'%(t,self.h))
    
    def estimate_error(self):
        
        temp = 1./self.h*(self.E[0]*self._Z[:self._2leny]+self.E[1]*self._Z[self._2leny:2*self._2leny]+self.E[2]*self._Z[2*self._2leny:3*self._2leny])
        temp = N.dot(self.M,temp)
        
        self._scaling = self._scaling/self.h**(self.index-1)#hfac
        
        scal = self._scaling#/self.h
        err_v = N.linalg.solve(self._U1,N.linalg.solve(self._L1,N.linalg.solve(self._P1,self._f0+temp)))
        err = N.linalg.norm(err_v/scal)
        err = max(err/N.sqrt(self._2leny),1.e-10)

        if (self._rejected or self._first) and err >= 1.: #If the step was rejected, use the more expensive error estimation
            self.statistics["nfcn"] += 1
            err_v = self._ode_f(self._tc,N.append(self._yc,self._ydc)+err_v)
            err_v = N.linalg.solve(self._U1,N.linalg.solve(self._L1,N.linalg.solve(self._P1,err_v+temp)))
            err = N.linalg.norm(err_v/scal)
            err = max(err/N.sqrt(self._2leny),1.e-10)
            
        return err
    
    def interpolate(self, t, k=0):
        """
        Calculates the continuous output from Radau5.
        """
        leny = self._2leny
        s = (t-self._newt)/self._oldh
        Z = self._col_poly
        
        diff = s*(Z[:leny]+(s-self.C[1,0]+1.)*(Z[leny:2*leny]+(s-self.C[0,0]+1.)*Z[2*leny:3*leny]))
        
        yout  = self._yc + diff[:self._leny]
        ydout = self._ydc+ diff[self._leny:]
        
        if k==0:
            return yout
        elif k==1:
            return ydout
        else:
            raise Implicit_ODE_Exception('Unknown value of k. Should be either 0 or 1')

    def jacobian(self, t, y, yd):
        """
        Calculates the Jacobian, either by an approximation or by the user
        defined (jac specified in the problem class).
        """
        self._curjac = True #The jacobian is up to date
        self._needLU = True #A new LU-decomposition is needed
        self._needjac = False #A new jacobian is not needed
        
        q = N.append(y,yd)
        
        if self.usejac: #Retrieve the user-defined jacobian
            cjac = self.problem.jac(t,y,yd)
        else:           #Calculate a numeric jacobian
            delt = N.array([(self._eps*max(abs(yi),1.e-5))**0.5 for yi in q])*N.identity(self._2leny) #Calculate a disturbance
            Fdelt = N.array([self._ode_f(t,q+e) for e in delt]) #Add the disturbance (row by row) 
            grad = ((Fdelt-self._ode_f(t,q)).T/delt.diagonal()).T
            cjac = N.array(grad).T
            self.statistics["njacfcn"] += 1+self._2leny #Add the number of function evaluations

        self.statistics["njac"] += 1 #add the number of jacobian evaluation
        return cjac
    
    def adjust_stepsize(self, err, predict=False):
        
        fac = min(self.safe, self.safe*(2.*self.newt+1.)/(2.*self.newt+self._curiter))
        quot = max(1./self.fac2,min(1./self.fac1,(err**0.25)/fac))        
        hnormal = self.h/quot
        
        if predict:
            if not self._first:
                facgus = (self._hacc/self.h)*(err**2/self._olderr)**0.25/self.safe
                facgus = max(1./self.fac2,min(1./self.fac1,facgus))
                quot = max(quot,facgus)
                h = self.h/quot
            else:
                h = hnormal
            self._hacc = self.h
        else:
            h = hnormal
        
        qt = h/self.h
        
        if (qt >= self.quot1) and (qt <= self.quot2):
            h = self.h
        
        if h > self.maxh:
            h = self.maxh
        
        if self._first and err>=1.0:
            self._hhfac = 0.1
            h = self.h*self._hhfac
        else:
            self._hhfac = h/self.h
        
        if h < self._eps:
            raise Implicit_ODE_Exception('Step-size to small at %e with h = %e'%(self._tc,self.h))
    
        return h
    
    def _collocation_pol(self, Z, col_poly, leny):

        col_poly[2*leny:3*leny] = Z[:leny] / self.C[0,0]
        col_poly[leny:2*leny]   = ( Z[:leny] - Z[leny:2*leny] ) / (self.C[0,0]-self.C[1,0])
        col_poly[:leny]         = ( Z[leny:2*leny] -Z[2*leny:3*leny] ) / (self.C[1,0]-1.)
        col_poly[2*leny:3*leny] = ( col_poly[leny:2*leny] - col_poly[2*leny:3*leny] ) / self.C[1,0]
        col_poly[leny:2*leny]   = ( col_poly[leny:2*leny] - col_poly[:leny] ) / (self.C[0,0]-1.)
        col_poly[2*leny:3*leny] =   col_poly[leny:2*leny]-col_poly[2*leny:3*leny]
        
        return col_poly
    
    def calc_start_values(self):
        """
        Calculate newton starting values.
        """
        if self._first:
            Z = N.zeros(self._2leny*3)
            W = N.zeros(self._2leny*3)
        else:
            Z = self._Z
            cq = self.C*self.h/self._oldh#self._oldoldh#self._oldh
            newtval = self._col_poly
            leny = self._2leny
            
            Z[:leny]        = cq[0,0]*(newtval[:leny]+(cq[0,0]-self.C[1,0]+1.)*(newtval[leny:2*leny]+(cq[0,0]-self.C[0,0]+1.)*newtval[2*leny:3*leny]))
            Z[leny:2*leny]  = cq[1,0]*(newtval[:leny]+(cq[1,0]-self.C[1,0]+1.)*(newtval[leny:2*leny]+(cq[1,0]-self.C[0,0]+1.)*newtval[2*leny:3*leny]))
            Z[2*leny:3*leny]= cq[2,0]*(newtval[:leny]+(cq[2,0]-self.C[1,0]+1.)*(newtval[leny:2*leny]+(cq[2,0]-self.C[0,0]+1.)*newtval[2*leny:3*leny]))
            
            W = N.dot(self.T2,Z)
            
        return Z, W
    
    def _load_parameters(self):
        
        #Parameters
        A = N.zeros([3,3])
        A[0,0] = (88.-7.*N.sqrt(6.))/360.0
        A[0,1] = (296.-169.*N.sqrt(6.))/1800.0
        A[0,2] = (-2.0+3.0*N.sqrt(6.))/225.0
        A[1,0] = (296.0+169.0*N.sqrt(6.))/1800.0
        A[1,1] = (88.+7.*N.sqrt(6.))/360.0
        A[1,2] = (-2.-3.*N.sqrt(6.))/225.0
        A[2,0] = (16.0-N.sqrt(6.))/36.0
        A[2,1] = (16.0+N.sqrt(6.))/36.0
        A[2,2] = (1.0/9.0)
        
        C = N.zeros([3,1])
        C[0,0]=(4.0-N.sqrt(6.0))/10.0
        C[1,0]=(4.0+N.sqrt(6.0))/10.0
        C[2,0]=1.0
        
        B = N.zeros([1,3])
        B[0,0]=(16.0-N.sqrt(6.0))/36.0
        B[0,1]=(16.0+N.sqrt(6.0))/36.0
        B[0,2]=1.0/9.0
        
        E = N.zeros(3)
        E[0] = -13.0-7.*N.sqrt(6.)
        E[1] = -13.0+7.0*N.sqrt(6.)
        E[2] = -1.0
        E = 1.0/3.0*E
        
        M = N.array([[1.,0.],[0.,0.]])
        
        Ainv = N.linalg.inv(A)
        [eig, T] = N.linalg.eig(Ainv)
        eig = N.array([eig[2],eig[0],eig[1]])
        J = N.diag(eig)

        self._alpha = eig[1]
        self._beta  = eig[2]
        self._gamma = eig[0].real
        
        temp0 = T[:,0].copy()
        temp1 = T[:,1].copy()
        temp2 = T[:,2].copy()
        T[:,0] = temp2
        T[:,1] = temp0
        T[:,2] = temp1
        Tinv = N.linalg.inv(T)
        
        I = N.eye(self._2leny)
        M = N.kron(M,N.eye(self._leny))
        I3 = N.eye(3)
        T1 = N.kron(J,M)
        T2 = N.kron(Tinv,I)
        T3 = N.kron(T,I)
        
        self.A = A
        self.B = B
        self.C = C
        self.I = I
        self.E = E
        self.M = M
        self.T1 = T1
        self.T2 = T2
        self.T3 = T3
        self.I3 = I3
        self.EIG = eig

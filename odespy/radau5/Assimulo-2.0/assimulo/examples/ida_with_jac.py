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
from assimulo.solvers.sundials import IDA
from assimulo.problem import Implicit_Problem
import nose

def run_example(with_plots=True):
    
    #Defines the residual
    def f(t,y,yd):
        
        res_0 = yd[0]-y[2]
        res_1 = yd[1]-y[3]
        res_2 = yd[2]+y[4]*y[0]
        res_3 = yd[3]+y[4]*y[1]+9.82
        #res_4 = y[0]**2+y[1]**2-1
        res_4 = y[2]**2+y[3]**2-y[4]*(y[0]**2+y[1]**2)-y[1]*9.82

        return N.array([res_0,res_1,res_2,res_3,res_4])
    
    #Defines the jacobian
    def jac(c,t,y,yd):
        jacobian = N.zeros([len(y),len(y)])
        
        #Derivative
        jacobian[0,0] = 1*c
        jacobian[1,1] = 1*c
        jacobian[2,2] = 1*c
        jacobian[3,3] = 1*c
        
        #Differentiated
        jacobian[0,2] = -1
        jacobian[1,3] = -1
        jacobian[2,0] = y[4]
        jacobian[3,1] = y[4]
        jacobian[4,0] = y[0]*2*y[4]*-1
        jacobian[4,1] = y[1]*2*y[4]*-1-9.82
        jacobian[4,2] = y[2]*2
        jacobian[4,3] = y[3]*2
        
        #Algebraic
        jacobian[2,4] = y[0]
        jacobian[3,4] = y[1]
        jacobian[4,4] = -(y[0]**2+y[1]**2)
        
        return jacobian
        
    #The initial conditons
    y0 = [1.0,0.0,0.0,0.0,5] #Initial conditions
    yd0 = [0.0,0.0,0.0,-9.82,0.0] #Initial conditions
    
    #Create an Assimulo implicit problem
    imp_mod = Implicit_Problem(f,y0,yd0)
    
    #Sets the options to the problem
    imp_mod.jac = jac #Sets the jacobian
    imp_mod.algvar = [1.0,1.0,1.0,1.0,0.0] #Set the algebraic components
    imp_mod.name = 'Test Jacobian'
    
    #Create an Assimulo implicit solver (IDA)
    imp_sim = IDA(imp_mod) #Create a IDA solver
    
    #Sets the paramters
    imp_sim.atol = 1e-6 #Default 1e-6
    imp_sim.rtol = 1e-6 #Default 1e-6
    imp_sim.suppress_alg = True #Suppres the algebraic variables on the error test
    
    #Let Sundials find consistent initial conditions by use of 'IDA_YA_YDP_INIT'
    imp_sim.make_consistent('IDA_YA_YDP_INIT')
    
    #Simulate
    t, y, yd = imp_sim.simulate(5,1000) #Simulate 5 seconds with 1000 communication points
    
    #Basic tests
    nose.tools.assert_almost_equal(y[-1][0],0.9401995, places=4)
    nose.tools.assert_almost_equal(y[-1][1],-0.34095124, places=4)
    nose.tools.assert_almost_equal(yd[-1][0], -0.88198927, places=4)
    nose.tools.assert_almost_equal(yd[-1][1], -2.43227069, places=4)
    
    #Plot
    if with_plots:
        P.plot(t,y)
        P.show()


if __name__=='__main__':
    run_example()
    

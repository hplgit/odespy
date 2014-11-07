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
from assimulo.solvers.sundials import CVode
from assimulo.problem import Explicit_Problem


def run_example(with_plots=True):
    
    #Defines the rhs
    def f(t,y):
        yd_0 = y[1]
        yd_1 = -9.82

        return N.array([yd_0,yd_1])
    
    #Defines the jacobian*vector product
    def jacv(t,y,fy,v):
        j = N.array([[0,1.],[0,0]])
        return N.dot(j,v)
    
    y0 = [1.0,0.0] #Initial conditions
    
    #Defines an Assimulo explicit problem
    exp_mod = Explicit_Problem(f,y0)
    
    exp_mod.jacv = jacv #Sets the jacobian
    exp_mod.name = 'Example using the Jacobian Vector product'
    
    exp_sim = CVode(exp_mod) #Create a CVode solver
    
    #Set the parameters
    exp_sim.iter = 'Newton' #Default 'FixedPoint'
    exp_sim.discr = 'BDF' #Default 'Adams'
    exp_sim.atol = 1e-5 #Default 1e-6
    exp_sim.rtol = 1e-5 #Default 1e-6
    exp_sim.linear_solver = 'SPGMR' #Change linear solver
    #exp_sim.options["usejac"] = False
    
    #Simulate
    t, y = exp_sim.simulate(5, 1000) #Simulate 5 seconds with 1000 communication points
    
    #Basic tests
    nose.tools.assert_almost_equal(y[-1][0],-121.75000000,4)
    nose.tools.assert_almost_equal(y[-1][1],-49.100000000)
    
    #Plot
    if with_plots:
        P.plot(t,y)
        P.show()

if __name__=='__main__':
    run_example()

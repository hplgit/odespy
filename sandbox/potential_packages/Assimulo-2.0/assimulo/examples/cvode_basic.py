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
#import pylab as P
import nose
from assimulo.solvers.sundials import CVode
from assimulo.problem import Explicit_Problem

def run_example(with_plots=True):
    global exp_mod, exp_sim
    
    #Define the rhs
    def f(t,y):
        ydot = -y[0]
        return N.array([ydot])
    
    #Define an Assimulo problem
    exp_mod = Explicit_Problem(f, y0=4)
    exp_mod.name = 'Simple CVode Example'
    
    #Define an explicit solver
    exp_sim = CVode(exp_mod) #Create a CVode solver
    
    #Sets the parameters
    exp_sim.iter  = 'Newton' #Default 'FixedPoint'
    exp_sim.discr = 'BDF' #Default 'Adams'
    exp_sim.atol = [1e-4] #Default 1e-6
    exp_sim.rtol = 1e-4 #Default 1e-6

    #Simulate
    t1, y1 = exp_sim.simulate(5,100) #Simulate 5 seconds
    t2, y2 = exp_sim.simulate(7) #Simulate 2 seconds more
    
    #Basic test
    nose.tools.assert_almost_equal(y2[-1], 0.00347746, 5)
    
    #Plot
    #if with_plots:
    #    P.plot(t1, y1, color="b")
    #    P.plot(t2, y2, color="b")
    #    P.show()

if __name__=='__main__':
    run_example()

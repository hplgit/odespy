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
from assimulo.solvers.runge_kutta import *
from assimulo.problem import Explicit_Problem

def run_example(with_plots=True):
        
    #Defines the rhs
    def f(t,y):
        ydot = -y[0]
        return N.array([ydot])

    #Define an Assimulo problem
    exp_mod = Explicit_Problem(f, 4.0)
    exp_mod.name = 'Simple Explicit Example'
    
    exp_sim = RungeKutta34(exp_mod) #Create a RungeKutta34 solver
    exp_sim.inith = 0.1 #Sets the initial step, default = 0.01
    
    #Simulate
    t, y = exp_sim.simulate(5) #Simulate 5 seconds
    
    #Basic test
    nose.tools.assert_almost_equal(y[-1],0.02695199,5)
    
    #Plot
    if with_plots:
        P.plot(t,y)
        P.show()

if __name__=='__main__':
    run_example()

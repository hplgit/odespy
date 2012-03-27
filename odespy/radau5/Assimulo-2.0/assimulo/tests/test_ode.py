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

import nose
from assimulo import testattr
from assimulo.ode import *
from assimulo.problem import Explicit_Problem
from assimulo.exception import *

class Test_ODE:
    
    def setUp(self):
        self.problem = Explicit_Problem(y0=4.0)
        self.simulator = ODE(self.problem)
    
    @testattr(stddist = True)
    def test_init(self):
        """
        This tests the functionality of the method __init__.
        """
        assert self.simulator.verbosity == NORMAL
        assert self.simulator.continuous_output == False
    
    @testattr(stddist = True)
    def test_verbosity(self):
        """
        This tests the functionality of the property verbosity.
        """
        nose.tools.assert_raises(AssimuloException, self.simulator._set_verbosity, 'Test')
        nose.tools.assert_raises(AssimuloException, self.simulator._set_verbosity, [1, 31])
        nose.tools.assert_raises(AssimuloException, self.simulator._set_verbosity, [1])
        
        self.simulator.verbosity=1
        assert self.simulator.verbosity==1
        assert self.simulator.options["verbosity"] == 1
        self.simulator.verbosity=4
        assert self.simulator.verbosity==4
        assert self.simulator.options["verbosity"] == 4
        
    @testattr(stddist = True)    
    def test_continuous_output(self):
        """
        This tests the functionality of the property continuous_output.
        """
        assert self.simulator.continuous_output == False #Test the default value
        
        self.simulator.continuous_output = True
        assert self.simulator.continuous_output == True
        assert self.simulator.options["continuous_output"] == True

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
from assimulo.exception import *
from assimulo.examples import *

class Test_Examples:
    
    @testattr(stddist = True)
    def test_cvode_gyro(self):
        cvode_gyro.run_example(with_plots=False)
    
    @testattr(stddist = True)
    def test_cvode_basic(self):
        cvode_basic.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_cvode_with_disc(self):
        cvode_with_disc.run_example(with_plots=False)
    
    @testattr(stddist = True)
    def test_cvode_with_initial_sensitivity(self):
        cvode_with_initial_sensitivity.run_example(with_plots=False)
    
    @testattr(stddist = True)
    def test_cvode_with_jac(self):
        cvode_with_jac.run_example(with_plots=False)
    
    @testattr(stddist = True)
    def test_cvode_with_jac_spgmr(self):
        cvode_with_jac_spgmr.run_example(with_plots=False)
    
    @testattr(stddist = True)
    def test_cvode_with_parameters(self):
        cvode_with_parameters.run_example(with_plots=False)

    @testattr(stddist = True)
    def test_euler_basic(self):
        euler_basic.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_rungekutta4_basic(self):
        rungekutta4_basic.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_rungekutta34_basic(self):
        rungekutta34_basic.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_ida_with_disc(self):
        ida_with_disc.run_example(with_plots=False)
    
    @testattr(stddist = True)
    def test_ida_with_initial_sensitivity(self):
        ida_with_initial_sensitivity.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_ida_with_jac(self):
        ida_with_jac.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_ida_with_parameters(self):
        ida_with_parameters.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_radau5ode_vanderpol(self):
        radau5ode_vanderpol.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_radau5dae_vanderpol(self):
        radau5dae_vanderpol.run_example(with_plots=False)

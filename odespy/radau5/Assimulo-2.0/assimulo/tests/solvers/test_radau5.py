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
from assimulo.solvers.radau5 import *
from assimulo.problem import Explicit_Problem
from assimulo.problem import Implicit_Problem
from assimulo.exception import *
from assimulo.lib.radau_core import Radau_Exception, Radau_Common


class Test_Explicit_Radau5:
    """
    Tests the explicit Radau solver.
    """
    def setUp(self):
        """
        This sets up the test case.
        """
        def f(t,y):
            eps = 1.e-6
            my = 1./eps
            yd_0 = y[1]
            yd_1 = my*((1.-y[0]**2)*y[1]-y[0])
            
            return N.array([yd_0,yd_1])
        
        def jac(t,y):
            eps = 1.e-6
            my = 1./eps
            J = N.zeros([2,2])
            
            J[0,0]=0.
            J[0,1]=1.
            J[1,0]=my*(-2.*y[0]*y[1]-1.)
            J[1,1]=my*(1.-y[0]**2)
            
            return J
        
        #Define an Assimulo problem
        y0 = [2.0,-0.6] #Initial conditions
        
        exp_mod = Explicit_Problem(f,y0)
        exp_mod_t0 = Explicit_Problem(f,y0,1.0)
        
        exp_mod.jac = jac
        self.mod = exp_mod
            
        #Define an explicit solver
        self.sim = Radau5ODE(exp_mod) #Create a Radau5 solve
        self.sim_t0 = Radau5ODE(exp_mod_t0)
        
        #Sets the parameters
        self.sim.atol = 1e-4 #Default 1e-6
        self.sim.rtol = 1e-4 #Default 1e-6
        self.sim.inith = 1.e-4 #Initial step-size
        self.sim.usejac = False
    
    @testattr(stddist = True)
    def test_time_event(self):
        f = lambda t,y: [1.0]
        global tnext
        global nevent
        tnext = 0.0
        nevent = 0
        def time_events(t,y,sw):
            global tnext,nevent
            events = [1.0, 2.0, 2.5, 3.0]
            for ev in events:
                if t < ev:
                    tnext = ev
                    break
                else:
                    tnext = None
            nevent += 1
            return tnext
            
        def handle_event(solver, event_info):
            solver.y+= 1.0
            global tnext
            nose.tools.assert_almost_equal(solver.t, tnext)
            assert event_info[0] == []
            assert event_info[1] == True
    
        exp_mod = Explicit_Problem(f,0.0)
        exp_mod.time_events = time_events
        exp_mod.handle_event = handle_event
        
        #CVode
        exp_sim = Radau5ODE(exp_mod)
        exp_sim(5.,100)
        
        assert nevent == 5
    
    @testattr(stddist = True)
    def test_init(self):
        
        #Test both y0 in problem and not.
        sim = Radau5ODE(self.mod)
        
        assert sim._leny == 2
    
    @testattr(stddist = True)
    def test_collocation_polynomial(self):
        """
        This tests the functionality of the collocation polynomial (communication points)
        """
        self.sim.continuous_output = False
        
        self.sim.simulate(2.,200) #Simulate 2 seconds
        
        assert self.sim.statistics["nsteps"] < 300
        
        #nose.tools.assert_almost_equal(self.sim.y[-2][0], 1.71505001, 4)
        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], 1.7061680350, 4)
        
        self.sim.continuous_output = True
        self.sim.reset()
        self.sim.simulate(2.,200) #Simulate 2 seconds
        
        assert self.sim.statistics["nsteps"] < 300

        #nose.tools.assert_almost_equal(self.sim.y[-2][0], 1.71505001, 4)
        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], 1.7061680350, 4)
        
        self.sim_t0.simulate(3.)
        nose.tools.assert_almost_equal(self.sim_t0.t_sol[0], 1.0000000, 4)
        nose.tools.assert_almost_equal(self.sim_t0.t_sol[-1], 3.0000000, 4)
        nose.tools.assert_almost_equal(self.sim_t0.y_sol[-1][0], 1.7061680350, 4)
        
    @testattr(stddist = True)
    def test_simulation(self):
        """
        This tests the Radau5 with a simulation of the van der pol problem.
        """
        self.sim.simulate(2.) #Simulate 2 seconds
        
        assert self.sim.statistics["nsteps"] < 300

        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], 1.7061680350, 4)
    
    @testattr(stddist = True)    
    def test_simulation_ncp(self):
        """
        Test a simulation with ncp.
        """
        self.sim.continuous_output = True
        
        self.sim.simulate(1.0, 200) #Simulate 1 second
        assert len(self.sim.t_sol) == 201
        
        self.sim.reset()
        self.sim.continuous_output = False
        
        self.sim.simulate(1.0, 200) #Simulate 1 second
        assert len(self.sim.t_sol) == 201
    
    @testattr(stddist = True)
    def test_usejac(self):
        """
        This tests the usejac property.
        """
        self.sim.usejac = True
        
        self.sim.simulate(2.) #Simulate 2 seconds

        assert self.sim.statistics["njacfcn"] == 0
        
        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], 1.7061680350, 4)

    @testattr(stddist = True)
    def test_thet(self):
        """
        This tests a negative value of thet.
        """
        self.sim.thet = -1
        self.sim.simulate(2.) #Simulate 2 seconds

        assert self.sim.statistics["nsteps"] == self.sim.statistics["njac"]
    
    @testattr(stddist = True)
    def test_maxh(self):
        """
        This tests the maximum step length.
        """
        self.sim.maxh = 0.01
        self.sim.simulate(0.5)
        assert max(N.diff(self.sim.t_sol))-N.finfo('double').eps <= 0.01
        
    @testattr(stddist = True)
    def test_newt(self):
        """
        This tests the maximum number of newton iterations.
        """
        self.sim.newt = 10
        self.sim.simulate(1.0)
        
        assert self.sim.statistics["nniterfail"] == 1
    
    @testattr(stddist = True)
    def test_safe(self):
        """
        This tests the safety factor in the step-size prediction.
        """
        self.sim.safe = 0.99
        self.sim.simulate(1.0)
        assert self.sim.statistics["nsteps"] < 150
        
    @testattr(stddist = True)
    def test_reset_statistics(self):
        """
        Tests that the statistics are reset.
        """
        self.sim.simulate(1.0)
        steps = self.sim.statistics["nsteps"]
        
        self.sim.reset()
        self.sim.simulate(1.0)
        
        assert self.sim.statistics["nsteps"] < steps*1.5
        
    @testattr(stddist = True)
    def test_atol(self):
        """
        This test the absolute tolerance.
        """
        self.sim.simulate(1.0)
        
        steps = self.sim.statistics["nsteps"]
        
        self.sim.reset()
        
        self.sim.rtol = 1e-8
        self.sim.atol = 1e-8
        
        self.sim.simulate(1.0)
        steps2 = self.sim.statistics["nsteps"]
        
        assert steps2 > steps
        
        self.sim.reset()
        self.sim.atol = [1e-8, 1e-8]
        
        steps3 = self.sim.statistics["nsteps"]
        
        assert steps3==steps2
        
        nose.tools.assert_raises(Radau_Exception, self.sim._set_atol, [1e-6,1e-6,1e-6])
    
class Test_Implicit_Radau5:
    """
    Tests the implicit Radau solver.
    """
    def setUp(self):
        """
        This sets up the test case.
        """
        #Define the residual
        def f(t,y,yd):
            eps = 1.e-6
            my = 1./eps
            yd_0 = y[1]
            yd_1 = my*((1.-y[0]**2)*y[1]-y[0])
            
            res_0 = yd[0]-yd_0
            res_1 = yd[1]-yd_1
            
            return N.array([res_0,res_1])
        
        y0 = [2.0,-0.6] #Initial conditions
        yd0 = [-.6,-200000.]
        
        #Define an Assimulo problem
        self.mod = Implicit_Problem(f,y0,yd0)
        self.mod_t0 = Implicit_Problem(f,y0,yd0,1.0)
            
        #Define an explicit solver
        self.sim = Radau5DAE(self.mod) #Create a Radau5 solve
        self.sim_t0 = Radau5DAE(self.mod_t0)
        
        #Sets the parameters
        self.sim.atol = 1e-4 #Default 1e-6
        self.sim.rtol = 1e-4 #Default 1e-6
        self.sim.inith = 1.e-4 #Initial step-size
    
    @testattr(stddist = True)
    def test_time_event(self):
        f = lambda t,y,yd: y-yd
        global tnext
        global nevent
        tnext = 0.0
        nevent = 0
        def time_events(t,y,yd,sw):
            global tnext,nevent
            events = [1.0, 2.0, 2.5, 3.0]
            for ev in events:
                if t < ev:
                    tnext = ev
                    break
                else:
                    tnext = None
            nevent += 1
            print tnext
            return tnext
            
        def handle_event(solver, event_info):
            #solver.y+= 1.0
            global tnext
            nose.tools.assert_almost_equal(solver.t, tnext)
            assert event_info[0] == []
            assert event_info[1] == True
    
        exp_mod = Implicit_Problem(f,0.0,0.0)
        exp_mod.time_events = time_events
        exp_mod.handle_event = handle_event
        
        #CVode
        exp_sim = Radau5DAE(exp_mod)
        exp_sim.verbosity = 0
        exp_sim(5.,100)
        
        assert nevent == 5
    
    @testattr(stddist = True)
    def test_init(self):
        """
        This tests the functionality of Radau5 Implicit Init.
        """
        #Test both y0 in problem and not.

        sim = Radau5DAE(self.mod)
        
        assert sim._leny == 2
    
    @testattr(stddist = True)
    def test_thet(self):
        """
        This tests a negative value of thet.
        """
        self.sim.thet = -1
        self.sim.simulate(.5) #Simulate 2 seconds

        assert self.sim.statistics["nsteps"] == self.sim.statistics["njac"]
        
    @testattr(stddist = True)    
    def test_simulation(self):
        """
        Test a simulation of the vanderpol equations.
        """
        #Simulate
        self.sim.simulate(2.) #Simulate 2 seconds
        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], 1.706272, 3)
        
        self.sim.reset()
        
        self.sim.continuous_output = True
        
        #Simulate
        self.sim.simulate(2.) #Simulate 2 seconds
        nose.tools.assert_almost_equal(self.sim.y_sol[-1][0], 1.706272, 3)
        
        self.sim_t0.simulate(3.)
        nose.tools.assert_almost_equal(self.sim_t0.t_sol[0], 1.0000000, 4)
        nose.tools.assert_almost_equal(self.sim_t0.t_sol[-1], 3.0000000, 4)
        nose.tools.assert_almost_equal(self.sim_t0.y_sol[-1][0], 1.7061680350, 4)
    
    @testattr(stddist = True)    
    def test_simulation_ncp(self):
        """
        Test a simulation with ncp.
        """
        self.sim.continuous_output = True
        
        self.sim.simulate(1.0, 200) #Simulate 1 second
        assert len(self.sim.t_sol) == 201
        
        self.sim.reset()
        self.sim.continuous_output = False
        
        self.sim.simulate(1.0, 200) #Simulate 1 second
        assert len(self.sim.t_sol) == 201
    
    @testattr(stddist = True)
    def test_maxh(self):
        """
        Tests implicit radau with maxh.
        """
        self.sim.maxh = 0.01
        self.sim.simulate(0.5)
        assert max(N.diff(self.sim.t_sol))-N.finfo('double').eps <= 0.01

class Test_Radau_Common:
    """
    Tests the common attributes of the Radau solvers.
    """
    def setUp(self):
        """
        This sets up the test case.
        """

        f = lambda t,y:[1.0,2.0]

        #Define an Assimulo problem
        y0 = [2.0,-0.6] #Initial conditions
        exp_mod = Explicit_Problem(f,y0)
        self.sim = Radau5ODE(exp_mod)
    
    @testattr(stddist = True)    
    def test_fac1(self):
        """
        This tests the functionality of the property fac1.
        """
        self.sim.fac1 = 0.01
        assert self.sim.fac1 == 0.01
        self.sim.fac1 = 0.001
        assert self.sim.fac1 == 0.001
        
        nose.tools.assert_raises(Radau_Exception, self.sim._set_fac1, 'Test')
        nose.tools.assert_raises(Radau_Exception, self.sim._set_fac1, [-1.0])
    
    @testattr(stddist = True)
    def test_fac2(self):
        """
        This tests the functionality of the property fac2.
        """
        self.sim.fac2 = 0.01
        assert self.sim.fac2 == 0.01
        self.sim.fac2 = 0.001
        assert self.sim.fac2 == 0.001
        
        nose.tools.assert_raises(Radau_Exception, self.sim._set_fac2, 'Test')
        nose.tools.assert_raises(Radau_Exception, self.sim._set_fac2, [-1.0])
    
    @testattr(stddist = True)
    def test_fnewt(self):
        """
        This tests the functionality of the property fnewt.
        """
        self.sim.fnewt = 0.01
        assert self.sim.fnewt == 0.01
        self.sim.fnewt = 0.001
        assert self.sim.fnewt == 0.001
        
        nose.tools.assert_raises(Radau_Exception, self.sim._set_fnewt, 'Test')
        nose.tools.assert_raises(Radau_Exception, self.sim._set_fnewt, [-1.0])
    
    @testattr(stddist = True)
    def test_h(self):
        """
        This tests the functionality of the property h.
        """
        self.sim.h = 0.01
        assert self.sim.h == 0.01
        self.sim.h = 0.001
        assert self.sim.h == 0.001
    
    @testattr(stddist = True)
    def test_initial_step(self):
        """
        This tests the functionality of the property initial step.
        """
        self.sim.inith = 0.01
        assert self.sim.inith == 0.01
        self.sim.inith = 0.001
        assert self.sim.inith == 0.001
        
        nose.tools.assert_raises(Radau_Exception, self.sim._set_initial_step, 'Test')
        nose.tools.assert_raises(Radau_Exception, self.sim._set_initial_step, [-1.0])
    
    @testattr(stddist = True)
    def test_newt(self):
        """
        This tests the functionality of the property newt.
        """
        self.sim.newt = 1
        assert self.sim.newt == 1
        self.sim.newt = 10
        assert self.sim.newt == 10
        self.sim.newt = 9.8
        assert self.sim.newt == 9
        
        nose.tools.assert_raises(Radau_Exception, self.sim._set_newt, 'Test')
        nose.tools.assert_raises(Radau_Exception, self.sim._set_newt, [-1.0])
    
    @testattr(stddist = True)
    def test_quot1(self):
        """
        This tests the functionality of the property quot1.
        """
        self.sim.quot1 = 0.01
        assert self.sim.quot1 == 0.01
        self.sim.quot1 = 0.001
        assert self.sim.quot1 == 0.001
        
        nose.tools.assert_raises(Radau_Exception, self.sim._set_quot1, 'Test')
        nose.tools.assert_raises(Radau_Exception, self.sim._set_quot1, [-1.0])
    
    @testattr(stddist = True)    
    def test_quot2(self):
        """
        This tests the functionality of the property quot2.
        """
        self.sim.quot2 = 0.01
        assert self.sim.quot2 == 0.01
        self.sim.quot2 = 0.001
        assert self.sim.quot2 == 0.001
        
        nose.tools.assert_raises(Radau_Exception, self.sim._set_quot2, 'Test')
        nose.tools.assert_raises(Radau_Exception, self.sim._set_quot2, [-1.0])
    
    @testattr(stddist = True)
    def test_safe(self):
        """
        This tests the functionality of the property safe.
        """
        self.sim.safe = 0.01
        assert self.sim.safe == 0.01
        self.sim.safe = 0.001
        assert self.sim.safe == 0.001
        
        nose.tools.assert_raises(Radau_Exception, self.sim._set_safe, 'Test')
        nose.tools.assert_raises(Radau_Exception, self.sim._set_safe, [-1.0])
    
    @testattr(stddist = True)
    def test_thet(self):
        """
        This tests the functionality of the property thet.
        """
        self.sim.thet = 0.01
        assert self.sim.thet == 0.01
        self.sim.thet = 0.001
        assert self.sim.thet == 0.001
        
        nose.tools.assert_raises(Radau_Exception, self.sim._set_thet, 'Test')
        nose.tools.assert_raises(Radau_Exception, self.sim._set_thet, [-1.0])
    
    @testattr(stddist = True)
    def test_usejac(self):
        """
        This tests the functionality of the property usejac.
        """
        self.sim.usejac = True
        assert self.sim.usejac == True
        self.sim.usejac = False
        assert self.sim.usejac == False
        self.sim.usejac = 1
        assert self.sim.usejac == True
        self.sim.usejac = []
        assert self.sim.usejac == False


import nose
import numpy as N
from assimulo.problem_algebraic import *
from assimulo.kinsol import *


class Test_KINSOL:
    
    def setUp(self):
        
        """
        sets up the test case
        """
          
        class Prob1(ProblemAlgebraic):
            f = lambda self, x:N.array([x[0]-1.0, x[1]-2.0,x[2]-3.0])
            _x0 = [0.0,0.0,0.0]
            
            def set_x0(self,x0):
                self._x0 = x0
                
            def get_x0(self):
                return self._x0
            
            def print_var_info(self,i):
                print "No specific info to print"
        self.p1 = Prob1()
        self.solve_p1 = KINSOL(self.p1)
        
        class Prob_no_f(ProblemAlgebraic):
            _x0 = [0.0,0.0,0.0]
            
            def set_x0(self,_x0):
                self._x0 = x0
                
            def get_x0(self):
                return self._x0
            
            def print_var_info(self,i):
                print "No specific info to print"
            
        self.no_f = Prob_no_f()
        #self.solve_no_f = KINSOL(self.no_f)
        
        class Prob_no_x0(ProblemAlgebraic):
            f = lambda self, x:N.array([x[0]-1.0, x[1]-2.0,x[2]-3.0])
            
            def set_x0(self,x0):
                self._x0 =x0
                
            def get_x0(self):
                return self._x0
            def print_var_info(self,i):
                print "No specific info to print"
            
        self.no_x0 = Prob_no_f()
        #self.solve_no_x0 = KINSOL(self.no_x0)
        
        class Prob_no_getx0(ProblemAlgebraic):
            f = lambda self, x:N.array([x[0]-1.0, x[1]-2.0,x[2]-3.0])
            _x0 = [0.0,0.0,0.0]
            
            def set_x0(self,x0):
                self._x0 = x0
                
            def print_var_info(self,i):
                print "No specific info to print"
                
        self.no_getx0 = Prob_no_getx0()
        #self.solve_no_getx0 = KINSOL(self.no_getx0)
        
        class Prob_Jac(ProblemAlgebraic):
            f = lambda self, x:N.array([x[0]-1.0, x[1]-2.0,x[2]-3.0])
            _x0 = [0.0,0.0,0.0]
            jac_called = False
            
            def jac(self,x):
                self.jac_called = True
                print "Jacobian called"
                return N.eye(3)
            
            def set_x0(self,x0):
                self._x0 = x0
                
            def get_x0(self):
                return self._x0
            
            def print_var_info(self,i):
                print "No specific info to print"
            
        self.pb_jac = Prob_Jac()
        self.solve_pb_jac = KINSOL(self.pb_jac)
        
        class Prob_Const(ProblemAlgebraic):
            f = lambda self, x:N.array([-(x[0]-1.0)**2 +4, x[1]-2.0,x[2]-3.0])
            _x0 = [0.0,1.0,1.0]
            _const = None
            
            def set_x0(self,x0):
                self._x0 = x0
                
            def get_x0(self):
                return self._x0
            
            def set_constraints(self,const):
                self._const = const
                
            def get_constraints(self):
                return self._const
            
            def print_var_info(self,i):
                print self._x0[i], "constrained by ", self._const[i]
                
        self.pb_const = Prob_Const()
        self.solve_pb_const = KINSOL(self.pb_const)
        
        class Prob_BadConst(ProblemAlgebraic):
            f = lambda self, x:N.array([-(x[0]-1.0)**2 +4, x[1]-2.0,x[2]-3.0, x[3]+1,x[4]+2])
            _x0 = [-1.0,-1.0,0.0,1.0,1.0]
            _const = N.array([-2.0,-1.0,0.0,1.0,2.0])
            
            def set_x0(self,x0):
                self._x0 = x0
                
            def get_x0(self):
                return self._x0
            
            def set_constraints(self,const):
                self._const = const
                
            def get_constraints(self):
                return self._const
            
            def print_var_info(self,i):
                print self._x0[i], "constrained by ", self._const[i]
                
        self.pb_badconst = Prob_BadConst()
        
        class Prob_Sparse(ProblemAlgebraic):
            f = lambda self, x:N.array([-(x[0]-1.0)**2 +4, x[1]-2.0,x[2]-3.0, x[3]+1,x[4]+2])
            _x0 = [-1.0,-1.0,0.0,1.0,1.0]
            _const = N.array([-2.0,-1.0,0.0,1.0,2.0])
            
            def set_x0(self,x0):
                self._x0 = x0
                
            def get_x0(self):
                return self._x0
            
            def set_constraints(self,const):
                self._const = const
                
            def get_constraints(self):
                return self._const
            
            def sparse_jac(self):
                pass
            
            def print_var_info(self,i):
                print self._x0[i], "constrained by ", self._const[i]
                
        self.pb_sparse = Prob_Sparse()
        self.solve_pb_sparse = KINSOL(self.pb_sparse)
        
          
    def test_solve(self):
        """
        Test if solve works in kinsol.py
        """
        
        # tests for problems without all member attributes/methods
        nose.tools.assert_raises(KINSOL_Exception,KINSOL, self.no_f)
        nose.tools.assert_raises(KINSOL_Exception,KINSOL, self.no_x0)
        nose.tools.assert_raises(KINSOL_Exception,KINSOL, self.no_getx0)
        
        # test solver for simple problem
        res = self.solve_p1.solve()
        nose.tools.assert_almost_equal(res[0],1.0,5)
        nose.tools.assert_almost_equal(res[1],2.0,5)
        nose.tools.assert_almost_equal(res[2],3.0,5)
        
        self.p1.set_x0('a string')
        nose.tools.assert_raises(KINSOL_Exception,KINSOL, self.p1)
        
        self.p1.set_x0(5)
        nose.tools.assert_raises(KINSOL_Exception,KINSOL, self.p1)
        
        self.p1.set_x0(0.6)
        nose.tools.assert_raises(KINSOL_Exception,KINSOL, self.p1)
        
    def test_jac_usage(self):
        """
        Tests if user-supplied jacobians are implemented correctly in kinsol.py
        """
        # test is supplied jacobian is called
        try:
            self.solve_pb_jac.set_jac_usage(True)
            self.solve_pb_jac.solve()
        except :
            pass
        
        nose.tools.assert_true(self.pb_jac.jac_called,'Jacobian not used although use_jac = true')
        
        self.pb_jac.jac_called = False
        try:
            self.solve_pb_jac.set_jac_usage(False)
            self.solve_pb_jac.solve()
        except :
            pass
        
        nose.tools.assert_false(self.pb_jac.jac_called,'Jacobian used although use_jac = false')
        
    def test_constraints_usage(self):
        """
        Tests if constraints are implemented correctly in kinsol.py
        """
        res =self.solve_pb_const.solve()
        nose.tools.assert_almost_equal(res[0],-1.0,5)
        nose.tools.assert_almost_equal(res[1],2.0,5)
        nose.tools.assert_almost_equal(res[2],3.0,5)
        
        self.pb_const.set_constraints(5)
        nose.tools.assert_raises(KINSOL_Exception,KINSOL,self.pb_const)
        
        self.pb_const.set_constraints('a')
        nose.tools.assert_raises(KINSOL_Exception,KINSOL,self.pb_const)
        
        self.pb_const.set_constraints(True)
        nose.tools.assert_raises(KINSOL_Exception,KINSOL,self.pb_const)
        
        self.pb_const.set_constraints([1.,1.,1.])
        nose.tools.assert_raises(KINSOL_Exception,KINSOL,self.pb_const)
        
        self.pb_const.set_constraints(N.ones(2))
        nose.tools.assert_raises(KINSOL_Exception,KINSOL,self.pb_const)
        
        self.pb_const.set_constraints(N.ones(3,dtype = int))
        nose.tools.assert_raises(KINSOL_Exception,KINSOL,self.pb_const)
        
        self.pb_const.set_constraints(-N.ones(3,dtype = float))
        nose.tools.assert_raises(KINSOL_Exception,KINSOL,self.pb_const)
        
        self.pb_const.set_constraints(N.ones(3))
        print "const: ",self.pb_const._const
        self.pb_const.set_x0(N.array([1.5,1.0,1.0]))
        self.solve_pb_const = KINSOL(self.pb_const)
        
        res2 = self.solve_pb_const.solve()
        nose.tools.assert_almost_equal(res2[0],3.0,5)
        nose.tools.assert_almost_equal(res2[1],2.0,5)
        nose.tools.assert_almost_equal(res2[2],3.0,5)
        
    def test_verbosity_usage(self):
        """
        test if the setting of printout level works
        """
        # test default value
        nose.tools.assert_equal(self.solve_p1.verbosity,0)
        
        # Test for faulty input
        nose.tools.assert_raises(KINSOL_Exception,self.solve_p1.set_verbosity,1.0)
        nose.tools.assert_raises(KINSOL_Exception,self.solve_p1.set_verbosity,True)
        nose.tools.assert_raises(KINSOL_Exception,self.solve_p1.set_verbosity,'a')
        nose.tools.assert_raises(KINSOL_Exception,self.solve_p1.set_verbosity,4)
        nose.tools.assert_raises(KINSOL_Exception,self.solve_p1.set_verbosity,-1)
        nose.tools.assert_raises(KINSOL_Exception,self.solve_p1.set_verbosity,N.ones(3))
        nose.tools.assert_raises(KINSOL_Exception,self.solve_p1.set_verbosity,N.ones(3,dtype = int))
        
        # Test if set correctly
        self.solve_p1.set_verbosity(1)
        nose.tools.assert_equal(self.solve_p1.verbosity,1)
        
        self.solve_p1.set_verbosity(2)
        nose.tools.assert_equal(self.solve_p1.verbosity,2)
        
        self.solve_p1.set_verbosity(3)
        nose.tools.assert_equal(self.solve_p1.verbosity,3)
        
        self.solve_p1.set_verbosity(0)
        nose.tools.assert_equal(self.solve_p1.verbosity,0)
        
    def test_constraints_check(self):
        """
        test if the cosntraint checking works properly
        """
        
        bad_solver = KINSOL(self.pb_badconst)
        
        self.pb_badconst.set_x0([1.0,-1.0,0.0,1.0,1.0])
        nose.tools.assert_raises(KINSOL_Exception,KINSOL,self.pb_badconst)
        self.pb_badconst.set_x0([-1.0,1.0,0.0,1.0,1.0])
        nose.tools.assert_raises(KINSOL_Exception,KINSOL,self.pb_badconst)
        self.pb_badconst.set_x0([-1.0,-1.0,0.0,-1.0,1.0])
        nose.tools.assert_raises(KINSOL_Exception,KINSOL,self.pb_badconst)
        self.pb_badconst.set_x0([-1.0,-1.0,0.0,1.0,-1.0])
        nose.tools.assert_raises(KINSOL_Exception,KINSOL,self.pb_badconst)
        self.pb_badconst.set_x0([0.0,-1.0,0.0,1.0,1.0])
        nose.tools.assert_raises(KINSOL_Exception,KINSOL,self.pb_badconst)
        self.pb_badconst.set_x0([-1.0,-1.0,0.0,1.0,0.0])
        nose.tools.assert_raises(KINSOL_Exception,KINSOL,self.pb_badconst)
        
        self.pb_badconst.set_x0([-1.0,-1.0,-1.0,1.0,1.0])
        bad_solver = KINSOL(self.pb_badconst)
        self.pb_badconst.set_x0([-1.0,-1.0,1.0,1.0,1.0])
        bad_solver = KINSOL(self.pb_badconst)
        
    def test_sparsity_settings(self):
        """
        test if sparsity settings are set correctly
        """
        
        nose.tools.assert_raises(KINSOL_Exception,self.solve_pb_jac.set_sparsity,True)
        nose.tools.assert_raises(KINSOL_Exception,self.solve_pb_jac.set_sparsity,False)
        
        self.solve_pb_sparse.set_sparsity(True)
        nose.tools.assert_true(self.solve_pb_sparse.use_sparse)
        
        self.solve_pb_sparse.set_sparsity(False)
        nose.tools.assert_false(self.solve_pb_sparse.use_sparse)
        
        
        
        
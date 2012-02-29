# Author: Liwei Wang
"""
"""
from solvers import *
import numpy as np

class Pyds(Solver):
    '''
    Wrapper for ODE solvers in package PyDSTOOL.

    In this wrapper, we focus on solution of ODE problems. Other functionalities
    in PyDSTool (like analysis, continous trajectories) are neglected.
    
    If jacobian matrix is not supplied by users, PyDSTool.DIFF() would be used
    to estimate jacobian matrix approximately.
    
    This is an imcomplete version. Only Vode_ODESYSTEM(), which is the main
    solver in current PYDSTOOL, are implemented.    

    '''
    
    _optional_parameters = Solver._optional_parameters + ['jac', ]
    _name_pydstool = None

    def initialize(self):
        try:
            import PyDSTool
        except ImportError:
            raise ImportError,'''
        PyDSTool is not installed - required for solvers from PyDSTool'''

    def solve(self, time_points, terminate=None):
        # Common parts as superclass
        if terminate is None:
            terminate = lambda u, t, step_no: False

        self.t = np.asarray(time_points)
        self.initialize_for_solve()

        N = self.t.size - 1  # no of intervals
        self.validate_data()

        # As a main designing priciple of PyDSTool, most of data structures in
        # PyDSTool are index-free. That is, the numerical data are stored mainly
        # through Python dictionaries with string keys.
                
        # Start setting for PyDSTool
        import PyDSTool
        neq, f, u0 = self.neq, self.f, self.U0
        
        # Initialize variables as trajectories in PyDSTOOL
        # Each item of u has to be defined separately: y0,y1,y2,...
        name_list = ['y%d' % i for i in range(neq)]
        y, t, ydot = [PyDSTool.Var(name) for name in name_list], \
            PyDSTool.Var('t'), []
        
        # f has to be wrapped from f(y,t) to f(y0,y1,y2,...,t)
        f_wrap = eval('lambda *args: f(args[:-1], args[-1])')     
        #f_wrap = eval('lambda *args: f(args[:-1], args[-1])', locals()) # Error
        #f_wrap = lambda *args: f(args[:-1], args[-1])  # Error!

        # Define differiential functions in PyDSTOOL, item by item
        string2 = ','.join(['y[%d]' % i for i in range(neq)])  
        # y[0],y[1],y[2],...
        ydot = [eval('PyDSTool.Fun(f_wrap(%s,t)[%d],[%s],\"ydot%d\")' \
                         % (string2, i,string2,i)) for i in range(neq)]
        
        # Jacobian matrix 
        if getattr(self,'jac') is None:
            # apply Diff() to calculate jacobian matrix approximately.
            # Diff will return a QuantSpecct object
            F = eval('PyDSTool.Fun(f_wrap(%s,t),[%s],\"F\")' \
                         % (string2,string2))
            JAC = eval('PyDSTool.Fun(PyDSTool.Diff(F,[%s]),[t,%s],\"JAC\")' \
                         % (string2,string2))
        else:
            jac = self.jac
            # Wrap user-supplied jacobian function in the same manner as f
            #jac_wrap = lambda *args: jac(args[1:],args[0])    # Error
            jac_wrap = eval('lambda *args: jac(args[1:],args[0])')  
            JAC = eval('PyDSTool.Fun(jac_wrap(t,%s),[t,%s],\"JAC\")' \
                       % (string2,string2) )

        # Settings in PyDSTOOL
        DSargs = PyDSTool.args(name='pydstest',checklevel=2)
        # Function set is {JAC, ydot[0],ydot[1],...}
        string3 = ','.join(['ydot[%d]' % i for i in range(neq)])  
        DSargs.fnspecs = eval('[JAC,%s]' % string3)
        # Variable set is {y[i]:ydot[i](y[0],y[1]...)} for any i from 0 to neq-1
        string4 = ','.join(['y[%d]: ydot[%d](%s)' % (i, i, string2) \
                                for i in range(neq)])  
        DSargs.varspecs = eval('{%s}' % string4)
        # Time domain
        DSargs.tdomain = [time_points[0],time_points[-1]]
        # Initial status {y[0]:self.U0[0],y[1]:self.U0[1],...}
        string5 = ','.join(['y[%d]:%g' % (i, self.U0[i]) for i in range(neq)])  
        DSargs.ics = eval('{%s}' % string5)
        # Optional parameters
        DSargs.algparams = getattr(self,'params_pydstool',{})

        # Start computation
        test = eval(self._name_pydstool)(DSargs)
        result = test.compute('test','c')   # A trajectory returned
        # Extract and set values at points in time range
        self.u = np.asarray([[result(time)['y%d' % j] for j in range(neq)] \
                             for time in time_points])
        return self.u,self.t

class Vode_pyds(Pyds):
    '''
    Wrapper for Vode_ODESystem in PyDSTool.
    Applying vode.f. 

    Input parameters:       
       f         : Function to define the right side of equations. f(u,t)=u"(t)
       f_args    : Extra arguments for f.
       f_kwargs  : Extra keyword-arguments for f.
       complex_valued: Flag to indicate complex datatype.
       jac       : Function to define jacobian matrix.  
       init_step : Fixed step size for time mesh.,
       strictdt  : Boolean determining whether to evenly space time mesh
                   (default=False), or to use exactly dt spacing.
       stiff     : Boolean to activate the BDF method, otherwise Adams method
                   used. Default False.
       use_special: Switch for using special times,
       specialtimes: List of special times to use during iteration
       
    '''
    _optional_parameters = Pyds._optional_parameters + \
        ['init_step','strictdt','stiff','use_special','specialtimes']
    
    # Corresponding name in PyDSTool
    _name_pydstool = 'PyDSTool.Vode_ODEsystem'
    # Legal parameters in PyDSTool.Vode_ODEsystem
    _params_pydstool = ['init-step','strictdt','stiff','use_special',\
                       'specialtimes']

    def initialize_for_solve(self):
        if not hasattr(self,'init_step'):
            self.init_step = self.t[1] - self.t[0]
        self.params_pydstool = dict(\
            (key,value) for key,value in self.__dict__.items() \
            if key in self._params_pydstool)    

if __name__ == '__main__':

    f = lambda (u00,u11),t: [u11,-u00]
    method = Vode_pyds(f)
    method.set_initial_condition([0.,1.])
    u,t = method.solve(np.linspace(0.,10.,50))
    print u
    import scitools.std as st
    st.plot(t,u[:,0])
    print max(u[:,0]-np.sin(t))

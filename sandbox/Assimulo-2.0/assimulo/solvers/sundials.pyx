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
cimport numpy as N
from numpy cimport PyArray_DATA

N.import_array()

import numpy.linalg
import traceback 
 
from assimulo.exception import *

from assimulo.explicit_ode cimport Explicit_ODE
from assimulo.implicit_ode cimport Implicit_ODE

cimport sundials_includes as Sun

#Various C includes transfered to namespace
from sundials_includes cimport N_Vector, realtype, N_VectorContent_Serial, DENSE_COL
from sundials_includes cimport memcpy, N_VNew_Serial, DlsMat
from sundials_includes cimport malloc, free, realtype, N_VCloneVectorArray_Serial
from sundials_includes cimport N_VConst_Serial, N_VDestroy_Serial

include "constants.pxi" #Includes the constants (textual include)
include "sundials_constants.pxi" #Sundials related constants
include "sundials_callbacks.pxi"


cdef class IDA(Implicit_ODE):
    """
    This class provides a connection to the Sundials 
    (https://computation.llnl.gov/casc/sundials/main.html) solver IDA.
    
    IDA is a variable-order, variable-step multi-step algorithm for 
    solving differential algebraic equations of the form,
    
    .. math::
    
        F(t,y,\dot{y})=0, \quad y(t_0) = y_0, \quad \dot{y}(t_0) = \dot{y}_0.
    
    IDA includes the Backward Differentiation Formulas (BDFs).
    """
    cdef void* ida_mem
    cdef ProblemData pData      #A struct containing information about the problem
    cdef N_Vector yTemp, ydTemp
    cdef N_Vector *ySO
    cdef N_Vector *ydSO
    cdef object f
    #cdef public dict statistics
    cdef object pt_root, pt_fcn, pt_jac, pt_jacv, pt_sens
    cdef public N.ndarray yS0
    
    def __init__(self, problem):
        Implicit_ODE.__init__(self, problem) #Calls the base class
        
        self.pData = ProblemData()
        
        #Populate the ProblemData
        self.set_problem_data()
        
        #Solver options
        self.options["atol"] = N.array([1.0e-6]*self.problem_info["dim"])        #The absolute tolerance
        self.options["rtol"] = 1.0e-6        #The relative tolerance
        self.options["lsoff"] = False        #Turn on or off the linesearch algorithm
        self.options["tout1"] = 0.0001       #Direction of initial calculation
        self.options["suppress_alg"] = False #Turn on or off the local error test on algebraic variables
        self.options["suppress_sens"] = False #Turn on or off the local error test on the sensitivity variables
        self.options["maxsteps"] = 10000     #Maximum number of steps
        self.options["maxh"] = 0.0           #Maximum step-size
        self.options["maxord"] = 5           #Maximum order of method
        self.options["maxcorS"] = 3          #Maximum number of nonlinear iteration for sensitivity variables
        self.options["usejac"]   = True if self.problem_info["jac_fcn"] else False
        self.options["usesens"] = True if self.problem_info["dimSens"] > 0 else False
        self.options["inith"] = 0.0          #Initial step-size
        self.options["algvar"] = N.array([1.0]*self.problem_info["dim"])
        self.options["sensmethod"] = 'STAGGERED'
        self.options["dqtype"] = "CENTERED"
        self.options["dqrhomax"] = 0.0
        self.options["pbar"] = [1]*self.problem_info["dimSens"]
        
        #Statistics
        self.statistics["nfevals"]    = 0 #Function evaluations
        self.statistics["nsteps"]     = 0 #Number of steps
        self.statistics["netfails"]   = 0 #Number of error test failures
        self.statistics["nlinsetups"] = 0
        self.statistics["nncfails"]   = 0 #Nonlinear fails
        self.statistics["nniters"]    = 0 #Nonlinear iterations
        self.statistics["ngevals"]    = 0 #Root evaluations
        self.statistics["njevals"]    = 0 #Jacobian evaluations
        self.statistics["nfevalsLS"]  = 0 #Function evaluations due to Jac
        self.statistics["nfSevals"]   = 0 #Number of sensitivity evaluations
        self.statistics["nfevalsS"]   = 0 #Number of function evaluations due to finite sensitivity approximation
        self.statistics["nSetfails"]  = 0 #Number of error test failures
        self.statistics["nSniters"]   = 0 #Number of sensitivity nonlinear iterations
        self.statistics["nSncfails"]  = 0 #Number of sensitivity convergence failures
        
        #Solver support
        self.supports["one_step_mode"] = True
        self.supports["interpolated_output"] = True
        self.supports["interpolated_sensitivity_output"] = True
        self.supports["state_events"] = True
        
        #Get options from Problem
        if hasattr(problem, 'pbar'):
            self.pbar = problem.pbar
        elif hasattr(problem, 'p0'):
            self.pbar = N.array([N.abs(x) if N.abs(x) > 0 else 1.0 for x in self.p0])
        if hasattr(problem, 'algvar'):
            self.algvar = problem.algvar
        if hasattr(problem, 'yS0'):
            self.yS0 = problem.yS0
        
        
    cdef set_problem_data(self):
        #Sets the residual or rhs
        self.pt_fcn = self.problem.res
        self.pData.RHS = <void*>self.pt_fcn#<void*>self.problem.f
        self.pData.dim = self.problem_info["dim"] 
        self.pData.memSize = self.pData.dim*sizeof(realtype)
        
        #Set the ndarray to the problem struct
        #self.yTemp   = N.zeros(self.pData.dim, dtype=N.float, order='c')
        #self.ydTemp  = N.zeros(self.pData.dim, dtype=N.float, order='c')
        #self.pData.y  = <void*>self.yTemp
        #self.pData.yd = <void*>self.ydTemp 
        
        if self.problem_info["state_events"] is True: #Sets the root function
            self.pt_root = self.problem.state_events
            self.pData.ROOT = <void*>self.pt_root#<void*>self.problem.state_events
            self.pData.dimRoot = self.problem_info["dimRoot"]
            self.pData.memSizeRoot = self.pData.dimRoot*sizeof(realtype) 
    
        if self.problem_info["jac_fcn"] is True: #Sets the jacobian 
            self.pt_jac = self.problem.jac
            self.pData.JAC = <void*>self.pt_jac#<void*>self.problem.jac
            self.pData.memSizeJac = self.pData.dim*self.pData.dim*sizeof(realtype)
        
        if self.problem_info["jacv_fcn"] is True: #Sets the jacobian times vector
            self.pt_jacv = self.problem.jacv
            self.pData.JACV = <void*>self.pt_jacv#<void*>self.problem.jacv
            
        if self.problem_info["sens_fcn"] is True: #Sets the sensitivity function
            self.pt_sens = self.problem.sens
            self.pData.SENS = <void*>self.pt_sens#<void*>self.problem.sens
         
        if self.problem_info["dimSens"] > 0: #Sensitivity parameters (does not need the sensitivity function)
            self.pData.dimSens = self.problem_info["dimSens"]   
            self.pData.p    = <realtype*> malloc(self.problem_info["dimSens"]*sizeof(realtype))
            self.pData.pbar = <realtype*> malloc(self.problem_info["dimSens"]*sizeof(realtype))
        else:
            self.pData.dimSens = 0
        
        self.pData.verbose = 2
    
    def __dealloc__(self):
        
        if self.yTemp != NULL:
            #Deallocate N_Vector
            N_VDestroy_Serial(self.yTemp)
            
        if self.ydTemp != NULL:
            #Deallocate N_Vector
            N_VDestroy_Serial(self.ydTemp)
        
        if self.ida_mem != NULL: 
            #Free Memory
            Sun.IDAFree(&self.ida_mem)
    
    cpdef state_event_info(self):
        """
        Returns the event info.
        """
        cdef int* c_info 
        cdef flag
        
        # Allocate memory for the event_info_ vector and initialize to zeros
        c_info = <int*> malloc(self.pData.dimRoot*sizeof(int))
        for k in range(self.pData.dimRoot):
            c_info[k] = 0
        
        # Fetch data on which root functions that became zero and store in class
        flag = Sun.IDAGetRootInfo(self.ida_mem, c_info)
        if flag < 0:
            raise IDAError(flag)
        
        #event_info = PyArray_SimpleNew(1,&self.pData.dimRoot,NPY_INT)
        event_info = [0]*self.pData.dimRoot
        
        for k in range(self.pData.dimRoot):
            event_info[k] = c_info[k]
        
        # Remember to deallocate
        free(c_info)
        
        return event_info
    
    cpdef initialize(self):
        
        #Initialize storing of sensitivyt result in handle_result
        if self.problem_info['step_events'] or self.options['continuous_output']:
            self.problem._sensitivity_result = 1
        
        #Reset statistics
        for k in self.statistics.keys():
            self.statistics[k] = 0
        
        self.initialize_ida()
    
    cdef initialize_ida(self):
        cdef int flag #Used for return
        cdef realtype ZERO = 0.0

        self.yTemp  = arr2nv(self.y)
        self.ydTemp = arr2nv(self.yd)
        
        #Updates the switches
        if self.problem_info["switches"]:
            self.pData.sw = <void*>self.sw
        
        if self.pData.dimSens > 0:
            #Create the initial matrices
            self.ySO  = N_VCloneVectorArray_Serial(self.pData.dimSens, self.yTemp)
            self.ydSO = N_VCloneVectorArray_Serial(self.pData.dimSens, self.ydTemp)
            
            #Filling the start vectors
            for i in range(self.pData.dimSens):
                 N_VConst_Serial(ZERO,  self.ySO[i]);
                 N_VConst_Serial(ZERO, self.ydSO[i]); 
                 if self.yS0 != None:
                    for j in range(self.pData.dim):
                        (<N_VectorContent_Serial>self.ySO[i].content).data[j] = self.yS0[i,j]

        if self.ida_mem == NULL: #The solver is not initialized
        
            self.ida_mem = Sun.IDACreate() #Create solver
            if self.ida_mem == NULL:
                raise IDAError(IDA_MEM_FAIL)
            
            #Specify the residual and the initial conditions to the solver
            flag = Sun.IDAInit(self.ida_mem, ida_res, self.t, self.yTemp, self.ydTemp)
            if flag < 0:
                raise IDAError(flag, self.t)
                
            #Specify the use of the internal dense linear algebra functions.
            flag = Sun.IDADense(self.ida_mem, self.pData.dim)
            if flag < 0:
                raise IDAError(flag, self.t)
            
            #Specify the root function to the solver
            if self.pData.ROOT != NULL:
                flag = Sun.IDARootInit(self.ida_mem, self.pData.dimRoot, ida_root)
                if flag < 0:
                    raise IDAError(flag,self.t)
            
            #Specify the error handling
            flag = Sun.IDASetErrHandlerFn(self.ida_mem, ida_err, <void*>self.pData)
            if flag < 0:
                raise IDAError(flag, self.t)
                
            if self.pData.dimSens > 0:
                flag = Sun.IDASensInit(self.ida_mem, self.pData.dimSens, IDA_STAGGERED if self.options["sensmethod"] == "STAGGERED" else IDA_SIMULTANEOUS, NULL, self.ySO, self.ydSO)
                if flag < 0:
                    raise IDAError(flag, self.t)
            
        else: #The solver needs to be reinitialized
            
            #Reinitialize
            flag = Sun.IDAReInit(self.ida_mem, self.t, self.yTemp, self.ydTemp)
            if flag < 0:
                raise IDAError(flag, self.t)
                
            if self.pData.dimSens > 0:
                flag = Sun.IDASensReInit(self.ida_mem, IDA_STAGGERED if self.options["sensmethod"] == "STAGGERED" else IDA_SIMULTANEOUS, self.ySO, self.ydSO)
                if flag < 0:
                    raise IDAError(flag, self.t)
        
        #Specify the jacobian to the solver
        if self.pData.JAC != NULL and self.options["usejac"]:
            
            flag = Sun.IDADlsSetDenseJacFn(self.ida_mem, ida_jac)
            if flag < 0:
                raise IDAError(flag,t0)
        else:
            flag = Sun.IDADlsSetDenseJacFn(self.ida_mem, NULL)
            if flag < 0:
                raise IDAError(flag,t0)
        
        #Set the user data
        flag = Sun.IDASetUserData(self.ida_mem, <void*>self.pData)
        if flag < 0:
            raise IDAError(flag, self.t)
    
    cdef initialize_sensitivity_options(self):
        """
        Sets the sensitivity information.
        """
        cdef int flag
        
        #Sets the parameters to the userdata object.
        for i in range(self.pData.dimSens):
            self.pData.p[i] = self.p[i]
        #Sets the pbar to the userdata object.
        for i in range(self.pData.dimSens):
            self.pData.pbar[i] = self.options["pbar"][i]
        
        #Specify problem parameter information for sensitivity calculations
        flag = Sun.IDASetSensParams(self.ida_mem, self.pData.p, self.pData.pbar, NULL)
        if flag < 0:
            raise IDAError(flag, self.t)
        
        #Specify the difference quotient strategy
        flag = Sun.IDASetSensDQMethod(self.ida_mem, IDA_CENTERED if self.options["dqtype"]=="CENTERED" else IDA_FORWARD, self.options["dqrhomax"])
        if flag<0:
            raise IDAError(flag, self.t)
        
        #Specify the error control strategy
        flag = Sun.IDASetSensErrCon(self.ida_mem, self.options["suppress_sens"]==False)
        if flag < 0:
            raise IDAError(flag, self.t)
        
        #Specify the maximum number of nonlinear solver iterations
        flag = Sun.IDASetSensMaxNonlinIters(self.ida_mem, self.options["maxcorS"])
        if flag < 0:
            raise IDAError(flag, self.t)
        
        #Estimate the sensitivity  ----SHOULD BE IMPROVED with IDASensSVTolerances ...
        flag = Sun.IDASensEEtolerances(self.ida_mem)
        if flag < 0:
            raise IDAError(flag, self.t)
        
        #Should the sensitivities be calculated this time around?
        if self.options["usesens"] == False:
            flag = Sun.IDASensToggleOff(self.ida_mem)
            if flag < 0:
                raise IDAError(flag, self.t)
    
    
    cdef initialize_options(self):
        """
        Updates the simulation options.
        """
        cdef flag
        
        #Maximum order
        flag = Sun.IDASetMaxOrd(self.ida_mem, self.options["maxord"])
        if flag < 0:
            raise IDAError(flag)
            
        #Initial step
        flag = Sun.IDASetInitStep(self.ida_mem, self.options["inith"])
        if flag < 0:
            raise IDAError(flag)
            
        #Maximum step
        flag = Sun.IDASetMaxStep(self.ida_mem, self.options["maxh"])
        if flag < 0:
            raise IDAError(flag)
            
        #Maximum Number of steps
        flag = Sun.IDASetMaxNumSteps(self.ida_mem, self.options["maxsteps"])
        if flag < 0:
            raise IDAError(flag)
        
        #Set the algebraic components and the differential
        flag = Sun.IDASetId(self.ida_mem, arr2nv(self.options["algvar"]))
        if flag < 0:
            raise IDAError(flag)
        
        #Suppress algebraic components on the error test
        flag = Sun.IDASetSuppressAlg(self.ida_mem, self.options["suppress_alg"])
        if flag < 0:
            raise IDAError(flag)
            
        #Set the tolerances
        flag = Sun.IDASVtolerances(self.ida_mem, self.options["rtol"], arr2nv(self.options["atol"]))
        if flag < 0:
            raise IDAError(flag)
            
        #Initialize sensitivity if any
        if self.pData.dimSens > 0:
            self.initialize_sensitivity_options()
    
    cpdef integrate(self,double t,N.ndarray[ndim=1, dtype=realtype] y,N.ndarray[ndim=1, dtype=realtype] yd,double tf,dict opts):
        cdef int flag, output_index, normal_mode
        cdef N_Vector yout, ydout
        cdef double tret = 0.0, tout
        cdef list tr = [], yr = [], ydr = []
        cdef N.ndarray output_list
        yout = arr2nv(y)
        ydout = arr2nv(yd)
        
        #Get options
        initialize   = opts["initialize"]
        
        #Initialize? 
        if initialize:
            self.initialize_ida()
            self.initialize_options()
        
        #Set stop time
        flag = Sun.IDASetStopTime(self.ida_mem, tf)
        if flag < 0:
            raise IDAError(flag, t)
        
        #Run in normal mode?
        normal_mode = 1 if opts["output_list"] != None else 0
        
        if normal_mode == 0: 
            #Integration loop
            while True:
                
                #Integration loop
                flag = Sun.IDASolve(self.ida_mem,tf,&tret,yout,ydout,IDA_ONE_STEP)
                if flag < 0:
                    raise IDAError(flag, tret)
                
                #Store results
                tr.append(tret)
                yr.append(nv2arr(yout))
                ydr.append(nv2arr(ydout))
                
                if flag == IDA_ROOT_RETURN: #Found a root
                    flag = ID_EVENT #Convert to Assimulo flags
                    self.store_statistics()
                    break
                    
                if flag == IDA_TSTOP_RETURN: #Reached tf
                    flag = ID_COMPLETE
                    self.store_statistics()
                    break
        else:
            output_index = opts["output_index"]
            output_list  = opts["output_list"][output_index:]
            
            for tout in output_list:
                output_index += 1
                #Integration loop
                flag = Sun.IDASolve(self.ida_mem,tout,&tret,yout,ydout,IDA_NORMAL)
                if flag < 0:
                    raise IDAError(flag, tret)
                
                #Store results
                tr.append(tret)
                yr.append(nv2arr(yout))
                ydr.append(nv2arr(ydout))
                
                if flag == IDA_ROOT_RETURN: #Found a root
                    flag = ID_EVENT #Convert to Assimulo flags
                    self.store_statistics()
                    break
                    
                if flag == IDA_TSTOP_RETURN: #Reached tf
                    flag = ID_COMPLETE
                    self.store_statistics()
                    break
            else:
                flag = ID_COMPLETE
                self.store_statistics()
            
            opts["output_index"] = output_index
        
        #Deallocate
        N_VDestroy_Serial(yout)
        N_VDestroy_Serial(ydout)
        
        return flag, tr, yr, ydr
    
    
    cpdef step(self,double t,N.ndarray y, N.ndarray yd, double tf,dict opts):
        cdef int flag
        cdef N_Vector yout, ydout
        cdef double tret = 0.0
        cdef double tr
        cdef N.ndarray yr, ydr
        
        yout  = arr2nv(y)
        ydout = arr2nv(yd)
        
        #Get options
        initialize  = opts["initialize"]
        
        #Initialize?
        if initialize:
            self.initialize_ida()
            self.initialize_options()
        
        #Set stop time
        flag = Sun.IDASetStopTime(self.ida_mem, tf)
        if flag < 0:
            raise IDAError(flag, t)
        
        #Integration loop
        flag = Sun.IDASolve(self.ida_mem,tf,&tret,yout,ydout,IDA_ONE_STEP)
        if flag < 0:
            raise IDAError(flag, tret)
            
        #Store results
        tr  = tret
        yr  = nv2arr(yout)
        ydr = nv2arr(ydout)
        
        if flag == IDA_ROOT_RETURN: #Found a root
            flag = ID_EVENT #Convert to Assimulo flags
            self.store_statistics()
            
        if flag == IDA_TSTOP_RETURN: #Reached tf
            flag = ID_COMPLETE
            self.store_statistics()
        
        #Deallocate
        N_VDestroy_Serial(yout)
        N_VDestroy_Serial(ydout)
        
        return flag, tr, yr, ydr
    
    cpdef make_consistent(self, method):
        """
        Directs IDA to try to calculate consistent initial conditions.
            
            Parameters::
            
                method  
                        - 'IDA_YA_YDP_INIT'
                                - This tries to calculate the
                                  algebraic components of y and the differential
                                  components of yd given the differential components
                                  of y. The algebraic components of y must have been
                                  specified with the property 'algvar'. The property
                                  'IDA.tout1' is  used in the calculations. It 
                                  represents the the next output point.
                                
                        - 'IDA_Y_INIT' 
                                - This tries to calculate all components
                                  of y given yd.
            
        See SUNDIALS IDA documentation 4.5.4 for more details.
        """
        cdef double direction
        cdef int icopt, flag
        
        self.initialize()
        self.initialize_options()
        
        #Determine method
        if method == 'IDA_Y_INIT':
            icopt = IDA_Y_INIT
        elif method == 'IDA_YA_YDP_INIT':
            icopt = IDA_YA_YDP_INIT 
        else:
            raise Exception("The method is unknown.")
        
        direction = self.t+self.options["tout1"] #tout1 is needed for the solver to determine the direction of the integration
        
        if self.ida_mem == NULL: 
            raise Exception("IDA must be initialized.")
        
        #Set the options lsoff and calculate initial conditions
        flag = Sun.IDASetLineSearchOffIC(self.ida_mem, self.options["lsoff"])
        flag = Sun.IDACalcIC(self.ida_mem, icopt, direction)
        
        if flag == IDA_SUCCESS: #Gets the calculated values
            flag = Sun.IDAGetConsistentIC(self.ida_mem, self.yTemp, self.ydTemp)
        
        #Set the calculated values to the current ones
        self.y  = nv2arr(self.yTemp)
        self.yd = nv2arr(self.ydTemp)
        
        return [flag, self.y, self.yd]
    
    cpdef N.ndarray interpolate(self,double t,int k = 0):
        """
        Calls the internal IDAGetDky for the interpolated values at time t.
        t must be within the last internal step. k is the derivative of y which
        can be from zero to the current order.
        """
        cdef flag
        cdef N.ndarray res
        cdef N_Vector dky=N_VNew_Serial(self.pData.dim)
        
        flag = Sun.IDAGetDky(self.ida_mem, t, k, dky)
        
        if flag < 0:
            raise IDAError(flag, t)
        
        res = nv2arr(dky)
        
        N_VDestroy_Serial(dky) #Deallocate
        
        return res
        
    cpdef interpolate_sensitivity(self,double t, int k = 0, int i=-1):
        """
        This method calls the internal method IDAGetSensDky which computes the k-th derivatives
        of the interpolating polynomials for the sensitivity variables at time t.
        
            Parameters::
                    
                    t
                        - Specifies the time at which sensitivity information is requested. The time
                          t must fall within the interval defined by the last successful step taken
                          by IDAS.
                    
                    k   
                        - The order of derivatives.
                        
                    i
                        - Specifies the sensitivity derivative vector to be returned (0<=i<=Ns)
                        
            Return::
            
                    A matrix containing the Ns vectors or a vector if i is specified.
        """
        cdef N_Vector dkyS=N_VNew_Serial(self.pData.dim)
        cdef flag
        cdef N.ndarray res
        
        if i==-1:
            
            matrix = []
            
            for x in xrange(self.pData.dimSens):
                flag = Sun.IDAGetSensDky1(self.ida_mem, t, k, x, dkyS)
                
                if flag<0:
                    raise IDAError(flag, t)
                
                matrix += [nv2arr(dkyS)]
            
            N_VDestroy_Serial(dkyS)
            
            return np.array(matrix)
        else:
            flag = Sun.IDAGetSensDky1(self.ida_mem, t, k, i, dkyS)
            
            if flag <0:
                raise IDAError(flag, t)
            
            res = nv2arr(dkyS)
            
            N_VDestroy_Serial(dkyS)
            
            return res
            
    def _set_lsoff(self, lsoff):
        try:
            self.options["lsoff"] = bool(lsoff)
        except:
            Exception('Unkown input to lsoff, must be a boolean.')
        
    def _get_lsoff(self):
        """
        Boolean value to turn OFF Sundials LineSearch when calculating
        initial conditions.
            
            Parameters::
            
                lsoff   
                        - Default 'False'. False indicates the use of
                          linesearch.

                        - Should be a boolean.
                        
                            Example:
                                lsoff = True
        """
        return self.options["lsoff"]
    
    lsoff = property(_get_lsoff, _set_lsoff)
    
    def _set_initial_step(self, initstep):
        try:
            self.options["inith"] = float(initstep)
        except (ValueError, TypeError):
            raise Exception('The initial step must be an integer or float.')
        
    def _get_initial_step(self):
        """
        This determines the initial step-size to be used in the integration.
        
            Parameters::
            
                initstep    
                            - Default '0.0', which result in that the
                              the initial step is approximated.
                            
                            - Should be float.
                            
                                Example:
                                    initstep = 0.01
        """
        return self.options["inith"]
        
    inith = property(_get_initial_step,_set_initial_step)
    
    def _set_tout1(self, tout1):
        try:
            self.options["tout1"] = float(tout1)
        except (ValueError,TypeError):
            raise Exception('tout1 must be an integer or float.')
        
    def _get_tout1(self):
        """
        Sets the value used in the internal Sundials function
        for determine initial conditions. This value is needed
        in order to manipulate the existing machinery to be used
        in determining the initial conditions.
        
            Parameters::
            
                tout1       
                            - Default '0.001'.
                            
                            - Should be a float.
                            
                                Example:
                                    tout1 = 0.01
        """
        return self.options["tout1"]

    tout1=property(_get_tout1,_set_tout1)
    
    def _set_suppress_alg(self,suppress_alg):
        try:
            self.options["suppress_alg"] = bool(suppress_alg)
        except:
            raise Exception("Unkown input to suppress_alg, must be a boolean.")

    def _get_suppress_alg(self):
        """
        A Boolean flag which indicates that the error-tests are 
        suppressed on algebraic variables. The algebraic variables
        are defined by setting the property 'algvar'.
        
            Parameters::
            
                suppress_alg    
                                - Default 'False'.
                
                                - Should be a boolean.
                                
                                    Example:
                                        suppress_alg = True
                                        
        See SUNDIALS IDA documentation 4.5.7 'IDASetSuppressAlg' 
        for more details.
        """
        return self.options["suppress_alg"]    

    suppress_alg=property(_get_suppress_alg,_set_suppress_alg)
    
    def _set_suppress_sens(self,suppress_sens):
        try:
            self.options["suppress_sens"] = bool(suppress_sens)
        except:
            raise Exception("Unkown input to suppress_sens, must be a boolean.")

    def _get_suppress_sens(self):
        """
        A Boolean flag which indicates that the error-tests are 
        suppressed on the sensitivity variables.
        
            Parameters::
            
                suppress_sens    
                                - Default 'False'.
                
                                - Should be a boolean.
                                
                                    Example:
                                        suppress_alg = True
                                        
        See SUNDIALS IDA documentation 5.2.6 'IDASetSensErrCon' 
        for more details.
        """
        return self.options["suppress_sens"]    

    suppress_sens=property(_get_suppress_sens,_set_suppress_sens)
    
    def _set_atol(self,atol):
        
        self.options["atol"] = N.array(atol,dtype=N.float) if len(N.array(atol,dtype=N.float).shape)>0 else N.array([atol],dtype=N.float)
    
        if len(self.options["atol"]) == 1:
            self.options["atol"] = self.options["atol"]*N.ones(self.pData.dim)
        elif len(self.options["atol"]) != self.pData.dim:
            raise Exception("atol must be of length one or same as the dimension of the problem.")
        if (self.options["atol"]<0.0).any():
            raise Exception("atol must be positive.")
    
    def _get_atol(self):
        """
        Defines the absolute tolerance(s) that is to be used by the solver.
        Can be set differently for each variable.
        
            Parameters::
            
                atol    
                        - Default '1.0e-6'.
                
                        - Should be a positive float or a numpy vector
                          of floats.
                        
                            Example:
                                atol = [1.0e-4, 1.0e-6]
        
        See SUNDIALS IDA documentation 4.5.2 for more details.
        """
        return self.options["atol"]
    
    atol=property(_get_atol,_set_atol)
    
    def _set_rtol(self,rtol):
        
        try:
            rtol = float(rtol)
        except (ValueError, TypeError):
            raise Exception('Relative tolerance must be a (scalar) float.')
        if rtol <= 0.0:
            raise Exception('Relative tolerance must be a positive (scalar) float.')
        
        self.options["rtol"] = rtol
    
    def _get_rtol(self):
        """
        Defines the relative tolerance that is to be used by the solver.
        
            Parameters::
            
                rtol    
                        - Default '1.0e-6'.
                
                        - Should be a positive float.
                        
                            Example:
                                rtol = 1.0e-4
                                
        """
        return self.options["rtol"]
        
    rtol=property(_get_rtol,_set_rtol)
    
    def _set_max_ord(self,maxord):
        try:
            maxord = int(maxord)
        except ValueError:
            raise Exception("The maximal order must be an integer.")
        
        if maxord > 5 or maxord < 1:
            self.options["maxord"] = 5 if maxord > 5 else 1
            self.log_message("The maximal order must be between 1-5.", NORMAL)
        else:
            self.options["maxord"] = maxord
            
    
    def _get_max_ord(self):
        """
        This determines the maximal order that is be used by the solver.
        
            Parameters::
            
                maxord  
                        - For the BDF method the maximum order is 5. 
                          'maxord' can be set in an interval from 1 to 
                          the maximum order allowed.
                
                        - Should be an integer.
                        
                            Example:
                                maxord = 3
    
        
        An input value greater than the maximal order will result in the 
        maximum value.
        """
        return self.options["maxord"]

    maxord=property(_get_max_ord,_set_max_ord)
    
    def _set_max_cor_S(self,maxcorS):
        try:
            self.options["maxcorS"] = int(maxcorS)
        except:
            raise Exception("The maximum nonlinear sensitivity iterations must be a positiv integer.")
    
    def _get_max_cor_S(self):
        """
        This detmines the maximum number of nonlinear iterations for the
        sensitivity variables.
        
            Parameters::
            
                maxcorS
                        - Default 3
                        
                        - Should be an integer
                        
        For more information see SUNDIALS IDAS documentation 5.2.6.
        """
        return self.options["maxcorS"]
    
    maxcorS=property(_get_max_cor_S,_set_max_cor_S)
    
    def _set_max_steps(self, maxsteps):
        if not isinstance(maxsteps,int):
            raise Exception('The maximum number of steps must be an integer.')
        if maxsteps < 1:
            raise Exception('The maximum number of steps must be a positive integer.')
        self.options["maxsteps"] = maxsteps
    
    def _get_max_steps(self):
        """
        Determines the maximum number of steps the solver is allowed
        to take to finish the simulation.
        
            Parameters::
            
                maxsteps    
                            - Default '10000'.
                            
                            - Should be an integer.
                            
                                Example:
                                    maxsteps = 1000.
                                    
        """
        return self.options["maxsteps"]

    maxsteps = property(_get_max_steps, _set_max_steps)
    
    def _set_max_h(self,max_h):
        try:
            self.options["maxh"] = float(max_h)
        except:
            raise Exception("Maximal stepsize must be a (scalar) float.")
    
    def _get_max_h(self):
        """
        Defines the maximal step-size that is to be used by the solver.
        
            Parameters::
            
                maxh    
                        - Default '0', which indicates that the maximal
                          step-size is infinity.
                          
                        - Should be a float.
                        
                            Example:
                                maxh = 0.01
                                
        """
        return self.options["maxh"]
        
    maxh=property(_get_max_h,_set_max_h)
    
    def _set_algvar(self,algvar):
        self.options["algvar"] = N.array(algvar,dtype=N.float) if len(N.array(algvar,dtype=N.float).shape)>0 else N.array([algvar],dtype=N.float)
        
        if len(self.options["algvar"]) != self.pData.dim:
            raise Exception('When setting the algebraic variables, the' \
                                ' vector must be of the same size as the problem dimension.')
        
    def _get_algvar(self):
        """
        A list for defining which variables are differential and
        which are algebraic.
        This list is used when excluding algebraic variables from the error test
        by setting suppress_alg=True  and it is used, when computing consistent initial 
        values using the method make_consistency
        
            Parameters::
            
                algvar  
                        - The value True(1.0) indicates a differential
                          variable and the value False(0.0) indicates an
                          algebraic variable.
                          
                        - Should be a list or a numpy vector (ndarray)
                        
                            Example:
                                algvar = [1.0, 0.0, 1.0]
                                algvar = [True, False, True]
                                algvar = [1,0,1]
                                
        """
        return self.options["algvar"]    

    algvar=property(_get_algvar,_set_algvar)
    
    def _set_usejac(self, jac):
        self.options["usejac"] = bool(jac)
    
    def _get_usejac(self):
        """
        This sets the option to use the user defined jacobian. If a
        user provided jacobian is implemented into the problem the
        default setting is to use that jacobian. If not, an
        approximation is used.
        
            Parameters::
            
                usejac  
                        - True - use user defined jacobian
                          False - use an approximation
                    
                        - Should be a boolean.
                        
                            Example:
                                usejac = False
        """
        return self.options["usejac"]
    
    usejac = property(_get_usejac,_set_usejac)
    
    def _set_usesens(self, sens):
        self.options["usesens"] = bool(sens)
    
    def _get_usesens(self):
        """
        This options activates or deactivates the sensitivity 
        calculations.
        
            Parameters::
            
                usejac  
                        - True -  Activate sensitivity calculations
                          False - Deactivate sensitivity calculations
                    
                        - Should be a boolean.
                        
                            Example:
                                usesens = False
        """
        return self.options["usesens"]
    
    usesens = property(_get_usesens,_set_usesens)
    
    def _set_sensitivity_method(self, ism):
        if not isinstance(ism, str):
            raise Exception('Sensitivity method must be string.')
        
        if ism.upper() == 'SIMULTANEOUS':
            self.options["sensmethod"] = 'SIMULTANEOUS' 
        elif ism.upper() == 'STAGGERED':
            self.options["sensmethod"] = 'STAGGERED'
        else:
            raise Exception('Sensitivity method must be either "SIMULTANEOUS" or "STAGGERED".')
        
    def _get_sensitivity_method(self):
        """
        Specifies the sensitivity solution method. Can be either,
        'SIMULTANEOUS' or 'STAGGERED'.
            
        Parameters::
            
            ism
                    - A string of either 'SIMULTANEOUS' or 'STAGGERED'
                    - Default 'STAGGERED'
                        
        Returns::
            
            The current value of sensmethod (string)
        
        See SUNDIALS documentation '(IDA/CVode)SensInit'
        """
        return self.options["sensmethod"]
        
    sensmethod = property(_get_sensitivity_method, _set_sensitivity_method)
    
    def _set_dqtype(self, dqtype):
        if not isinstance(dqtype, str):
            raise Exception('DQtype must be string.')
        
        if dqtype.upper() == 'CENTERED':
            self.options["dqtype"] = "CENTERED"
        elif dqtype.upper() == 'FORWARD':
            self.options["dqtype"] = "FORWARD"
        else:
            raise Exception('DQtype must be either "CENTERED" or "FORWARD".')
            
    def _get_dqtype(self):
        """
        Specifies the difference quotient type in the sensitivity calculations.
        Can be either, 'CENTERED' or 'FORWARD'.
        
        Parameters::
            
            DQtype 
                    - A string of either 'CENTERED' or 'FORWARD'
                    - Default 'CENTERED'
                        
        Returns::
            
            The current value of DQtype.
        
        See SUNDIALS documentation '(IDA/CVode)SetSensDQMethod' 
        """
        return self.options["dqtype"]
    
    dqtype = property(_get_dqtype, _set_dqtype)
    
    def _set_dqrhomax(self, dqrhomax):
        try:
            self.options["dqrhomax"] = float(dqrhomax)
        except (TypeError, ValueError):
            raise Exception('DQrhomax must be convertable to a float.')
        if self.options["dqrhomax"] < 0.0:
            raise Exception("DQrhomax must be positive.")
            
    def _get_dqrhomax(self):
        """
        Specifies the selection parameters used in deciding switching between a simultaneous
        or separate approximation of the two terms in the sensitivity residual.
        
            Parameters::
            
                DQrhomax
                        - A postive float.
                        - Default 0.0
                        
            Returns::
            
                The current value of DQrhomax (float)
        
        See SUNDIALS documentation '(IDA/CVode)SetSensDQMethod' 
        """
        return self.options["dqrhomax"]
    
    dqrhomax = property(_get_dqrhomax, _set_dqrhomax)
    
    def _set_pbar(self, pbar):
        if len(pbar) != self.problem_info['dimSens']:
            raise Exception('pbar must be of equal length as the parameters.')
        
        self.options["pbar"] = pbar
    
    def _get_pbar(self):
        """
        Specifies the order of magnitude for the parameters. This is useful if IDAS is
        to estimate tolerances for the sensitivity solution vectors.
        
            Parameters::
            
                pbar
                        - An array of positive floats equal to the number of parameters.
                        - Default absolute values of the parameters.
                        
            Returns::
            
                The current value of pbar.
                
        See SUNDIALS documentation '(IDA/CVode)SetSensParams'
        """
        return self.options["pbar"]
    
    pbar = property(_get_pbar, _set_pbar)
    
    cdef void store_statistics(self):
        """
        Retrieves and stores the statistics.
        """
        cdef long int nsteps = 0, nrevals = 0, nlinsetups = 0, netfails = 0
        cdef long int nniters = 0, nncfails = 0, ngevals = 0
        cdef long int nSniters = 0, nSncfails = 0, njevals = 0, nrevalsLS = 0
        cdef long int nfSevals = 0, nfevalsS = 0, nSetfails = 0, nlinsetupsS = 0
        cdef int klast, kcur
        cdef realtype hinused, hlast, hcur, tcur
        
        flag = Sun.IDAGetIntegratorStats(self.ida_mem, &nsteps, &nrevals, &nlinsetups, &netfails,
                                         &klast, &kcur, &hinused, &hlast, &hcur, &tcur)
        flag = Sun.IDAGetNonlinSolvStats(self.ida_mem, &nniters, &nncfails)
        flag = Sun.IDAGetNumGEvals(self.ida_mem, &ngevals)
        flag = Sun.IDADlsGetNumJacEvals(self.ida_mem, &njevals)
        flag = Sun.IDADlsGetNumResEvals(self.ida_mem, &nrevalsLS)
        
        self.statistics["nsteps"] += nsteps
        self.statistics["nfevals"] += nrevals
        self.statistics["netfails"] += netfails
        self.statistics["nniters"] += nniters
        self.statistics["nncfails"] += nncfails
        self.statistics["ngevals"] += ngevals
        self.statistics["njevals"] += njevals
        self.statistics["nfevalsLS"] += nrevalsLS
        
        #If sensitivity    
        if self.pData.dimSens > 0:
            flag = Sun.IDAGetSensStats(self.ida_mem, &nfSevals, &nfevalsS, &nSetfails, &nlinsetupsS)
            flag = Sun.IDAGetSensNonlinSolvStats(self.ida_mem, &nSniters, &nSncfails)
            
            self.statistics["nfSevals"]   += nfSevals
            self.statistics["nfevalsS"]   += nfevalsS
            self.statistics["nSetfails"]  += nSetfails
            self.statistics["nSniters"]   += nSniters
            self.statistics["nSncfails"]  += nSncfails
    
    def print_statistics(self, verbose=NORMAL):
        """
        Prints the run-time statistics for the problem.
        """
        self.log_message('Final Run Statistics: %s \n' % self.problem.name,        verbose)
        
        self.log_message(' Number of Steps                          : '+ str(self.statistics["nsteps"]),         verbose)               
        self.log_message(' Number of Function Evaluations           : '+ str(self.statistics["nfevals"]),        verbose)
        self.log_message(' Number of Jacobian Evaluations           : '+ str(self.statistics["njevals"]),        verbose)
        self.log_message(' Number of F-Eval During Jac-Eval         : '+ str(self.statistics["nfevalsLS"]),      verbose)
        self.log_message(' Number of Root Evaluations               : '+ str(self.statistics["ngevals"]),        verbose)
        self.log_message(' Number of Error Test Failures            : '+ str(self.statistics["netfails"]),       verbose)
        self.log_message(' Number of Newton Iterations              : '+ str(self.statistics["nniters"]),        verbose)
        self.log_message(' Number of Newton Convergence Failures    : '+ str(self.statistics["nncfails"]),       verbose)
        
        if self.problem_info['dimSens'] > 0: #Senstivity calculations is on
            self.log_message('\nSensitivity Statistics:\n', verbose)
            self.log_message(' Number of Sensitivity Calculations             : ' + str(self.statistics["nfSevals"]), verbose)
            self.log_message(' Number of F-Evals Due to Finite Approximation  : ' + str(self.statistics["nfevalsS"]), verbose)
            self.log_message(' Number of Local Error Test Failures            : ' + str(self.statistics["nSetfails"]),verbose)
            self.log_message(' Number of Newton Iterations                    : ' + str(self.statistics["nSniters"]), verbose)
            self.log_message(' Number of Newton Convergance Failures          : ' + str(self.statistics["nSncfails"]),verbose)
            
            self.log_message('\nSensitivity options:\n' , verbose)
            self.log_message(' Method                   : ' + str(self.options["sensmethod"]), verbose)
            self.log_message(' Difference quotient type : ' + str(self.options["dqtype"]), verbose)
            self.log_message(' Suppress Sens            : ' + str(self.options["suppress_sens"]), verbose)
   
        self.log_message('\nSolver options:\n',                                       verbose)
        self.log_message(' Solver                  : IDA (BDF)',                      verbose)
        self.log_message(' Maxord                  : ' + str(self.options["maxord"]), verbose)
        self.log_message(' Suppress Alg            : ' + str(self.options["suppress_alg"]), verbose)
        self.log_message(' Tolerances (absolute)   : ' + str(self.options["atol"]),   verbose)
        self.log_message(' Tolerances (relative)   : ' + str(self.options["rtol"]),   verbose)
        self.log_message('',                                                          verbose)

cdef class CVode(Explicit_ODE):
    """
    This class provides a connection to the Sundials 
    (https://computation.llnl.gov/casc/sundials/main.html) solver CVode.
    
    CVode is a variable-order, variable-step multi-step algorithm for 
    solving ordinary differential equations of the form,
    
    .. math::
    
        \dot{y} = f(t,y), \quad y(t_0) = y_0.
    
    CVode includes the Backward Differentiation Formulas (BDFs) which 
    are suitable for stiff problems and also the Adams-Moulton formulas 
    for non-stiff systems.
    """
    cdef void* cvode_mem
    cdef ProblemData pData      #A struct containing information about the problem
    cdef N_Vector yTemp, ydTemp
    cdef N_Vector *ySO
    cdef object f
    #cdef public dict statistics
    cdef object pt_root, pt_fcn, pt_jac, pt_jacv, pt_sens
    cdef public N.ndarray yS0
    
    def __init__(self, problem):
        Explicit_ODE.__init__(self, problem) #Calls the base class

        self.pData = ProblemData()
        
        #Populate the ProblemData
        self.set_problem_data()
        
        #Solver options
        self.options["atol"] = N.array([1.0e-6]*self.problem_info["dim"])        #The absolute tolerance
        self.options["rtol"] = 1.0e-6        #The relative tolerance
        self.options["maxh"] = 0.0           #Maximum step-size
        self.options["inith"] = 0.0          #Initial step-size
        self.options["maxord"] = 5        #Maximum order allowed
        self.options["usejac"]   = True if (self.problem_info["jac_fcn"] or self.problem_info["jacv_fcn"]) else False
        self.options["usesens"] = True if self.problem_info["dimSens"] > 0 else False
        self.options["maxsteps"] = 10000  #Maximum number of steps
        self.options["sensmethod"] = 'STAGGERED'
        self.options["linear_solver"] = "DENSE"
        self.options["iter"] = "Newton"
        self.options["discr"] = "BDF"
        self.options["suppress_sens"] = False #Turn on or off the local error test on the sensitivity variables
        self.options["maxcorS"] = 3          #Maximum number of nonlinear iteration for sensitivity variables
        self.options["dqtype"] = "CENTERED"
        self.options["dqrhomax"] = 0.0
        self.options["pbar"] = [1]*self.problem_info["dimSens"]
        
        self.options["maxkrylov"] = 5
        self.options["precond"] = PREC_NONE
        
        #Statistics
        self.statistics["nfevals"]    = 0 #Function evaluations
        self.statistics["nsteps"]     = 0 #Number of steps
        self.statistics["netfails"]   = 0 #Number of error test failures
        self.statistics["nlinsetups"] = 0
        self.statistics["nncfails"]   = 0 #Nonlinear fails
        self.statistics["nniters"]    = 0 #Nonlinear iterations
        self.statistics["ngevals"]    = 0 #Root evaluations
        self.statistics["njevals"]    = 0 #Jacobian evaluations
        self.statistics["nfevalsLS"]  = 0 #Function evaluations due to Jac
        self.statistics["njvevals"]   = 0 #Number of Jacobian*Vector eval
        self.statistics["nfSevals"]   = 0 #Number of sensitivity evaluations
        self.statistics["nfevalsS"]   = 0 #Number of function evaluations due to finite sensitivity approximation
        self.statistics["nSetfails"]  = 0 #Number of error test failures
        self.statistics["nSniters"]   = 0 #Number of sensitivity nonlinear iterations
        self.statistics["nSncfails"]  = 0 #Number of sensitivity convergence failures
        
        #Solver support
        self.supports["one_step_mode"] = True
        self.supports["interpolated_output"] = True
        self.supports["interpolated_sensitivity_output"] = True
        self.supports["state_events"] = True
        
        #Get options from Problem
        if hasattr(problem, 'pbar'):
            self.pbar = problem.pbar
        elif hasattr(problem, 'p0'):
            self.pbar = N.array([N.abs(x) if N.abs(x) > 0 else 1.0 for x in self.problem.p0])
        if hasattr(problem, 'yS0'):
            self.yS0 = problem.yS0
    
    def __dealloc__(self):
        
        if self.yTemp != NULL:
            #Deallocate N_Vector
            N_VDestroy_Serial(self.yTemp)
        
        if self.cvode_mem != NULL:
            #Free Memory
            Sun.CVodeFree(&self.cvode_mem)
    
    cdef set_problem_data(self):
        
        #Sets the residual or rhs
        self.pt_fcn = self.problem.rhs
        self.pData.RHS = <void*>self.pt_fcn#<void*>self.problem.f
        self.pData.dim = self.problem_info["dim"]
        self.pData.memSize = self.pData.dim*sizeof(realtype)
        
        #Set the ndarray to the problem struct
        #self.yTemp   = N.zeros(self.pData.dim, dtype=N.float, order='c')
        #self.ydTemp  = N.zeros(self.pData.dim, dtype=N.float, order='c')
        #self.pData.y  = <void*>self.yTemp
        #self.pData.yd = <void*>self.ydTemp
        
        if self.problem_info["state_events"] is True: #Sets the root function
            self.pt_root = self.problem.state_events
            self.pData.ROOT = <void*>self.pt_root#<void*>self.problem.state_events
            self.pData.dimRoot = self.problem_info["dimRoot"]
            self.pData.memSizeRoot = self.pData.dimRoot*sizeof(realtype)

        if self.problem_info["jac_fcn"] is True: #Sets the jacobian
            self.pt_jac = self.problem.jac
            self.pData.JAC = <void*>self.pt_jac#<void*>self.problem.jac
            self.pData.memSizeJac = self.pData.dim*self.pData.dim*sizeof(realtype)
        
        if self.problem_info["jacv_fcn"] is True: #Sets the jacobian times vector
            self.pt_jacv = self.problem.jacv
            self.pData.JACV = <void*>self.pt_jacv#<void*>self.problem.jacv  
            
        if self.problem_info["sens_fcn"] is True: #Sets the sensitivity function
            self.pt_sens = self.problem.sens
            self.pData.SENS = <void*>self.pt_sens#<void*>self.problem.sens
           
        if self.problem_info["dimSens"] > 0: #Sensitivity parameters (does not need the sensitivity function)
            self.pData.dimSens = self.problem_info["dimSens"]    
            self.pData.p = <realtype*> malloc(self.problem_info["dimSens"]*sizeof(realtype))
            self.pData.pbar = <realtype*> malloc(self.problem_info["dimSens"]*sizeof(realtype))
        else:
            self.pData.dimSens = 0
            
        self.pData.verbose = 2
    
    cdef initialize_cvode(self):
        cdef int flag #Used for return
        cdef realtype ZERO = 0.0
        
        self.yTemp = arr2nv(self.y)
        
        if self.pData.dimSens > 0:
            #Create the initial matrices
            self.ySO  = N_VCloneVectorArray_Serial(self.pData.dimSens, self.yTemp)
            
            #Filling the start vectors
            for i in range(self.pData.dimSens):
                 N_VConst_Serial(ZERO,  self.ySO[i]);
                 if self.yS0 != None:
                    for j in range(self.pData.dim):
                        (<N_VectorContent_Serial>self.ySO[i].content).data[j] = self.yS0[i,j]


        #Updates the switches
        if self.problem_info["switches"]:
            self.pData.sw = <void*>self.sw
            
        if self.cvode_mem == NULL: #The solver is not initialized
            
            #Create the solver
            self.cvode_mem= Sun.CVodeCreate(CV_BDF if self.options["discr"] == "BDF" else CV_ADAMS, CV_NEWTON if self.options["iter"] == "Newton" else CV_FUNCTIONAL)
            if self.cvode_mem == NULL:
                raise CVodeError(CV_MEM_FAIL)
            
            #Specify the residual and the initial conditions to the solver
            flag = Sun.CVodeInit(self.cvode_mem, cv_rhs, self.t, self.yTemp)
            if flag < 0:
                raise CVodeError(flag, self.t)
                
            #Specify the root function to the solver
            if self.problem_info["state_events"]:
                flag = Sun.CVodeRootInit(self.cvode_mem, self.pData.dimRoot, cv_root)
                if flag < 0:
                    raise CVodeError(flag, self.t)
                    
            #Specify the error handling
            flag = Sun.CVodeSetErrHandlerFn(self.cvode_mem, cv_err, <void*>self.pData)
            if flag < 0:
                raise CVodeError(flag, self.t)
                
            #Sensitivity
            if self.pData.dimSens > 0:
                flag = Sun.CVodeSensInit(self.cvode_mem, self.pData.dimSens, CV_STAGGERED if self.options["sensmethod"] == "STAGGERED" else CV_SIMULTANEOUS, NULL, self.ySO)
                if flag < 0:
                    raise CVodeError(flag, self.t)
            
        else: #The solver needs to be reinitialized
            #Reinitialize
            flag = Sun.CVodeReInit(self.cvode_mem, self.t, self.yTemp)
            if flag < 0:
                raise CVodeError(flag, self.t)
            
            #Sensitivity
            if self.pData.dimSens > 0:
                flag = Sun.CVodeSensReInit(self.cvode_mem, CV_STAGGERED if self.options["sensmethod"] == "STAGGERED" else CV_SIMULTANEOUS, self.ySO)
                if flag < 0:
                    raise CVodeError(flag, self.t)
            
            
        #Set the user data
        flag = Sun.CVodeSetUserData(self.cvode_mem, <void*>self.pData)
        if flag < 0:
            raise CVodeError(flag, self.t)
            
    
    cpdef N.ndarray interpolate(self,double t,int k = 0):
        """
        Calls the internal CVodeGetDky for the interpolated values at time t.
        t must be within the last internal step. k is the derivative of y which
        can be from zero to the current order.
        """
        cdef flag
        cdef N.ndarray res
        cdef N_Vector dky=N_VNew_Serial(self.pData.dim) #Allocates a new N_Vector
        
        flag = Sun.CVodeGetDky(self.cvode_mem, t, k, dky)
        
        if flag < 0:
            raise CVodeError(flag, t)
        
        res = nv2arr(dky)
        
        #Deallocate N_Vector
        N_VDestroy_Serial(dky)
        
        return res
        
    cpdef N.ndarray interpolate_sensitivity(self, realtype t, int k = 0, int i=-1):
        """
        This method calls the internal method CVodeGetSensDky which computes the k-th derivatives
        of the interpolating polynomials for the sensitivity variables at time t.
        
            Parameters::
                    
                    t
                        - Specifies the time at which sensitivity information is requested. The time
                          t must fall within the interval defined by the last successful step taken
                          by CVodeS.
                    
                    k   
                        - The order of derivatives.
                        
                    i
                        - Specifies the sensitivity derivative vector to be returned (0<=i<=Ns)
                        
            Return::
            
                    A matrix containing the Ns vectors or a vector if i is specified.
        """
        cdef N_Vector dkyS=N_VNew_Serial(self.pData.dimSens)
        cdef int flag
        cdef N.ndarray res
        
        if i==-1:
            
            matrix = []
            
            for x in range(self.pData.dimSens):
                flag = Sun.CVodeGetSensDky1(self.cvode_mem, t, k, x, dkyS)
                if flag<0:
                    raise CVodeError(flag, t)
                
                matrix += [nv2arr(dkyS)]
            
            N_VDestroy_Serial(dkyS)
            
            return N.array(matrix)
        else:
            flag = Sun.CVodeGetSensDky1(self.cvode_mem, t, k, i, dkyS)
            if flag <0:
                raise CVodeError(flag, t)
            
            res = nv2arr(dkyS)
            
            N_VDestroy_Serial(dkyS)
            
            return res
    
    cpdef initialize(self):
        
        #Initialize storing of sensitivyt result in handle_result
        if self.problem_info['step_events'] or self.options['continuous_output']:
            self.problem._sensitivity_result = 1
        
        #Reset statistics
        for k in self.statistics.keys():
            self.statistics[k] = 0
        
        self.initialize_cvode()
    
    cpdef step(self,double t,N.ndarray y,double tf,dict opts):
        cdef int flag
        cdef N_Vector yout
        cdef double tret = 0.0
        cdef double tr
        cdef N.ndarray yr
        
        yout = arr2nv(y)
        
        #Get options
        initialize  = opts["initialize"]
        output_list = opts["output_list"]        
        
        #Initialize?
        if initialize:
            self.initialize_cvode()
            self.initialize_options()
        
        #Set stop time
        flag = Sun.CVodeSetStopTime(self.cvode_mem, tf)
        if flag < 0:
            raise CVodeError(flag, t)
        
        #Integration loop
        flag = Sun.CVode(self.cvode_mem,tf,yout,&tret,CV_ONE_STEP)
        if flag < 0:
            raise CVodeError(flag, tret)
            
        #Store results
        tr = tret
        yr = nv2arr(yout)
        
        if flag == CV_ROOT_RETURN: #Found a root
            flag = ID_EVENT #Convert to Assimulo flags
            self.store_statistics()
            
        if flag == CV_TSTOP_RETURN: #Reached tf
            flag = ID_COMPLETE
            self.store_statistics()
        
        #Deallocate
        N_VDestroy_Serial(yout)
                
        return flag, tr, yr
    
    cpdef integrate(self,double t,N.ndarray[ndim=1, dtype=realtype] y,double tf,dict opts):
        cdef int flag, output_index, normal_mode
        cdef N_Vector yout
        cdef double tret = 0.0, tout
        cdef list tr = [], yr = []
        cdef N.ndarray output_list
        
        yout = arr2nv(y)
        
        #Get options
        initialize   = opts["initialize"]
        
        #Initialize? 
        if initialize:
            self.initialize_cvode() 
            self.initialize_options()
        
        #Set stop time
        flag = Sun.CVodeSetStopTime(self.cvode_mem, tf)
        if flag < 0:
            raise CVodeError(flag, t)
        
        #Run in normal mode?
        normal_mode = 1 if opts["output_list"] != None else 0
        
        if normal_mode == 0: 
            #Integration loop
            while True:
                    
                flag = Sun.CVode(self.cvode_mem,tf,yout,&tret,CV_ONE_STEP)
                if flag < 0:
                    raise CVodeError(flag, tret)
                
                #Store results
                tr.append(tret)
                yr.append(nv2arr(yout))
                
                if flag == CV_ROOT_RETURN: #Found a root
                    flag = ID_EVENT #Convert to Assimulo flags
                    self.store_statistics()
                    break
                if flag == CV_TSTOP_RETURN: #Reached tf
                    flag = ID_COMPLETE
                    self.store_statistics()
                    break
        else:
            output_index = opts["output_index"]
            output_list  = opts["output_list"][output_index:]
            
            for tout in output_list:
                output_index += 1
                flag = Sun.CVode(self.cvode_mem,tout,yout,&tret,CV_NORMAL)
                if flag < 0:
                    raise CVodeError(flag, tret)
                
                #Store results
                tr.append(tret)
                yr.append(nv2arr(yout))
                
                if flag == CV_ROOT_RETURN: #Found a root
                    flag = ID_EVENT #Convert to Assimulo flags
                    self.store_statistics()
                    break
                if flag == CV_TSTOP_RETURN: #Reached tf
                    flag = ID_COMPLETE
                    self.store_statistics()
                    break
            else:
                flag = ID_COMPLETE
                self.store_statistics()
            
        
            opts["output_index"] = output_index
        
        #Deallocate
        N_VDestroy_Serial(yout)
        
        return flag, tr, yr
    
    cpdef state_event_info(self):
        """
        Returns the event info.
        """
        cdef int* c_info
        cdef flag
        
        # Allocate memory for the event_info_ vector and initialize to zeros
        c_info = <int*> malloc(self.pData.dimRoot*sizeof(int))
        for k in range(self.pData.dimRoot):
            c_info[k] = 0
        
        # Fetch data on which root functions that became zero and store in class
        flag = Sun.CVodeGetRootInfo(self.cvode_mem, c_info)
        if flag < 0:
            raise CVodeError(flag)
        
        #event_info = PyArray_SimpleNew(1,&self.pData.dimRoot,NPY_INT)
        event_info = [0]*self.pData.dimRoot
        
        for k in range(self.pData.dimRoot):
            event_info[k] = c_info[k]
        
        # Remember to deallocate
        free(c_info)
        
        return event_info
    
    cpdef initialize_sensitivity_options(self):
        cdef int flag
        
        #Sets the parameters to the userdata object.
        for i in range(self.pData.dimSens):
            self.pData.p[i] = self.p[i]
        #Sets the pbar to the userdata object.
        for i in range(self.pData.dimSens):
            self.pData.pbar[i] = self.options["pbar"][i]
        
        #Problem parameter information
        flag = Sun.CVodeSetSensParams(self.cvode_mem, self.pData.p, self.pData.pbar, NULL)
        if flag < 0:
            raise CVodeError(flag, self.t)
        
        #Difference quotient strategy
        flag = Sun.CVodeSetSensDQMethod(self.cvode_mem, CV_CENTERED if self.options["dqtype"]=="CENTERED" else CV_FORWARD, self.options["dqrhomax"])
        if flag < 0:
            raise CVodeError(flag, self.t)
        
        #Maximum number of nonlinear iterations
        flag = Sun.CVodeSetSensMaxNonlinIters(self.cvode_mem, self.options["maxcorS"])
        if flag < 0:
            raise CVodeError(flag, self.t)
        
        #Specify the error control strategy
        flag = Sun.CVodeSetSensErrCon(self.cvode_mem, self.options["suppress_sens"]==False)
        if flag < 0:
            raise CVodeError(flag, self.t)
        
        #Estimate the sensitivity
        flag = Sun.CVodeSensEEtolerances(self.cvode_mem)
        if flag < 0:
            raise CVodeError(flag, self.t)
        
        #Should the sensitivities be calculated this time around?
        if self.options["usesens"] == False:
            flag = Sun.CVodeSensToggleOff(self.cvode_mem)
            if flag < 0:
                raise CVodeError(flag, self.t)
    
    cpdef initialize_options(self):
        """
        Updates the simulation options.
        """
        cdef flag
        
        #Choose a linear solver if and only if NEWTON is choosen
        if self.options["linear_solver"] == 'DENSE' and self.options["iter"] == "Newton":
            #Specify the use of the internal dense linear algebra functions.
            flag = Sun.CVDense(self.cvode_mem, self.pData.dim)
            if flag < 0:
                raise CVodeError(flag)
                
            #Specify the jacobian to the solver
            if self.pData.JAC != NULL and self.options["usejac"]:
                flag = Sun.CVDlsSetDenseJacFn(self.cvode_mem, cv_jac)
                if flag < 0:
                    raise CVodeError(flag)
            else:
                flag = Sun.CVDlsSetDenseJacFn(self.cvode_mem, NULL)
                if flag < 0:
                    raise CVodeError(flag)
                    
        elif self.options["linear_solver"] == 'SPGMR' and self.options["iter"] == "Newton":
            #Specify the use of CVSPGMR linear solver.
            flag = Sun.CVSpgmr(self.cvode_mem, self.options["precond"], self.options["maxkrylov"])
            if flag < 0:
                raise CVodeError(flag)
                
            #Specify the jacobian times vector function
            if self.pData.JACV != NULL and self.options["usejac"]:
                flag = Sun.CVSpilsSetJacTimesVecFn(self.cvode_mem, cv_jacv)
                if flag < 0:
                    raise CVodeError(flag)
            else:
                flag = Sun.CVSpilsSetJacTimesVecFn(self.cvode_mem, NULL)
                if flag < 0:
                    raise CVodeError(flag)
        else: #Functional Iteration choosen.
            pass #raise CVodeError(100,t0) #Unknown error message

        #Maximum order
        flag = Sun.CVodeSetMaxOrd(self.cvode_mem, int(self.options["maxord"]))
        if flag < 0:
            raise CVodeError(flag)
            
        #Initial step
        flag = Sun.CVodeSetInitStep(self.cvode_mem, self.options["inith"])
        if flag < 0:
            raise CVodeError(flag)
        
        #Maximum step
        flag = Sun.CVodeSetMaxStep(self.cvode_mem, self.options["maxh"])
        if flag < 0:
            raise CVodeError(flag)
            
        #Maximum Number of steps
        flag = Sun.CVodeSetMaxNumSteps(self.cvode_mem, self.options["maxsteps"])
        if flag < 0:
            raise CVodeError(flag)
        
        #Tolerances
        flag = Sun.CVodeSVtolerances(self.cvode_mem, self.options["rtol"], arr2nv(self.options["atol"]))
        if flag < 0:
            raise CVodeError(flag)
            
        #Initialize sensitivity if any
        if self.pData.dimSens > 0:
            self.initialize_sensitivity_options()
    
    def _set_discr_method(self,discr='Adams'):
        
        if discr.upper() =='BDF':
            self.options["discr"] = "BDF"
            self.options["maxord"] = 5
                
        elif discr.upper() =='ADAMS':

            self.options["discr"] = "Adams"
            self.options["maxord"] = 12
        else:
            raise Exception('Discretization method must be either Adams or BDF')
            
        #Free Memory as we need another CVode memory object
        Sun.CVodeFree(&self.cvode_mem)
            
    def _get_discr_method(self):
        """
        This determines the discretization method.
        
            Parameters::
            
                discr   
                        - Default 'Adams', which indicates the use
                          of the Adams method. Can also be set to
                          'BDF' which indicates the use of the BDF
                          method.
                
                    Example:
                        discr = 'BDF'
        
        Note:: 
        
            Automatically sets the maximum order to 5 in the BDF case
            and to 12 in the Adams case. If necessary, change the maximum
            order After setting the discretization method.
        
        See SUNDIALS CVODE documentation 2.1 for more details.
        """
        return self.options["discr"]

    discr= property(_get_discr_method,_set_discr_method)
    
    def _set_iter_method(self,iter='FixedPoint'):

        if iter.upper()=='NEWTON':
            self.options["iter"] = "Newton"
        elif iter.upper()=='FIXEDPOINT':
            self.options["iter"] = "FixedPoint"
        else:
            raise Exception('Iteration method must be either FixedPoint or Newton')
            
        #Free Memory as we need another CVode memory object
        Sun.CVodeFree(&self.cvode_mem)
    
    def _get_iter_method(self):
        """
        This determines the iteration method that is be used by the
        solver.
        
            Parameters::
            
                iter    
                        - Default 'FixedPoint', which indicates the
                          use of a fixedpoint iteration method. Can
                          also be set to 'Newton' which indicates
                          the use of a Newton method.
                          
                            Example:
                                iter = 'Newton'
        
        See SUNDIALS CVODE documentation 2.1 for more details.
        """
        return self.options["iter"]
        
    iter = property(_get_iter_method,_set_iter_method)   
    
    def _set_atol(self,atol):
        
        self.options["atol"] = N.array(atol,dtype=N.float) if len(N.array(atol,dtype=N.float).shape)>0 else N.array([atol],dtype=N.float)
    
        if len(self.options["atol"]) == 1:
            self.options["atol"] = self.options["atol"]*N.ones(self.pData.dim)
        elif len(self.options["atol"]) != self.pData.dim:
            raise Exception("atol must be of length one or same as the dimension of the problem.")
        if (self.options["atol"]<0.0).any():
            raise Exception("atol must be positive.")
    
    def _get_atol(self):
        """
        Defines the absolute tolerance(s) that is to be used by the solver.
        Can be set differently for each variable.
        
            Parameters::
            
                atol    
                        - Default '1.0e-6'.
                
                        - Should be a positive float or a numpy vector
                          of floats.
                        
                            Example:
                                atol = [1.0e-4, 1.0e-6]
        
        See SUNDIALS IDA documentation 4.5.2 for more details.
        """
        return self.options["atol"]
    
    atol = property(_get_atol,_set_atol)
    
    def _set_rtol(self,rtol):
        
        try:
            rtol = float(rtol)
        except (ValueError, TypeError):
            raise Exception('Relative tolerance must be a (scalar) float.')
        if rtol <= 0.0:
            raise Exception('Relative tolerance must be a positive (scalar) float.')
        
        self.options["rtol"] = rtol
    
    def _get_rtol(self):
        """
        Defines the relative tolerance that is to be used by the solver.
        
            Parameters::
            
                rtol    
                        - Default '1.0e-6'.
                
                        - Should be a positive float.
                        
                            Example:
                                rtol = 1.0e-4
                                
        """
        return self.options["rtol"]
        
    rtol=property(_get_rtol,_set_rtol)
    
    def _set_max_ord(self,maxord):
        try:
            maxord = int(maxord)
        except ValueError:
            raise Exception("The maximal order must be an integer.")
        
        if self.options["discr"] == "Adams":
            if maxord > 12 or maxord < 1:
                self.options["maxord"] = 12 if maxord > 12 else 1
                self.log_message("The maximal order must be between 1-12.",NORMAL)
            else:
                self.options["maxord"] = maxord
        else:
            if maxord > 5 or maxord < 1:
                self.options["maxord"] = 5 if maxord > 5 else 1
                self.log_message("The maximal order must be between 1-5.", NORMAL)
            else:
                self.options["maxord"] = maxord
    
    def _get_max_ord(self):
        """
        This determines the maximal order that is be used by the solver.
        
            Parameters::
            
                maxord  
                        - Default '12', which is the maximum for the
                          Adams method, which is also default. For the
                          BDF method the maximum order is 5. 'maxord'
                          can be set in an interval from 1 to the
                          maximum order allowed.
                
                        - Should be an integer.
                        
                            Example:
                                maxord = 3
    
        
        An input value greater than the maximal order will result in the 
        maximum value.
        """
        return self.options["maxord"]

    maxord=property(_get_max_ord,_set_max_ord)
    
    def _set_linear_solver(self, lsolver):
        if lsolver.upper() == "DENSE" or lsolver.upper() == "SPGMR":
            self.options["linear_solver"] = lsolver.upper()
        else:
            raise Exception('The linear solver must be either "DENSE" or "SPGMR".')
        
    def _get_linear_solver(self):
        """
        Specifies the linear solver to be used.
        
            Parameters::
            
                linearsolver
                        - Default 'DENSE'. Can also be 'SPGMR'.
        """
        return self.options["linear_solver"]
    
    linear_solver = property(_get_linear_solver, _set_linear_solver)
    
    def _set_initial_step(self, initstep):
        try:
            self.options["inith"] = float(initstep)
        except (ValueError, TypeError):
            raise Exception('The initial step must be an integer or float.')
        
    def _get_initial_step(self):
        """
        This determines the initial step-size to be used in the integration.
        
            Parameters::
            
                initstep    
                            - Default '0.0', which result in that the
                              the initial step is approximated.
                            
                            - Should be float.
                            
                                Example:
                                    initstep = 0.01
        """
        return self.options["inith"]
        
    inith = property(_get_initial_step,_set_initial_step)
    
    def _set_max_steps(self, maxsteps):
        if not isinstance(maxsteps,int):
            raise Exception('The maximum number of steps must be an integer.')
        if maxsteps < 1:
            raise Exception('The maximum number of steps must be a positive integer.')
        self.options["maxsteps"] = maxsteps
    
    def _get_max_steps(self):
        """
        Determines the maximum number of steps the solver is allowed
        to take to finish the simulation.
        
            Parameters::
            
                maxsteps    
                            - Default '10000'.
                            
                            - Should be an integer.
                            
                                Example:
                                    maxsteps = 1000.
                                    
        """
        return self.options["maxsteps"]
    
    def _set_max_cor_S(self,maxcorS):
        try:
            self.options["maxcorS"] = int(maxcorS)
        except:
            raise Exception("The maximum nonlinear sensitivity iterations must be a positiv integer.")
    
    def _get_max_cor_S(self):
        """
        This detmines the maximum number of nonlinear iterations for the
        sensitivity variables.
        
            Parameters::
            
                maxcorS
                        - Default 3
                        
                        - Should be an integer
                        
        For more information see SUNDIALS IDAS documentation 5.2.6.
        """
        return self.options["maxcorS"]
    
    maxcorS=property(_get_max_cor_S,_set_max_cor_S)
    
    maxsteps = property(_get_max_steps, _set_max_steps)
    
    def _set_max_h(self,max_h):
        try:
            self.options["maxh"] = float(max_h)
        except:
            raise Exception("Maximal stepsize must be a (scalar) float.")
    
    def _get_max_h(self):
        """
        Defines the maximal step-size that is to be used by the solver.
        
            Parameters::
            
                maxh    
                        - Default '0', which indicates that the maximal
                          step-size is infinity.
                          
                        - Should be a float.
                        
                            Example:
                                maxh = 0.01
                                
        """
        return self.options["maxh"]
        
    maxh=property(_get_max_h,_set_max_h)
    
    def _set_usejac(self, jac):
        self.options["usejac"] = bool(jac)
    
    def _get_usejac(self):
        """
        This sets the option to use the user defined jacobian. If a
        user provided jacobian is implemented into the problem the
        default setting is to use that jacobian. If not, an
        approximation is used.
        
            Parameters::
            
                usejac  
                        - True - use user defined jacobian
                          False - use an approximation
                    
                        - Should be a boolean.
                        
                            Example:
                                usejac = False
        """
        return self.options["usejac"]
    
    usejac = property(_get_usejac,_set_usejac)
    
    def _set_usesens(self, sens):
        self.options["usesens"] = bool(sens)
    
    def _get_usesens(self):
        """
        This options activates or deactivates the sensitivity 
        calculations.
        
            Parameters::
            
                usejac  
                        - True -  Activate sensitivity calculations
                          False - Deactivate sensitivity calculations
                    
                        - Should be a boolean.
                        
                            Example:
                                usesens = False
        """
        return self.options["usesens"]
    
    usesens = property(_get_usesens,_set_usesens)
    
    def _set_sensitivity_method(self, ism):
        if not isinstance(ism, str):
            raise Exception('Sensitivity method must be string.')
        
        if ism.upper() == 'SIMULTANEOUS':
            self.options["sensmethod"] = 'SIMULTANEOUS' 
        elif ism.upper() == 'STAGGERED':
            self.options["sensmethod"] = 'STAGGERED'
        else:
            raise Exception('Sensitivity method must be either "SIMULTANEOUS" or "STAGGERED".')
        
    def _get_sensitivity_method(self):
        """
        Specifies the sensitivity solution method. Can be either,
        'SIMULTANEOUS' or 'STAGGERED'.
            
        Parameters::
            
            ism
                    - A string of either 'SIMULTANEOUS' or 'STAGGERED'
                    - Default 'STAGGERED'
                        
        Returns::
            
            The current value of sensmethod (string)
        
        See SUNDIALS documentation '(IDA/CVode)SensInit'
        """
        return self.options["sensmethod"]
        
    sensmethod = property(_get_sensitivity_method, _set_sensitivity_method)
    
    def _set_suppress_sens(self,suppress_sens):
        try:
            self.options["suppress_sens"] = bool(suppress_sens)
        except:
            raise Exception("Unkown input to suppress_sens, must be a boolean.")

    def _get_suppress_sens(self):
        """
        A Boolean flag which indicates that the error-tests are 
        suppressed on the sensitivity variables.
        
            Parameters::
            
                suppress_sens    
                                - Default 'False'.
                
                                - Should be a boolean.
                                
                                    Example:
                                        suppress_alg = True
                                        
        See SUNDIALS CVode documentation 5.2.6 'CVodeSetSensErrCon' 
        for more details.
        """
        return self.options["suppress_sens"]    

    suppress_sens=property(_get_suppress_sens,_set_suppress_sens)
    
    def _set_dqtype(self, dqtype):
        if not isinstance(dqtype, str):
            raise Exception('DQtype must be string.')
        
        if dqtype.upper() == 'CENTERED':
            self.options["dqtype"] = "CENTERED"
        elif dqtype.upper() == 'FORWARD':
            self.options["dqtype"] = "FORWARD"
        else:
            raise Exception('DQtype must be either "CENTERED" or "FORWARD".')
            
    def _get_dqtype(self):
        """
        Specifies the difference quotient type in the sensitivity calculations.
        Can be either, 'CENTERED' or 'FORWARD'.
        
        Parameters::
            
            DQtype 
                    - A string of either 'CENTERED' or 'FORWARD'
                    - Default 'CENTERED'
                        
        Returns::
            
            The current value of DQtype.
        
        See SUNDIALS documentation '(IDA/CVode)SetSensDQMethod' 
        """
        return self.options["dqtype"]
    
    dqtype = property(_get_dqtype, _set_dqtype)
    
    def _set_dqrhomax(self, dqrhomax):
        try:
            self.options["dqrhomax"] = float(dqrhomax)
        except (TypeError, ValueError):
            raise Exception('DQrhomax must be convertable to a float.')
        if self.options["dqrhomax"] < 0.0:
            raise Exception("DQrhomax must be positive.")
            
    def _get_dqrhomax(self):
        """
        Specifies the selection parameters used in deciding switching between a simultaneous
        or separate approximation of the two terms in the sensitivity residual.
        
            Parameters::
            
                DQrhomax
                        - A postive float.
                        - Default 0.0
                        
            Returns::
            
                The current value of DQrhomax (float)
        
        See SUNDIALS documentation '(IDA/CVode)SetSensDQMethod' 
        """
        return self.options["dqrhomax"]
    
    dqrhomax = property(_get_dqrhomax, _set_dqrhomax)
    
    def _set_max_krylov(self, maxkrylov):
        try:
            self.options["maxkrylov"] = int(maxkrylov)
        except:
            raise Exception("Maximum number of krylov dimension should be an integer.")
        if self.options["maxkrylov"] < 0:
            raise Exception("Maximum number of krylov dimension should be an positive integer.")
            
    def _get_max_krylov(self):
        """
        Specifies the maximum number of dimensions for the krylov subspace to be used.
        
            Parameters::
            
                    maxkrylov
                            - A positive integer.
                            - Default 0
            
            Returns::
            
                The current value of maxkrylov.
                
        See SUNDIALS documentation 'CVSpgmr'
        """
        return self.options["maxkrylov"]
    
    maxkrylov = property(_get_max_krylov, _set_max_krylov)
    
    def _set_pre_cond(self, precond):
        if precond.upper() == "PREC_NONE":
            self.options["precond"] = PREC_NONE
        elif precond.upper() == "PREC_LEFT":
            self.options["precond"] = PREC_LEFT
        elif precond.upper() == "PREC_RIGHT":
            self.options["precond"] = PREC_RIGHT
        elif precond.upper() == "PREC_BOTH":
            self.options["precond"] = PREC_BOTH
        else:
            raise Exception('Unknown input of precond. Should be either "PREC_NONE", "PREC_LEFT","PREC_RIGHT" or "PREC_BOTH"')
    def _get_pre_cond(self):
        """
        Specifies the preconditioning type.
        
            Parameters::
            
                    precond
                            - Should be either "PREC_NONE", "PREC_LEFT"
                              "PREC_RIGHT" or "PREC_BOTH"
                            - Default PREC_NONE
            
            Returns::
            
                The current value of precond (as string).
                
        See SUNDIALS documentation 'CVSpgmr'
        """
        if self.options["precond"] == PREC_NONE:
            return "PREC_NONE"
        elif self.options["precond"] == PREC_LEFT:
            return "PREC_LEFT"
        elif self.options["precond"] == PREC_RIGHT:
            return "PREC_RIGHT"
        elif self.options["precond"] == PREC_BOTH:
            return "PREC_BOTH"
    
    precond = property(_get_pre_cond, _set_pre_cond)
    
    def _set_pbar(self, pbar):
        if len(pbar) != self.problem_info['dimSens']:
            raise Exception('pbar must be of equal length as the parameters.')
        
        self.options["pbar"] = pbar
    
    def _get_pbar(self):
        """
        Specifies the order of magnitude for the parameters. This is useful if IDAS is
        to estimate tolerances for the sensitivity solution vectors.
        
            Parameters::
            
                pbar
                        - An array of positive floats equal to the number of parameters.
                        - Default absolute values of the parameters.
                        
            Returns::
            
                The current value of pbar.
                
        See SUNDIALS documentation '(IDA/CVode)SetSensParams'
        """
        return self.options["pbar"]
    
    pbar = property(_get_pbar, _set_pbar)
    
    cdef void store_statistics(self):
        """
        Retrieves and stores the statistics.
        """
        cdef long int nsteps = 0, njevals = 0, ngevals = 0, netfails = 0, nniters = 0, nncfails = 0
        cdef long int nSniters = 0, nSncfails = 0, nfevalsLS = 0, njvevals = 0, nfevals = 0
        cdef long int nfSevals = 0,nfevalsS = 0,nSetfails = 0,nlinsetupsS = 0, nlinsetups = 0
        cdef int qlast = 0, qcur = 0
        cdef realtype hinused = 0.0, hlast = 0.0, hcur = 0.0, tcur = 0.0

        if self.options["linear_solver"] == "SPGMR":
            flag = Sun.CVSpilsGetNumJtimesEvals(self.cvode_mem, &njvevals) #Number of jac*vector
            flag = Sun.CVSpilsGetNumRhsEvals(self.cvode_mem, &nfevalsLS) #Number of rhs due to jac*vector
        else:
            flag = Sun.CVDlsGetNumJacEvals(self.cvode_mem, &njevals) #Number of jac evals
            flag = Sun.CVDlsGetNumRhsEvals(self.cvode_mem, &nfevalsLS) #Number of res evals due to jac evals
            
        flag = Sun.CVodeGetNumGEvals(self.cvode_mem, &ngevals) #Number of root evals
        
        #Get all integrator statistics
        flag = Sun.CVodeGetIntegratorStats(self.cvode_mem, &nsteps, &nfevals, &nlinsetups, &netfails, &qlast,
                                       &qcur, &hinused, &hlast, &hcur, &tcur)
        
        flag = Sun.CVodeGetNonlinSolvStats(self.cvode_mem, &nniters, &nncfails) #Number of nonlinear iteration
                                                                            #Number of nonlinear conv failures
        
        self.statistics["nsteps"]    += nsteps
        self.statistics["nfevals"]   += nfevals
        self.statistics["netfails"]  += netfails
        self.statistics["nniters"]   += nniters
        self.statistics["nncfails"]  += nncfails
        self.statistics["nfevalsLS"] += nfevalsLS
        self.statistics["njvevals"]  += njvevals
        self.statistics["njevals"]   += njevals
        self.statistics["ngevals"]   += ngevals
        
        #If sensitivity    
        if self.pData.dimSens > 0:
            flag = Sun.CVodeGetSensStats(self.cvode_mem, &nfSevals, &nfevalsS, &nSetfails, &nlinsetupsS)
            flag = Sun.CVodeGetSensNonlinSolvStats(self.cvode_mem, &nSniters, &nSncfails)
            
            self.statistics["nfSevals"]   += nfSevals
            self.statistics["nfevalsS"]   += nfevalsS
            self.statistics["nSetfails"]  += nSetfails
            self.statistics["nSniters"]   += nSniters
            self.statistics["nSncfails"]  += nSncfails
                
    def print_statistics(self, verbose=NORMAL):
        """
        Should print the statistics.
        """
        self.log_message('Final Run Statistics: %s \n' % self.problem.name,        verbose)
        
        self.log_message(' Number of Steps                          : '+str(self.statistics["nsteps"]),          verbose)               
        self.log_message(' Number of Function Evaluations           : '+str(self.statistics["nfevals"]),         verbose)
        if self.options["linear_solver"] == "SPGMR":
            self.log_message(' Number of Jacobian*Vector Evaluations    : ' + str(self.statistics["njvevals"]),  verbose)
            self.log_message(' Number of F-Evals During Jac*Vec-Evals   : ' + str(self.statistics["nfevalsLS"]), verbose)
        else:     
            self.log_message(' Number of Jacobian Evaluations           : '+ str(self.statistics["njevals"]),    verbose)
            self.log_message(' Number of F-Eval During Jac-Eval         : '+ str(self.statistics["nfevalsLS"]),  verbose)
        self.log_message(' Number of Root Evaluations               : '+ str(self.statistics["ngevals"]),        verbose)
        self.log_message(' Number of Error Test Failures            : '+ str(self.statistics["netfails"]),       verbose)
        if self.options["iter"] == "FixedPoint":
            self.log_message(' Number of Functional Iterations          : '+ str(self.statistics["nniters"]),        verbose)
            self.log_message(' Number of Functional Convergence Failures: '+ str(self.statistics["nncfails"]),       verbose)
        else:
            self.log_message(' Number of Newton Iterations              : '+ str(self.statistics["nniters"]),        verbose)
            self.log_message(' Number of Newton Convergence Failures    : '+ str(self.statistics["nncfails"]),       verbose)
        
        
        if self.problem_info['dimSens'] > 0: #Senstivity calculations is on
            self.log_message('\nSensitivity Statistics:\n', verbose)
            self.log_message(' Number of Sensitivity Calculations             : ' + str(self.statistics["nfSevals"]), verbose)
            self.log_message(' Number of F-Evals Due to Finite Approximation  : ' + str(self.statistics["nfevalsS"]), verbose)
            self.log_message(' Number of Local Error Test Failures            : ' + str(self.statistics["nSetfails"]),verbose)
            self.log_message(' Number of Newton Iterations                    : ' + str(self.statistics["nSniters"]), verbose)
            self.log_message(' Number of Newton Convergance Failures          : ' + str(self.statistics["nSncfails"]),verbose)
            
            self.log_message('\nSensitivity options:\n' , verbose)
            self.log_message(' Method                   : ' + str(self.options["sensmethod"]), verbose)
            self.log_message(' Difference quotient type : ' + str(self.options["dqtype"]), verbose)
            self.log_message(' Suppress Sens            : ' + str(self.options["suppress_sens"]), verbose)
    
        self.log_message('\nSolver options:\n',                                      verbose)
        self.log_message(' Solver                  : CVode',                         verbose)
        self.log_message(' Linear Multistep Method : ' +self.options["discr"],       verbose)
        self.log_message(' Nonlinear Solver        : ' + self.options["iter"],       verbose)
        self.log_message(' Maxord                  : ' + str(self.options["maxord"]),verbose)
        self.log_message(' Tolerances (absolute)   : ' + str(self.options["atol"]),  verbose)
        self.log_message(' Tolerances (relative)   : ' + str(self.options["rtol"]),  verbose)
        self.log_message('',                                                         verbose)


class CVodeError(Exception):
    """
    Defines the CVodeError and provides the textual error message.
    """
    msg = { CV_TOO_MUCH_WORK     : 'The solver took max internal steps but could not reach tout.',
            CV_TOO_MUCH_ACC      : 'The solver could not satisfy the accuracy demanded by the user for some internal step.',
            CV_ERR_FAIL          : 'Error test failures occurred too many times during one internal time step or minimum step size was reached.',
            CV_CONV_FAIL         : 'Convergence test failures occurred too many times during one internal time step or minimum step size was reached.',
            CV_LINIT_FAIL        : 'The linear solvers initialization function failed.',
            CV_LSETUP_FAIL       : 'The linear solvers setup function failed in an unrecoverable manner.',
            CV_LSOLVE_FAIL       : 'The linear solvers solve function failed in an unrecoverable manner.',
            CV_RHSFUNC_FAIL      : 'The user-provided rhs function failed in an unrecoverable manner.',
            CV_FIRST_RHSFUNC_ERR : 'The right-hand side function failed at the first call.',
            CV_REPTD_RHSFUNC_ERR : 'The right-hand side function had repetead recoverable errors.',
            CV_UNREC_RHSFUNC_ERR : 'The right-hand side function had a recoverable error, but no recovery is possible.',
            CV_RTFUNC_FAIL       : 'The rootfinding function failed in an unrecoverable manner.',
            CV_MEM_FAIL          : 'A memory allocation failed.',
            CV_MEM_NULL          : 'The cvode_mem argument was NULL.',
            CV_ILL_INPUT         : 'One of the function inputs is illegal.',
            CV_NO_MALLOC         : 'The CVode memory block was not allocated by a call to CVodeMalloc.',
            CV_BAD_K             : 'The derivative order k is larger than the order used.',
            CV_BAD_T             : 'The time t is outside the last step taken.',
            CV_BAD_DKY           : 'The output derivative vector is NULL.',
            CV_TOO_CLOSE         : 'The output and initial times are too close to each other.'}
    
    def __init__(self, value, t = 0.0):
        self.value = value
        self.t = t
        
    def __str__(self): 
        try:
            return repr(self.msg[self.value]+' At time %f.'%self.t)    
        except KeyError:
            return repr('Sundials failed with flag %s. At time %f.'%(self.value, self.t))   



class IDAError(Exception):
    """
    Defines the IDAError and provides the textual error message.
    """
    msg = { IDA_TOO_MUCH_WORK    : 'The solver took max internal steps but could not reach tout.',
            IDA_TOO_MUCH_ACC     : 'The solver could not satisfy the accuracy demanded by the user for some internal step.',
            IDA_ERR_FAIL         : 'Error test failures occurred too many times during one internal time step or minimum step size was reached.',
            IDA_CONV_FAIL        : 'Convergence test failures occurred too many times during one internal time step or minimum step size was reached.',
            IDA_LINIT_FAIL       : 'The linear solvers initialization function failed.',
            IDA_LSETUP_FAIL      : 'The linear solvers setup function failed in an unrecoverable manner.',
            IDA_LSOLVE_FAIL      : 'The linear solvers solve function failed in an unrecoverable manner.',
            IDA_RES_FAIL         : 'The user-provided residual function failed in an unrecoverable manner.',
            IDA_REP_RES_FAIL     : 'The user-provided residual function repeatedly returned a recoverable error flag, but the solver was unable to recover.',
            IDA_RTFUNC_FAIL      : 'The rootfinding function failed in an unrecoverable manner.',
            IDA_CONSTR_FAIL      : 'The inequality constraints were violated and the solver was unable to recover.',
            IDA_FIRST_RES_FAIL   : 'The user-provided residual function failed recoverable on the first call.',
            IDA_LINESEARCH_FAIL  : 'The line search failed.',
            IDA_NO_RECOVERY      : 'The residual function, linear solver setup function or linear solver solve function had a recoverable failure. But IDACalcIC could not recover.',
            IDA_MEM_NULL         : 'The ida_mem argument was NULL.',
            IDA_MEM_FAIL         : 'A memory allocation failed.',
            IDA_ILL_INPUT        : 'One of the function inputs is illegal.',
            IDA_NO_MALLOC        : 'The IDA memory was not allocated by a call to IDAInit.',
            IDA_BAD_EWT          : 'Zero value of some error weight component.',
            IDA_BAD_K            : 'The k-th derivative is not available.',
            IDA_BAD_T            : 'The time t is outside the last step taken.',
            IDA_BAD_DKY          : 'The vector argument where derivative should be stored is NULL.',
            IDA_SRES_FAIL        : 'The user-provided sensitivity residual function failed in an unrecoverable manner.',
            IDA_REP_SRES_ERR     : 'The user-provided sensitivity residual function repeatedly returned a recoverable error flag, but the solver was unable to recover.',
            IDA_BAD_IS           : 'The sensitivity identifier is not valid.'}
    
    def __init__(self, value, t = 0.0):
        self.value = value
        self.t = t
        
    def __str__(self): 
        try:
            return repr(self.msg[self.value]+' At time %f.'%self.t)    
        except KeyError:
            return repr('Sundials failed with flag %s. At time %f.'%(self.value, self.t))

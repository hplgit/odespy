#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2010 Modelon AB
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

"""
Cython Wrapper for interfacing Python with KINSOL (Sundials Version 2.4.0)
Johan Ylikiiskilä Modelon AB

see also Jon Olav Vik: 
http://codespeak.net/pipermail/cython-dev/2009-June/005947.html

"""

from __future__ import division

include "sundials_kinsol_includes.pxi" # Includes fcts from other header files
include "sundials_kinsol_core.pxi" # Includes the constants (textual include)

cdef int kin_res(N_Vector xv, N_Vector fval, void *problem_data):
    """
    Residual fct called by KINSOL
    """

    cdef:
        int i
        ProblemData pData = <ProblemData>problem_data
        ndarray[realtype, ndim = 1, mode = 'c'] rhs # Container holding return from user fct
        
        ndarray x = nv2arr(xv)
        realtype* resptr = (<N_VectorContent_Serial>fval.content).data

    try:
        rhs = (<object>pData.RHS)(x)
        for i in range(pData.dim):
            
            resptr[i] = rhs[i]

        return KIN_SUCCESS
    except:
        return KIN_REC_ERR

cdef int kin_jac(int Neq, N_Vector xv, N_Vector fval, DlsMat Jacobian, 
                void *problem_data, N_Vector tmp1, N_Vector tmp2):
    """
    This method is used to connect the Assimulo.Problem.jac to the Sundials
    Jacobian function.
    """
    cdef:
        ProblemData pData = <ProblemData>problem_data
        realtype* col_i=DENSE_COL(Jacobian,0)
        ndarray x = nv2arr(xv)
    
    try:
 
        jac=(<object>pData.JAC)(x)
        #This needs further investigations:
        #memcpy(Jacobian.data,<realtype*>jac.data, pData.memSizeJac)
        
        for i in range(Neq):
            col_i = DENSE_COL(Jacobian, i)
            for j in range(Neq):
                col_i[j] = jac[j,i]

        return KINDLS_SUCCESS
    except:
        return KINDLS_JACFUNC_RECVR #Recoverable Error (See Sundials description)


class KINError(Exception):
    """
    Kinsol exception
    """
    
    msg = {KIN_MEM_FAIL : 'Memory allocation failed.',
                KIN_MEM_NULL : 'KINSOL not properly created.',
                KIN_ILL_INPUT: 'There was an input error.',
                KIN_NO_MALLOC: 'Memory not allocated, call KINSOL_init(...)',
                KIN_LINESEARCH_NONCONV: 'Linesearch could not find a step big enough.',
                KIN_MAXITER_REACHED: 'Max number of iterations reached',
                KIN_MXNEWT_5X_EXCEEDED: 'Max size of newton step exceeded 5 times.',
                KIN_LINESEARCH_BCFAIL: 'Linesearch beta-test failed, probably because of too poor progress.',
                KIN_LINSOLV_NO_RECOVERY: 'psolve encountered an error, but preconditioner is current',
                KIN_LINIT_FAIL: 'Init of linear solver failed.',
                KIN_LSETUP_FAIL: 'pset (preconditioner setup fct) encountered an error',
                KIN_LSOLVE_FAIL: 'Error in either psolve or in linear solver routine.',
                KIN_SYSFUNC_FAIL: 'Call to RHS failed.',
                KIN_FIRST_SYSFUNC_ERR: 'Call to RHS failed on first call',
                KIN_REPTD_SYSFUNC_ERR: 'Call to RHS failed multiple times.',
                KIN_STEP_LT_STPTOL: 'Scaled step length too small. Either an approximate solution or a local minimum is reached. Check value of residual.',
                234: 'KINSOL not compiled with SuperLU, sparse functionality not available.'}
    value = 0
    def __init__(self, value):
        self.value = value

    def __str__(self):
        try:
            return repr(self.msg[self.value])
        except KeyError:
            return repr('KINSOL failed with flag %s.'%(self.value))

cdef class KINSOL_wrap:

    cdef:
        void* solver 
        ProblemData pData # A struct containing problem data
        #ProblemData *ppData # Pointer to above mentioned struct
        N_Vector x_cur,x_init, x_scale, f_scale, con_nv
        booleantype noInitSetup, sparse
        realtype norm_of_res,reg_param
        ndarray con
        int print_level
        ndarray fscale

    def __cinit__(self):
        self.pData = ProblemData() #Create a new problem struct
        #self.ppData = &self.pData

    cpdef KINSOL_set_problem_info(self, RHS, dim, JAC= None):
        """
        Sets the problem info to the desired struct
        """

        # Sets residual function
        self.pData.RHS = <void*>RHS
        self.pData.dim = dim
        if JAC != None:
            self.pData.JAC = <void*>JAC
        else:
            self.pData.JAC = NULL

    def KINSOL_init(self,RHS,x0,dim, JAC = None, con = None,sparse = False, print_level = 0,norm_of_res = 0,reg_param=0,fscale=None):
        """
        Initializes solver
        """        
        cdef int flag
        # Set problem info
        self.noInitSetup = False
        self.print_level = print_level
        self.norm_of_res = norm_of_res
        self.reg_param = reg_param
        self.KINSOL_set_problem_info(RHS,dim,JAC)
        self.con = con
        self.sparse = sparse
        self.fscale = fscale


        # Create initial guess from the supplied numpy array
        # print "x0 got in KINSOL wrapper: ", x0
        self.x_cur = arr2nv(x0)
        self.x_init = arr2nv(x0)


        if self.solver == NULL: # solver runs for the first time
            
            # Create solver
            self.solver = KINCreate()
            if self.solver == NULL:
                raise KINError(KIN_MEM_FAIL)
            if self.print_level >0:
                print "KINSOL created"
            
            # Stop solver from performing precond setup
            flag = KINSetNoInitSetup(self.solver, self.noInitSetup)
            if flag < 0:
                raise KINError(flag)

            flag = KINSetPrintLevel(self.solver, self.print_level)
            if flag < 0:
                raise KINError(flag)

            if self.norm_of_res >0:
                flag = KINSetMaxNewtonStep(self.solver, self.norm_of_res)
                if flag < 0:
                    raise KINError(flag)
            else:
                flag = KINSetMaxNewtonStep(self.solver, 1e20)
                if flag < 0:
                    raise KINError(flag)

            # If the user has specified constraints, connect them
            if self.con != None:
                if self.print_level >0:
                    print "Applying constraints"
                self.con_nv = arr2nv(self.con)
                flag = KINSetConstraints(self.solver, self.con_nv)
                if flag < 0:
                    raise KINError(flag)
            
            # Allocate internal memory
            flag = KINInit(self.solver,kin_res,self.x_cur)
            if flag < 0:
                raise KINError(flag)
            if self.print_level >0:
                print "KINSOL initialized"

            # Link to linear solver
            if self.sparse:
                raise KINError(234)
            else:
                flag = KINPinv(self.solver,self.pData.dim)
                #flag = KINDense(self.solver,self.pData.dim)
            
            if flag < 0:
                raise KINError(flag)
            if self.print_level >0:
                print "Linear solver connected"

            # Set regularization parameter, if necessary
            if self.reg_param != 0.0:
                if self.sparse:
                    raise KINError(234)
                else:
                    flag = KINPinvSetRegParam(self.solver, self.reg_param)
                    
                if flag < 0:
                    raise KINError(flag)
                if self.print_level >0:
                    print "Reg param ", self.reg_param, " set"
        else:
            # If the user has specified constraints, connect them
            if con != None:
                if self.print_level >0:
                    print "Applying constraints"
                self.con_nv = arr2nv(con)
                flag = KINSetConstraints(self.solver, self.con_nv)
                if flag < 0:
                    raise KINError(flag)
            else:
                flag = KINSetConstraints(self.solver, NULL)
                if flag < 0:
                    raise KINError(flag)
        
        # If the user supplied a Jacobien, link it to the solver
        if self.pData.JAC != NULL:
            if sparse:
                raise KINError(234)
            else:
                flag = KINPinvSetJacFn(self.solver,kin_jac)
                #flag = KINDlsSetDenseJacFn(self.solver,kin_jac)
                
            if flag < 0:
                raise KINError(flag)
            if self.print_level >0:
                print "Jacobian supplied by user connected"
        else:
            if sparse:
                raise KINError(234)
            else:
                flag = KINPinvSetJacFn(self.solver,NULL)
                #flag = KINDlsSetDenseJacFn(self.solver,NULL)

            if flag < 0:
                raise KINError(flag)

        # Link the solver to the supplied problem data
        flag = KINSetUserData(self.solver,<void*>self.pData)
        if flag < 0:
            raise KINError(flag)
        if self.print_level >0:
            print "User data set"

        # Create scaling vectors filled with ones
        # since the problem is assumed to be scaled
        self.x_scale = arr2nv(np.ones(self.pData.dim))
        if self.fscale != None:
            self.f_scale = arr2nv(self.fscale)
        else:
            self.f_scale = arr2nv(np.ones(self.pData.dim))

    
    def KINSOL_solve(self,wo_LineSearch):
        """
        Function that should be called after a call to KINSOL_init
        solves the function supplied as RHS
        """
        
        cdef:
            int flag, DLSflag, lsflag, count,ind
            long int nb_iters
            ndarray x0
        if self.print_level >0:
            print "Calling solver..."
        
        if not wo_LineSearch:
            # Use LineSearch
            flag = KINSol(<void*>self.solver,self.x_cur,KIN_LINESEARCH,self.x_scale,self.f_scale)
            count = 0
            if flag > 1:
                # Probably stuck at a minimum
                self.x_cur = self.x_init
                flag = KINSol(<void*>self.solver,self.x_cur,KIN_NONE,self.x_scale,self.f_scale)
            
            if flag > 1:
                # Stuck at minimum again
                raise KINError(flag)
        else:
            # Don't use linesearch
            flag = -7

        while flag <0:
            # Get nbr of non-linear iterations
            lsflag = KINGetNumNonlinSolvIters(<void*>self.solver, &nb_iters)
            if lsflag != 0:
                raise KINError(lsflag)
            
            
            if (flag == -5 and nb_iters == 1):
                # Cannot take a step bigger than zero on the first iteration
                # Try with heuristic initial guess
                raise KINError(42)

            if (flag == -7 or flag == -6 or (flag == -5 and nb_iters > 1)):
                if count == 0:
                    # First time
                    print "Calling solver again as simple newton."
                    self.x_cur = self.x_init
                    DLSflag = KINSetMaxSetupCalls(self.solver, 1);
                    if DLSflag < 0:
                        raise KINError(DLSflag)
                    
                elif count == 1:
                    # second time
                    print "Turning of residual monitoring"
                    DLSflag = KINSetNoResMon(self.solver, True);
                    if DLSflag < 0:
                        raise KINError(DLSflag)
                    
                elif count == 3:
                    # fourth time, release the bound on the step
                    
                    print "Remove limit on stepsize"
                    DLSflag = KINSetMaxNewtonStep(self.solver, 1.0e34)
                    if DLSflag < 0:
                        raise KINError(DLSflag)

                count += 1
                if count >= 10:
                    if self.print_level >0:
                        print "Tried with simple newton ten times"
                    raise KINError(-7)

            else:

                print "Error in solve, flag: ", flag
                lsflag = KINPinvGetLastFlag(self.solver, &lsflag)
                #lsflag = KINDlsGetLastFlag(self.solver, &lsflag)
                if lsflag != 0:
                    if lsflag <0:
                        print "Last flag from Linear solver: ", lsflag
                    else:
                        print "Jacobian singular at column ", lsflag
                raise KINError(flag)

            flag = KINSol(<void*>self.solver,self.x_cur,KIN_NONE,self.x_scale,self.f_scale)
        if self.print_level >0:
            print "Problem solved, returning result"

        return nv2arr(self.x_cur)

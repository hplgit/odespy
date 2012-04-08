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
Cython Wrapper for interfacing Python with CVode and IDA (Sundials Version 2.4.0)
Claus Fuhrer,        Lund University        
Christian Andersson, Lund University        

see also Jon Olav Vik: 
http://codespeak.net/pipermail/cython-dev/2009-June/005947.html

"""
import numpy as N
cimport numpy as N

from numpy cimport NPY_DOUBLE, npy_intp, NPY_INT

#==============================================
#External definitions from Sundials headers
#==============================================

cdef extern from "sundials/sundials_types.h":
    ctypedef double realtype
    ctypedef bint booleantype # should be bool instead of bint, but there is a bug in Cython

#==============================================
# C headers
#==============================================
cdef extern from "string.h":
    void *memcpy(void *s1, void *s2, int n)
cdef extern from "stdlib.h":
    void *malloc(int size)
    void free(void *ptr)
    
#==============================================
#External definitions from Sundials headers
#==============================================

cdef extern from "sundials/sundials_nvector.h":
    cdef struct _generic_N_Vector:
        void* content
    ctypedef _generic_N_Vector* N_Vector
    
    
    
cdef extern from "nvector/nvector_serial.h":
    cdef struct _N_VectorContent_Serial:
        long int length
        booleantype own_data
        realtype* data
    ctypedef _N_VectorContent_Serial* N_VectorContent_Serial
    cdef N_Vector N_VMake_Serial(long int vec_length, realtype *v_data)
    N_Vector *N_VCloneVectorArray_Serial(int count, N_Vector w)
    N_Vector *N_VCloneVectorArrayEmpty_Serial(int count, N_Vector w)
    void N_VSetArrayPointer_Serial(realtype *v_data, N_Vector v)
    void N_VConst_Serial(realtype c, N_Vector z)
    N_Vector N_VNew_Serial(long int vec_length)
    void N_VDestroy_Serial(N_Vector v)
    void N_VPrint_Serial(N_Vector v)

#Struct for handling the Jacobian data
cdef extern from "sundials/sundials_direct.h":
    cdef struct _DlsMat:
        int type
        int M
        int N
        int ldim
        int mu
        int ml
        int s_mu
        realtype *data
        int ldata
        realtype **cols
    ctypedef _DlsMat* DlsMat
    cdef realtype* DENSE_COL(DlsMat A, int j)

cdef extern from "cvodes/cvodes.h":
    void* CVodeCreate(int lmm, int iter)
    ctypedef int (*CVRhsFn)(realtype t, N_Vector y, N_Vector ydot, void *f_data)
    int CVodeInit(void *cvode_mem, CVRhsFn f, realtype t0, N_Vector y0)
    int CVodeReInit(void *cvode_mem, realtype t0, N_Vector y0)
    void CVodeFree(void **cvode_mem)
    int CVode(void *cvode_mem, realtype tout, N_Vector yout, realtype *tret, int itask)
    
    #Functions for settings options
    int CVodeSetMaxOrd(void *cvode_mem, int maxord)
    int CVodeSetMaxNumSteps(void *cvode_mem, long int mxsteps)
    int CVodeSetMaxStep(void   *cvode_mem, realtype hmax)
    int CVodeSetInitStep(void  *cvode_mem, realtype hin)
    int CVodeSStolerances(void *cvode_mem, realtype reltol, realtype abstol)
    int CVodeSVtolerances(void *cvode_mem, realtype reltol, N_Vector abstol)
    int CVodeSetStopTime(void  *cvode_mem, realtype tstop)
    int CVodeSetUserData(void  *cvode_mem,void *user_data)
    
    #Functions for retrieving results
    int CVodeGetDky(void *cvode_mem, realtype t, int k, N_Vector dky)
    
    #Functions for error handling
    ctypedef int (*CVErrHandlerFn)(int error_code, char *module, char *function, char *msg,
                                   void *eh_data)
    int CVodeSetErrHandlerFn(void *cvode_mem, CVErrHandlerFn ehfun, void* eh_data)
    
    #Functions for discontinuity handling
    ctypedef int (*CVRootFn)(realtype tt, N_Vector yy, realtype *gout, void *user_data)
    int CVodeRootDirection(void *cvode_mem, int *rootdir)
    int CVodeSetNoInactiveRootWarn(void *cvode_mem)
    int CVodeRootInit(void *cvode_mem, int nrtfn, CVRootFn g)
    int CVodeGetRootInfo(void *cvode_mem, int *rootsfound)
    
    #Functions for retrieving statistics
    int CVodeGetLastOrder(void * cvode_mem,int *qlast)
    int CVodeGetLastStep(void * cvode_mem, realtype *hlast)
    int CVodeGetCurrentOrder(void * cvode_mem,int *qcurrent)
    int CVodeGetNumSteps(void *cvode_mem, long int *nsteps) #Number of steps
    int CVodeGetNumRhsEvals(void *cvode_mem, long int *nrevals) #Number of function evals
    int CVDlsGetNumJacEvals(void *cvode_mem, long int *njevals) #Number of jac evals
    int CVDlsGetNumRhsEvals(void *cvode_mem, long int *nrevalsLS) #Number of res evals due to jac evals
    int CVodeGetNumGEvals(void *cvode_mem, long int *ngevals) #Number of root evals
    int CVodeGetNumErrTestFails(void *cvode_mem, long int *netfails) #Number of local error test failures
    int CVodeGetNumNonlinSolvIters(void *cvode_mem, long int *nniters) #Number of nonlinear iteration
    int CVodeGetNumNonlinSolvConvFails(void *cvode_mem, long int *nncfails) #Number of nonlinear conv failures
    int CVodeGetNonlinSolvStats(void *cvode_mem, long int *nniters, long int *nncfails)
    int CVodeGetIntegratorStats(void* cvode_mem, long int *nsteps, long int *nfevals,
                                long int *nlinsetups, long int *netfails, int *qlast, int *qcur,
                                realtype *hinused, realtype *hlast, realtype *hcur, realtype *tcur)
    
    #Sensitivity methods
    ctypedef int (*CVSensRhsFn)(int Ns, realtype t, N_Vector y, N_Vector ydot, N_Vector *yS,
                                N_Vector *ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2)
    ctypedef int (*CVSensRhs1Fn)(int Ns, realtype t, N_Vector y, N_Vector ydot, int iS, N_Vector *yS,
                                N_Vector *ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2)
    int CVodeSensInit(void *cvode_mem, int Ns, int ism, CVSensRhsFn fS, N_Vector *ySO)
    int CVodeSensInit1(void *cvode_mem, int Ns, int ism, CVSensRhs1Fn fS1, N_Vector *ySO)
    int CVodeSensReInit(void *cvode_mem, int ism, N_Vector *ySO)
    int CVodeSensFree(void *cvode_mem)
    int CVodeSensToggleOff(void *cvode_mem)
    int CVodeSensSStolerances(void *cvode_mem, realtype reltolS, realtype *abstolS)
    int CVodeSensSVtolerances(void *cvode_mem, realtype reltolS, N_Vector *abstolS)
    int CVodeSensEEtolerances(void *cvode_mem)
    int CVodeGetSens(void *cvode_mem, realtype *tret, N_Vector *yS)
    int CVodeGetSensDky(void *cvode_mem, realtype t, int k, N_Vector *dkyS)
    int CVodeGetSens1(void *cvode_mem, realtype *tret, int iss, N_Vector yS)
    int CVodeGetSensDky1(void *cvode_mem, realtype t, int k, int iss, N_Vector dkyS)
    int CVodeSetSensParams(void *cvode_mem, realtype *p, realtype *pbar, int *plist)
    int CVodeSetSensDQMethod(void *cvode_mem, int DQtype, realtype DQrhomax)
    int CVodeSetSensErrCon(void *cvode_mem, booleantype errconS)
    int CVodeSetSensMaxNonlinIters(void *cvode_mem, int maxcorS)
    
    #Statistics
    int CVodeGetEstLocalErrors(void *cvode_mem, N_Vector ele)               #Estimated local errors
    int CVodeGetSensNumRhsEvals(void *cvode_mem, long int *nfSevals)
    int CVodeGetNumRhsEvalsSens(void *cvode_mem, long int *nfevalsS)
    int CVodeGetSensNumErrTestFails(void *cvode_mem, long int *nSetfails)
    int CVodeGetSensNumLinSolvSetups(void *cvode_mem, long int *nlinsetupsS)
    int CVodeGetSensStats(void *cvode_mem, long int *nfSevals, long int *nfevalsS,
                         long int *nSetfails, long int *nlinsetupsS)
    int CVodeGetSensErrWeights(void *cvode_mem, N_Vector *eSweight)
    int CVodeGetSensNumNonlinSolvIters(void *cvode_mem, long int *nSniters)
    int CVodeGetSensNumNonlinSolvConvFails(void *cvode_mem, long int *nSncfails)
    int CVodeGetSensNonlinSolvStats(void *cvode_mem, long int *nSniters, long int *nSncfails)
    int CVodeGetStgrSensNumNonlinSolvIters(void *cvode_mem, long int *nSTGR1niters)
    int CVodeGetStgrSensNumNonlinSolvConvFails(void *cvode_mem, long int *nSTGR1ncfails)
    
    
cdef extern from "cvodes/cvodes_dense.h":
    int CVDense(void *cvode_mem, long int n)
    ctypedef int (*CVDlsDenseJacFn)(int n, realtype t, N_Vector y, N_Vector fy, 
                   DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
    int CVDlsSetDenseJacFn(void *cvode_mem, CVDlsDenseJacFn djac)

cdef extern from "cvodes/cvodes_spgmr.h":
    int CVSpgmr(void *cvode_mem, int pretype, int max1)
    
cdef extern from "cvodes/cvodes_spils.h":
    ctypedef int (*CVSpilsJacTimesVecFn)(N_Vector v, N_Vector Jv, realtype t,
				    N_Vector y, N_Vector fy,
				    void *user_data, N_Vector tmp)
    int CVSpilsSetJacTimesVecFn(void *cvode_mem,  CVSpilsJacTimesVecFn jtv)
    int CVSpilsGetNumJtimesEvals(void *cvode_mem, long int *njvevals) #Number of jac*vector evals
    int CVSpilsGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS) #Number of res evals due to jacÄvector evals
    

cdef extern from "idas/idas.h":
    ctypedef int (*IDAResFn)(realtype tt, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data)
    void* IDACreate()
    int IDAInit(void* ida_mem, IDAResFn res, realtype t0, N_Vector y0, N_Vector yp0)
    int IDAReInit(void* ida_mem, realtype t0, N_Vector y0, N_Vector yp0)
    void IDAFree(void **ida_mem)
    int IDASolve(void* ida_mem, realtype tout,realtype  *tret, N_Vector yret, 
                            N_Vector ypret, int itask)
    
    #Functions for settings options
    int IDASStolerances(void *ida_mem, realtype reltol, realtype abstol)
    int IDASVtolerances(void *ida_mem, realtype reltol, N_Vector abstol)
    int IDASetSuppressAlg(void *ida_mem, booleantype suppressalg)
    int IDASetId(void *ida_mem, N_Vector id)
    int IDASetUserData(void *ida_mem,void *user_data)
    int IDASetInitStep(void *ida_mem, realtype hin)
    int IDASetStopTime(void *ida_mem, realtype tstop)
    int IDASetMaxErrTestFails(void *ida_mem, int maxnef)
    int IDASetMaxNumSteps(void *ida_mem, long int mxsteps)
    int IDASetMaxOrd(void *ida_mem, int maxord)
    int IDASetMaxStep(void* ida_mem, realtype hmax)
    
    #Functions for retrieving results
    int IDAGetDky(void *ida_mem, realtype t, int k, N_Vector dky)
    
    #Functions for error handling
    ctypedef int (*IDAErrHandlerFn)(int error_code, char *module, char *function, char *msg,
                                    void *eh_data)
    int IDASetErrHandlerFn(void *ida_mem,IDAErrHandlerFn ehfun, void* eh_data)
    
    
    #Functions for discontinuity handling
    ctypedef int (*IDARootFn)(realtype tt, N_Vector yy, N_Vector yp, realtype *gout, void *user_data)
    int IDASetRootDirection(void *ida_mem, int *rootdir)
    int IDASetNoInactiveRootWarn(void *ida_mem)
    int IDARootInit(void *ida_mem, int nrtfn, IDARootFn g)
    int IDAGetRootInfo(void *ida_mem, int *rootsfound)
    int IDACalcIC(void *ida_men, int icopt, realtype tout1)
    int IDAGetConsistentIC(void *ida_mem, N_Vector y0, N_Vector yp0)
    int IDASetLineSearchOffIC(void *ida_mem, booleantype lsoff)
    
    #Functions for retrieving statistics
    int IDAGetEstLocalErrors(void *ida_mem, N_Vector ele)               #Estimated local errors
    int IDAGetLastStep(void *ida_mem, realtype *hlast)
    int IDAGetLastOrder(void *ida_mem,int *qlast)                       #Last order used
    int IDAGetCurrentOrder(void *ida_mem,int *qcurrent)                 #Order that is about to be tried
    int IDAGetNumSteps(void *ida_mem, long int *nsteps)                 #Number of steps
    int IDAGetNumResEvals(void *ida_mem, long int *nrevals)             #Number of res evals
    int IDADlsGetNumJacEvals(void *ida_mem, long int *njevals)          #Number of jac evals
    int IDADlsGetNumResEvals(void *ida_mem, long int *nrevalsLS)        #Number of res evals due to jac evals
    int IDAGetNumGEvals(void *ida_mem, long int *ngevals)               #Number of root evals
    int IDAGetNumErrTestFails(void *ida_mem, long int *netfails)        #Number of local error test failures
    int IDAGetNumNonlinSolvIters(void *ida_mem, long int *nniters)      #Number of nonlinear iteration
    int IDAGetNumNonlinSolvConvFails(void *ida_mem, long int *nncfails) #Number of nonlinear conv failures
    int IDAGetNonlinSolvStats(void *ida_mem, long int *nniters, long int *nncfails)
    int IDAGetIntegratorStats(void* ida_mem,long int  *nsteps, long int *nrevals, 
                            long int *nlinsetups, long int *netfails, int *klast, 
                            int *kcur, realtype *hinused, realtype *hlast, 
                            realtype *hcur, realtype *tcur)
    
    #Start Sensitivities
    #===================
    ctypedef int (*IDASensResFn)(int Ns, realtype t, N_Vector yy, N_Vector yp, 
                                N_Vector *yS, N_Vector *ypS, N_Vector *resvalS, 
                                void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
    int IDASensInit(void *ida_mem, int Ns, int ism, IDASensResFn resS, N_Vector *ySO, N_Vector *ypSO)
    int IDASensReInit(void *ida_mem, int ism, N_Vector *ySO, N_Vector *ypSO)
    
    #Options
    int IDASensToggleOff(void *ida_mem)
    int IDASensSStolerances(void *ida_mem, realtype reltolS, realtype *abstolS)
    int IDASensSVtolerances(void *ida_mem, realtype reltolS, N_Vector *abstolS)
    int IDASensEEtolerances(void *ida_mem)
    
    #Results
    int IDAGetSens(void *ida_mem, realtype tret, N_Vector *yS)
    int IDAGetSensDky(void *ida_mem, realtype t, int k, N_Vector *dkyS)
    int IDAGetSensDky1(void *ida_mem, realtype t, int k, int i, N_Vector dkyS)
    
    #Options (optional)
    int IDASetSensParams(void *ida_mem, realtype *p, realtype *pbar, int *plist)
    int IDASetSensDQMethod(void *ida_mem, int DQtype, realtype DQrhomax)
    int IDASetSensErrCon(void *ida_mem, booleantype errconS)
    int IDASetSensMaxNonlinIters(void *ida_mem, int maxcorS)
    
    #Statistics
    int IDAGetSensNumResEvals(void *ida_mem, long int nfSevals)
    int IDAGetNumResEvalsSens(void *ida_mem, long int nfevalsS)
    int IDAGetSensNumErrTestFails(void *ida_mem, long int nSetfails)
    int IDAGetSensNumLinSolvSetups(void *ida_mem, long int nlinsetupsS)
    int IDAGetSensStats(void *ida_mem, long int *nfSevals, long int *nfevalsS, 
                        long int *nSetfails, long int *nlinsetupsS)
    int IDAGetSensNumNonlinSolvIters(void *ida_mem, long int nSniters)
    int IDAGetSeonsNumNonlinSolvConvFails(void *ida_mem, long int nSncfails)
    int IDAGetSensNonlinSolvStats(void *ida_mem, long int *nSniters, long int *nSncfails)
    
    #End Sensitivities
    #=================

cdef extern from "idas/idas_dense.h":
    int IDADense(void *ida_mem, long int n)
    ctypedef int (*IDADlsDenseJacFn)(int Neq, realtype tt, realtype cj, N_Vector yy, 
                   N_Vector yp, N_Vector rr, DlsMat Jac, void *user_data, 
                   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
    int IDADlsSetDenseJacFn(void *ida_mem, IDADlsDenseJacFn djac)

#=========================
# END SUNDIALS DEFINITIONS
#=========================

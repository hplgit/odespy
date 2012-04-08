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

import numpy as np
from numpy cimport NPY_DOUBLE, npy_intp, NPY_INT

#==============================================
#External definitions from numpy headers
#==============================================
cdef extern from "numpy/arrayobject.h":

    ctypedef int npy_intp 

    ctypedef extern class numpy.dtype [object PyArray_Descr]:
        cdef int type_num, elsize, alignment
        cdef char type, kind, byteorder, hasobject
        cdef object fields, typeobj
        
    ctypedef int intp
    ctypedef extern class numpy.ndarray [object PyArrayObject]:
        cdef char *data
        cdef int nd
        cdef intp *dimensions
        cdef intp *strides
        cdef int flags
        
    cdef object PyArray_SimpleNew(int nd, npy_intp* dims, int typenum)
    cdef object PyArray_SimpleNewFromData(int nd, npy_intp *dims,
                                           int typenum, void *data)
    void import_array() 
    void* PyArray_GetPtr(ndarray aobj, npy_intp* ind)
    void *PyArray_DATA(ndarray aobj)

import_array()

#==============================================
# C headers
#==============================================
cdef extern from "string.h":
    void *memcpy(void *s1, void *s2, int n)
cdef extern from "stdio.h":
    int printf (char *__format, ...)
cdef extern from "stdlib.h":
    ctypedef unsigned long size_t
    void *calloc(size_t nelem, size_t elsize)

#==============================================
#External definitions from Sundials headers
#==============================================

# import data types
# basic data types
cdef extern from "sundials/sundials_types.h":
    ctypedef double realtype
    ctypedef bint booleantype

# sundials handling of vectors, so far only serial,
# parallel to be implemented later?
cdef extern from "sundials/sundials_nvector.h":
    cdef struct _generic_N_Vector:
        void* content
    ctypedef _generic_N_Vector* N_Vector

cdef extern from "nvector/nvector_serial.h":
    cdef struct _N_VectorContent_Serial:
        long int length
        realtype* data
    ctypedef _N_VectorContent_Serial* N_VectorContent_Serial
    N_Vector N_VNew_Serial(long int vec_length)
    void N_VDestroy_Serial(N_Vector v)
    void N_VPrint_Serial(N_Vector v)
    cdef N_Vector N_VMake_Serial(long int vec_length, realtype *v_data)
    N_Vector *N_VCloneVectorArray_Serial(int count, N_Vector w)
    N_Vector *N_VCloneVectorArrayEmpty_Serial(int count, N_Vector w)
    void N_VSetArrayPointer_Serial(realtype *v_data, N_Vector v)
    void N_VConst_Serial(realtype c, N_Vector z)

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

# KINSOL functions and routines
cdef extern from "kinsol/kinsol.h":
    # user defined functions
    ctypedef int (*KINSysFn)(N_Vector uu, N_Vector fval, void *user_data )
    ctypedef void (*KINErrHandlerFn)(int error_code, char *module, char *function, char *msg, void *user_data)
    ctypedef void (*KINInfoHandlerFn)(char *module, char *function, char *msg, void *user_data)
    # initialization routines
    void *KINCreate()
    int KINInit(void *kinmem, KINSysFn func, N_Vector tmpl)

    # optional input spec. functions,
    # for specificationsdocumentation cf. kinsol.h line 218-449
    int KINSetErrHandlerFn(void *kinmem, KINErrHandlerFn ehfun, void *eh_data)
    int KINSetInfoHandlerFn(void *kinmem, KINInfoHandlerFn ihfun, void *ih_data)
    int KINSetUserData(void *kinmem, void *user_data)
    int KINSetPrintLevel(void *kinmemm, int printfl)
    int KINSetNumMaxIters(void *kinmem, long int mxiter)
    int KINSetNoInitSetup(void *kinmem, booleantype noInitSetup)
    int KINSetNoResMon(void *kinmem, booleantype noNNIResMon)
    int KINSetMaxSetupCalls(void *kinmem, long int msbset)
    int KINSetMaxSubSetupCalls(void *kinmem, long int msbsetsub)
    int KINSetEtaForm(void *kinmem, int etachoice)
    int KINSetEtaConstValue(void *kinmem, realtype eta)
    int KINSetEtaParams(void *kinmem, realtype egamma, realtype ealpha)
    int KINSetResMonParams(void *kinmem, realtype omegamin, realtype omegamax)
    int KINSetResMonConstValue(void *kinmem, realtype omegaconst)
    int KINSetNoMinEps(void *kinmem, booleantype noMinEps)
    int KINSetMaxNewtonStep(void *kinmem, realtype mxnewtstep)
    int KINSetMaxBetaFails(void *kinmem, long int mxnbcf)
    int KINSetRelErrFunc(void *kinmem, realtype relfunc)
    int KINSetFuncNormTol(void *kinmem, realtype fnormtol)
    int KINSetScaledStepTol(void *kinmem, realtype scsteptol)
    int KINSetConstraints(void *kinmem, N_Vector constraints)
    int KINSetSysFunc(void *kinmem, KINSysFn func)

    # solver routine
    int KINSol(void *kinmem, N_Vector uu, int strategy, N_Vector u_scale, N_Vector f_scale)

    # optional output routines.
    # Documentation see kinsol.h line 670-735
    int KINGetWorkSpace(void *kinmem, long int *lenrw, long int *leniw)
    int KINGetNumNonlinSolvIters(void *kinmem, long int *nniters)
    int KINGetNumFuncEvals(void *kinmem, long int *nfevals)
    int KINGetNumBetaCondFails(void *kinmem, long int *nbcfails) 
    int KINGetNumBacktrackOps(void *kinmem, long int *nbacktr)
    int KINGetFuncNorm(void *kinmem, realtype *fnorm)
    int KINGetStepLength(void *kinmem, realtype *steplength)
    char *KINGetReturnFlagName(int flag)

    # fuction used to deallocate memory used by KINSOL
    void KINFree(void **kinmem)

# linear solver modules
cdef extern from "kinsol/kinsol_dense.h":
    int KINDense(void *kinmem, int N)

cdef extern from "kinpinv.h":
    int KINPinv(void *kinmem, int N)

# functions used for supplying jacobian, and receiving info from linear solver
cdef extern from "kinsol/kinsol_direct.h":
    # user functions
    ctypedef int (*KINDlsDenseJacFn)(int N, N_Vector u, N_Vector fu, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2)
    
    # function used to link user functions to KINSOL
    int KINDlsSetDenseJacFn(void *kinmem, KINDlsDenseJacFn jac)

    # optional output fcts for linear direct solver
    int KINDlsGetWorkSpace(void *kinmem, long int *lenrwB, long int *leniwB)
    int KINDlsGetNumJacEvals(void *kinmem, long int *njevalsB)
    int KINDlsGetNumFuncEvals(void *kinmem, long int *nfevalsB)
    int KINDlsGetLastFlag(void *kinmem, int *flag)
    char *KINDlsGetReturnFlagName(int flag)

cdef extern from "kinsol_jmod.h":
    # user functions
    ctypedef int (*KINPinvJacFn)(int N, N_Vector u, N_Vector fu, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2)

    # function used to link user jacobian to KINSOL
    int KINPinvSetJacFn(void *kinmem, KINPinvJacFn jac)
    
    # functions used to set regularizationparameter
    int KINPinvSetRegParam(void *kinmem, realtype reg_p)

    # optional output fcts for linear direct solver
    int KINPinvGetWorkSpace(void *kinmem, long int *lenrwB, long int *leniwB)
    int KINPinvGetNumJacEvals(void *kinmem, long int *njevalsB)
    int KINPinvGetNumFuncEvals(void *kinmem, long int *nfevalsB)
    int KINPinvGetLastFlag(void *kinmem, int *flag)
    char *KINPinvGetReturnFlagName(int flag)
    
    
#=========================
# END SUNDIALS DEFINITIONS
#=========================

#===========================
# Problem information struct
#===========================
cdef class ProblemData:
    cdef:
        void *RHS          # Residual or right-hand-side
        void *JAC
        int dim            # Dimension of the problem

#=================
# Module functions
#=================

cdef inline N_Vector arr2nv(x):
    x=np.array(x)
    cdef long int n = len(x)
    cdef ndarray[double, ndim=1,mode='c'] ndx=x
    import_array()
    cdef void* data_ptr=PyArray_DATA(ndx)
    cdef N_Vector v=N_VNew_Serial(n)
    memcpy((<N_VectorContent_Serial>v.content).data, data_ptr, n*sizeof(double))
    return v
    
cdef inline nv2arr(N_Vector v):
    cdef long int n = (<N_VectorContent_Serial>v.content).length
    cdef realtype* v_data = (<N_VectorContent_Serial>v.content).data
    cdef long int i
    import_array()
    cdef ndarray x=np.empty(n)
    #cdef ndarray x = PyArray_SimpleNewFromData(1, &dims, NPY_DOUBLE,v_data)
    memcpy(x.data, v_data, n*sizeof(double))
    return x

cdef inline realtype2arr(realtype *data, int n):
    """ Create new numpy array from realtype*"""
    import_array()
    cdef ndarray[realtype, ndim=1, mode='c'] x=np.empty(n)
    memcpy(x.data, data, n*sizeof(double))
    return x

cdef inline arr2realtype(ndarray x, realtype* out_data, int n):
    """ Create new realtype* from numpy """
    cdef ndarray[double, ndim=1,mode='c'] ndx=x
    import_array()
    cdef void* data_ptr=PyArray_DATA(ndx)
    memcpy(out_data, data_ptr, n*sizeof(double))


cdef inline arr2int(ndarray x, int* out_data, int n):
    """ Create new int* from numpy array """
    cdef ndarray[int, ndim=1,mode='c'] ndx=x
    import_array()
    cdef void* data_ptr=PyArray_DATA(ndx)
    memcpy(out_data, data_ptr, n*sizeof(int))


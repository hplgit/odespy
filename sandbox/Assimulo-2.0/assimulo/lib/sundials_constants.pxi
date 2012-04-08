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
This file contains the constants from Sundials.
"""

#==========
# CVode
#==========

#CVode Input
#---------
# Main solver module
DEF CV_ADAMS              = 1
DEF CV_BDF                = 2
DEF CV_FUNCTIONAL         = 1
DEF CV_NEWTON             = 2
DEF CV_NORMAL             = 1
DEF CV_ONE_STEP           = 2
DEF CV_SIMULTANEOUS       = 1
DEF CV_STAGGERED          = 2
DEF CV_STAGGERED1         = 3
DEF CV_CENTERED           = 1
DEF CV_FORWARD            = 2
# Iterative solver module
DEF PREC_NONE             = 0
DEF PREC_LEFT             = 1
DEF PREC_RIGHT            = 2
DEF PREC_BOTH             = 3

#CVode Output
#----------
# Main solver module
DEF CV_SUCCESS            = 0
DEF CV_TSTOP_RETURN       = 1
DEF CV_REC_ERR            = 1   # Recoverable error.
DEF CV_ROOT_RETURN        = 2   # CVSolve succeeded and found one or more roots.
DEF CV_WARNING            = 99
DEF CV_TOO_MUCH_WORK      = -1
DEF CV_TOO_MUCH_ACC       = -2
DEF CV_ERR_FAIL           = -3
DEF CV_CONV_FAIL          = -4
DEF CV_LINIT_FAIL         = -5
DEF CV_LSETUP_FAIL        = -6
DEF CV_LSOLVE_FAIL        = -7
DEF CV_RHSFUNC_FAIL       = -8
DEF CV_FIRST_RHSFUNC_ERR  = -9
DEF CV_REPTD_RHSFUNC_ERR  = -10
DEF CV_UNREC_RHSFUNC_ERR  = -11
DEF CV_RTFUNC_FAIL        = -12 # The rootfinding function failed in an unrecoverable manner.
DEF CV_MEM_FAIL           = -20
DEF CV_MEM_NULL           = -21
DEF CV_ILL_INPUT          = -22
DEF CV_NO_MALLOC          = -23
DEF CV_BAD_K              = -24
DEF CV_BAD_T              = -25
DEF CV_BAD_DKY            = -26
DEF CV_TOO_CLOSE          = -27
# Linear solver module
DEF CVDLS_SUCCESS         = 0  #Successful function return
DEF CVDLS_MEM_NULL        = -1 #The cvode_mem argument was NULL
DEF CVDLS_LMEM_NULL       = -2
DEF CVDLS_ILL_INPUT       = -3
DEF CVDLS_MEM_FAIL        = -4 #A memory allocation request failed.
DEF CVDLS_JACFUNC_UNRECVR = -5
DEF CVDLS_JACFUNC_RECVR   = -6
# Iterative
DEF SPGMR_SUCCESS         = 0 #Success, Converged
DEF SPGMR_ATIMES_FAIL_REC = 5 #The Jacobian-vector function failed recoverably

# Sensitivity constants
DEF CV_SRHSFUNC_FAIL      = -41   #The sensitivity right-hand side function failed unrecoverable
DEF CV_FIRST_SRHSFUNC_ERR = -42   #The sensitivity right-hand side function failed at first call
DEF CV_REPTD_SRHSFUNC_ERR = -43   #Covergence tests occurred too many times.
DEF CV_UNREC_SRHSFUNC_ERR = -44   #The sensitivity right-hand side function had a recoverable error but was unable to recover
DEF CV_BAD_IS             = -45


#==========
# IDA
#==========

#IDA Input
#---------
# Main solver module
DEF IDA_NORMAL            = 1   # Solver returns at specified output time.
DEF IDA_ONE_STEP          = 2   # Solver returns after each successful step.
DEF IDA_SIMULTANEOUS      = 1   # Simultaneous corrector forward sensitivity method.
DEF IDA_STAGGERED         = 2   # Staggered corrector forward sensitivity method.
DEF IDA_CENTERED          = 1   # Central difference quotient approximation (2nd order) of the sensitivity RHS.
DEF IDA_FORWARD           = 2   # Forward difference quotient approximation (1st order) of the sensitivity RHS.
DEF IDA_YA_YDP_INIT       = 1   # See IDA Documentation 4.5.4
DEF IDA_Y_INIT            = 2   # See IDA Documentation 4.5.4
# Iterative solver module
DEF PREC_NONE             = 0
DEF PREC_LEFT             = 1
DEF MODIFIED_GS           = 1
DEF CLASSICAL_GS          = 2

#IDA Output
#----------
# Main solver module
DEF IDA_SUCCESS           = 0 # Successful function return.
DEF IDA_REC_ERR           = 1 # Recoverable error.
DEF IDA_TSTOP_RETURN      = 1 # IDASolve succeeded by reaching the specified stopping point.
DEF IDA_ROOT_RETURN       = 2 # IDASolve succeeded and found one or more roots.
DEF IDA_WARNING           = 99
DEF IDA_TOO_MUCH_WORK     = -1
DEF IDA_TOO_MUCH_ACC      = -2
DEF IDA_ERR_FAIL          = -3
DEF IDA_CONV_FAIL         = -4
DEF IDA_LINIT_FAIL        = -5
DEF IDA_LSETUP_FAIL       = -6
DEF IDA_LSOLVE_FAIL       = -7
DEF IDA_RES_FAIL          = -8
DEF IDA_REP_RES_FAIL      = -9
DEF IDA_RTFUNC_FAIL       = -10 # The rootfinding function failed in an unrecoverable manner.
DEF IDA_CONSTR_FAIL       = -11
DEF IDA_FIRST_RES_FAIL    = -12
DEF IDA_LINESEARCH_FAIL   = -13
DEF IDA_NO_RECOVERY       = -14
DEF IDA_MEM_NULL          = -20
DEF IDA_MEM_FAIL          = -21
DEF IDA_ILL_INPUT         = -22
DEF IDA_NO_MALLOC         = -23
DEF IDA_BAD_EWT           = -24
DEF IDA_BAD_K             = -25
DEF IDA_BAD_T             = -26
DEF IDA_BAD_DKY           = -27
DEF IDA_NO_QUAD           = -30
DEF IDA_QRHS_FAIL         = -31
DEF IDA_FIRST_QRHS_ERR    = -32
DEF IDA_REP_QRHS_ERR      = -33
DEF IDA_NO_SENS           = -40
DEF IDA_SRES_FAIL         = -41
DEF IDA_REP_SRES_ERR      = -42
DEF IDA_BAD_IS            = -43
DEF IDA_NO_QUADSENS       = -50
DEF IDA_QSRHS_FAIL        = -51
DEF IDA_FIRST_QSRHS_ERR   = -52
DEF IDA_REP_QSRHS_ERR     = -53
# Linear solver module
DEF IDADLS_SUCCESS         = 0
DEF IDADLS_MEM_NULL        = -1
DEF IDADLS_LMEM_NULL       = -2
DEF IDADLS_ILL_INPUT       = -3
DEF IDADLS_MEM_FAIL        = -4
DEF IDADLS_JACFUNC_UNRECVR = -5
DEF IDADLS_JACFUNC_RECVR   = -6
DEF IDADLS_NO_ADJ          = -101
DEF IDADLS_LMEMB_NULL      = -102
DEF IDASPILS_SUCCESS       = 0
DEF IDASPILS_MEM_NULL      = -1
DEF IDASPILS_LMEM_NULL     = -2
DEF IDASPILS_ILL_INPUT     = -3
DEF IDASPILS_MEM_FAIL      = -4
DEF IDASPILS_PMEM_NULL     = -5
DEF IDASPILS_NO_ADJ        = -101
DEF IDASPILS_LMEMB_NULL    = -102




DEF CV_RHS_IND        = 0   # Index to user data rhs handling
DEF CV_RHSF_IND       = 0   # Index to user data rhs
DEF CV_JAC_IND        = 1   # Index to user data jacobian
DEF CV_ROOT_IND       = 1   # Index to user data root handling
DEF CV_ROOTF_IND      = 0   # Index to user data root function
DEF CV_SW_IND         = 1   # Index to user data root switches

DEF IDA_NORMAL         = 1   # Solver returns at specified output time.
DEF IDA_ONE_STEP       = 2   # Solver returns after each successful step.
DEF IDA_RES_IND        = 0   # Index to user data residual handling
DEF IDA_SENS_IND       = 0   # Index to indicator for sensitivites
DEF IDA_RESF_IND       = 1   # Index to user data residual
DEF IDA_JAC_IND        = 2   # Index to user data jacobian
DEF IDA_ROOT_IND       = 1   # Index to user data root handling
DEF IDA_ROOTF_IND      = 0   # Index to user data root function
DEF IDA_SW_IND         = 1   # Index to user data root switches


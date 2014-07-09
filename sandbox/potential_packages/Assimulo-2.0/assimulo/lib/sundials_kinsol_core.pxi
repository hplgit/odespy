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

This file contains the constants from Sundials.
"""
# -----------------------------------------------------------------
# KINSOL return flags
# -----------------------------------------------------------------

DEF KIN_SUCCESS             =0
DEF KIN_INITIAL_GUESS_OK    =1
DEF KIN_STEP_LT_STPTOL      =2

DEF KIN_REC_ERR             =5

DEF KIN_WARNING             =99

DEF KIN_MEM_NULL            =-1
DEF KIN_ILL_INPUT           =-2
DEF KIN_NO_MALLOC           =-3
DEF KIN_MEM_FAIL            =-4
DEF KIN_LINESEARCH_NONCONV  =-5
DEF KIN_MAXITER_REACHED     =-6
DEF KIN_MXNEWT_5X_EXCEEDED  =-7
DEF KIN_LINESEARCH_BCFAIL   =-8
DEF KIN_LINSOLV_NO_RECOVERY =-9
DEF KIN_LINIT_FAIL          =-10
DEF KIN_LSETUP_FAIL         =-11
DEF KIN_LSOLVE_FAIL         =-12

DEF KIN_SYSFUNC_FAIL        =-13
DEF KIN_FIRST_SYSFUNC_ERR   =-14
DEF KIN_REPTD_SYSFUNC_ERR   =-15


# -----------------------------------------------------------------
# Enumeration for inputs to KINSetEtaForm (eta choice)
# -----------------------------------------------------------------
# KIN_ETACONSTANT : use constant value for eta (default value is
#                   0.1 but a different value can be specified via
#                   a call to KINSetEtaConstValue)
#
# KIN_ETACHOICE1 : use choice #1 as given in Eisenstat and Walker's
#                  paper of SIAM J.Sci.Comput.,17 (1996), pp 16-32,
#                  wherein eta is defined to be:
#
#              eta(k+1) = ABS(||F(u_k+1)||_L2-||F(u_k)+J(u_k)*p_k||_L2)
#                       ---------------------------------------------
#                                       ||F(u_k)||_L2
#
#                                                      1+sqrt(5)
#              eta_safe = eta(k)^ealpha where ealpha = ---------
#                                                          2
#
# KIN_ETACHOICE2 : use choice #2 as given in Eisenstat and Walker's
#                  paper wherein eta is defined to be:
#
#                                  [ ||F(u_k+1)||_L2 ]^ealpha
#              eta(k+1) = egamma * [ --------------- ]
#                                  [  ||F(u_k)||_L2  ]
#
#                  where egamma = [0,1] and ealpha = (1,2]
#
#              eta_safe = egamma*(eta(k)^ealpha)
#
#                  Note: The default values of the scalar
#                  coefficients egamma and ealpha (both required)
#                  are egamma = 0.9 and ealpha = 2.0, but the
#                  routine KINSetEtaParams can be used to specify
#                  different values.
#
# When using either KIN_ETACHOICE1 or KIN_ETACHOICE2, if
# eta_safe > 0.1 then the following safeguard is applied:
#
#  eta(k+1) = MAX {eta(k+1), eta_safe}
#
# The following safeguards are always applied when using either
# KIN_ETACHOICE1 or KIN_ETACHOICE2 so that eta_min <= eta <= eta_max:
#
#  eta(k+1) = MAX {eta(k+1), eta_min}
#  eta(k+1) = MIN {eta(k+1), eta_max}
#
# where eta_min = 1.0e-4 and eta_max = 0.9 (see KINForcingTerm).
# -----------------------------------------------------------------

DEF KIN_ETACHOICE1  =1
DEF KIN_ETACHOICE2  =2
DEF KIN_ETACONSTANT =3

# -----------------------------------------------------------------
# Enumeration for global strategy
# -----------------------------------------------------------------
# Choices are KIN_NONE and KIN_LINESEARCH.
# -----------------------------------------------------------------
  
DEF KIN_NONE       =0
DEF KIN_LINESEARCH =1

# -----------------------------------------------------------------
# KINDirect constants
# -----------------------------------------------------------------

DEF KINDLS_SUCCESS           =0

DEF KINDLS_MEM_NULL         =-1
DEF KINDLS_LMEM_NULL        =-2
DEF KINDLS_ILL_INPUT        =-3
DEF KINDLS_MEM_FAIL         =-4
DEF KINDLS_JACFUNC_UNRECVR  =-5
DEF KINDLS_JACFUNC_RECVR    =-6

 
# -----------------------------------------------------------------
# KINPinv constants
# -----------------------------------------------------------------

DEF KINPINV_SUCCESS         =  0
DEF KINPINV_MEM_NULL        = -1
DEF KINPINV_LMEM_NULL       = -2
DEF KINPINV_ILL_INPUT       = -3
DEF KINPINV_MEM_FAIL        = -4
DEF KINPINV_JACFUNC_UNRECVR = -5
DEF KINPINV_JACFUNC_RECVR   = -6

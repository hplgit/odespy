/*
 * Copyright (C) 2010 Modelon AB / Copyright (c) 2002, The Regents of the University of California.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * This file is a modification of the file kinsol_impl.h , revision 1.5,
 * from the SUNDIALS suite. The file is modified by:
 *
 * Johan Ylikiiskilä - johan.ylikiiskila@gmail.com
 *
 */

#ifndef _KINJMOD_IMPL_H
#define _KINJMOD_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "kinsol_jmod.h"


/*
 * -----------------------------------------------------------------
 * Types: KINPinvMemRec, KINPinvMem                             
 * -----------------------------------------------------------------
 * The type KINPinvMem is pointer to a KINPinvMemRec.
 * This structure contains KINPinv solver-specific data. 
 * -----------------------------------------------------------------
 */

typedef struct KINPinvMemRec {

  int d_type;              /* SUNDIALS_DENSE or SUNDIALS_BAND              */

  int d_n;                 /* problem dimension                            */

  int d_ml;                /* lower bandwidth of Jacobian                  */
  int d_mu;                /* upper bandwidth of Jacobian                  */ 
  int d_smu;               /* upper bandwith of M = MIN(N-1,d_mu+d_ml)     */

  booleantype d_jacDQ;     /* TRUE if using internal DQ Jacobian approx.   */
  KINPinvJacFn d_djac;     /* dense Jacobian routine to be called          */

  void *d_J_data;          /* J_data is passed to djac or bjac             */
    
  DlsMat d_J;              /* problem Jacobian                             */

  int *d_pivots;           /* pivot array for PM = LU                      */
  realtype *d_beta;
  realtype d_reg_param;    /* Regularization parameter                     */
  long int d_nje;          /* no. of calls to jac                          */
    
  long int d_nfeDQ;        /* no. of calls to F due to DQ Jacobian approx. */
    
  int d_last_flag;         /* last error return flag                       */

  DlsMat d_JTJ;
  booleantype d_regularized; /* Boolean set to true if problem is regularized*/
  booleantype d_redojac;
    
} *KINPinvMem;


/*
 * -----------------------------------------------------------------
 * Prototypes of internal functions
 * -----------------------------------------------------------------
 */

int kinPinvDQJac(int N,
		 N_Vector u, N_Vector fu,
		 DlsMat Jac, void *data,
		 N_Vector tmp1, N_Vector tmp2);

/*
 * -----------------------------------------------------------------
 * Error Messages
 * -----------------------------------------------------------------
 */

#define MSGD_KINMEM_NULL "KINSOL memory is NULL."
#define MSGD_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGD_MEM_FAIL    "A memory request failed."
#define MSGD_LMEM_NULL   "Linear solver memory is NULL."
#define MSGD_BAD_SIZES   "Illegal bandwidth parameter(s). Must have 0 <=  ml, mu <= N-1."
#define MSGD_JACFUNC_FAILED "The Jacobian routine failed in an unrecoverable manner."

#ifdef __cplusplus
}
#endif

#endif

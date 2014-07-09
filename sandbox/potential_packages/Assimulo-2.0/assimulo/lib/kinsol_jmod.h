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
 * This file is a modification of the file kinsol_direct.h , revision 1.5,
 * from the SUNDIALS suite. The file is modified by:
 *
 * Johan Ylikiiskilä - johan.ylikiiskila@gmail.com
 *
 */


#ifndef _KINJMOD_H
#define _KINJMOD_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <sundials/sundials_direct.h>
#include <sundials/sundials_nvector.h>

/*
 * =================================================================
 *              K I N P I N V    C O N S T A N T S
 * =================================================================
 */

/* 
 * -----------------------------------------------------------------
 * KINPINV return values 
 * -----------------------------------------------------------------
 */

#define KINPINV_SUCCESS           0
#define KINPINV_MEM_NULL         -1
#define KINPINV_LMEM_NULL        -2
#define KINPINV_ILL_INPUT        -3
#define KINPINV_MEM_FAIL         -4

/* Additional last_flag values */

#define KINPINV_JACFUNC_UNRECVR  -5
#define KINPINV_JACFUNC_RECVR    -6

/*
 * =================================================================
 *              F U N C T I O N   T Y P E S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Type: KINPinvJacFn
 * -----------------------------------------------------------------
 *
 * A dense Jacobian approximation function Jac must be of type 
 * KINDlsDenseJacFn. Its parameters are:
 *
 * N        - problem size.
 *
 * u        - current iterate (unscaled) [input]
 *
 * fu       - vector (type N_Vector) containing result of nonlinear
 *            system function evaluated at current iterate:
 *            fu = F(u) [input]
 *
 * J        - dense matrix (of type DlsMat) that will be loaded
 *            by a KINDlsDenseJacFn with an approximation to the
 *            Jacobian matrix J = (dF_i/dy_j).
 *
 * user_data   - pointer to user data - the same as the user_data
 *            parameter passed to KINSetFdata.
 *
 * tmp1, tmp2 - available scratch vectors (volatile storage)
 *
 * A KINPinvJacFn should return 0 if successful, a positive 
 * value if a recoverable error occurred, and a negative value if 
 * an unrecoverable error occurred.
 *
 * -----------------------------------------------------------------
 *
 * NOTE: The following are two efficient ways to load a dense Jac:         
 * (1) (with macros - no explicit data structure references)      
 *     for (j=0; j < Neq; j++) {                                  
 *       col_j = DENSE_COL(Jac,j);                                 
 *       for (i=0; i < Neq; i++) {                                
 *         generate J_ij = the (i,j)th Jacobian element           
 *         col_j[i] = J_ij;                                       
 *       }                                                        
 *     }                                                          
 * (2) (without macros - explicit data structure references)      
 *     for (j=0; j < Neq; j++) {                                  
 *       col_j = (Jac->data)[j];                                   
 *       for (i=0; i < Neq; i++) {                                
 *         generate J_ij = the (i,j)th Jacobian element           
 *         col_j[i] = J_ij;                                       
 *       }                                                        
 *     }                                                          
 * A third way, using the DENSE_ELEM(A,i,j) macro, is much less   
 * efficient in general.  It is only appropriate for use in small 
 * problems in which efficiency of access is NOT a major concern. 
 *                                                                
 * -----------------------------------------------------------------
 */
  
  
typedef int (*KINPinvJacFn)(int N,
			    N_Vector u, N_Vector fu, 
			    DlsMat J, void *user_data,
			    N_Vector tmp1, N_Vector tmp2);
  

/*
 * =================================================================
 *            E X P O R T E D    F U N C T I O N S 
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Optional inputs to the KINPinv linear solver
 * -----------------------------------------------------------------
 *
 * KINDlsSetDenseJacFn specifies the dense Jacobian approximation
 * routine to be used for a direct dense linear solver.
 *
 * By default, a difference quotient approximation, supplied with
 * the solver is used.
 *
 * The return value is one of:
 *    KINPINV_SUCCESS   if successful
 *    KINPINV_MEM_NULL  if the KINSOL memory was NULL
 *    KINPINV_LMEM_NULL if the linear solver memory was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int KINPinvSetJacFn(void *kinmem, KINPinvJacFn jac);
SUNDIALS_EXPORT int KINPinvSetRegParam(void *kinmem, realtype reg_p);

/*
 * -----------------------------------------------------------------
 * Optional outputs from a KINDLS linear solver
 * -----------------------------------------------------------------
 *
 * KINPinvGetWorkSpace    returns the real and integer workspace used
 *                       by the KINDLS linear solver.
 * KINPinvGetNumJacEvals  returns the number of calls made to the
 *                       Jacobian evaluation routine.
 * KINPinvGetNumFuncEvals returns the number of calls to the user's F
 *                       routine due to finite difference Jacobian
 *                       evaluation.
 * KINPinvGetLastFlag     returns the last error flag set by any of
 *                       the KINDLS interface functions.
 * KINPinvGetReturnFlagName returns the name of the constant
 *                       associated with a KINDLS return flag
 *
 * The return value of KINPinvGet* is one of:
 *    KINPINV_SUCCESS   if successful
 *    KINPINV_MEM_NULL  if the KINSOL memory was NULL
 *    KINPINV_LMEM_NULL if the linear solver memory was NULL
 * -----------------------------------------------------------------
 */


SUNDIALS_EXPORT int KINPinvGetWorkSpace(void *kinmem, long int *lenrwB, long int *leniwB);
SUNDIALS_EXPORT int KINPinvGetNumJacEvals(void *kinmem, long int *njevalsB);
SUNDIALS_EXPORT int KINPinvGetNumFuncEvals(void *kinmem, long int *nfevalsB);
SUNDIALS_EXPORT int KINPinvGetLastFlag(void *kinmem, int *flag);
SUNDIALS_EXPORT char *KINPinvGetReturnFlagName(int flag);

#ifdef __cplusplus
}
#endif

#endif

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
 * This file is a modification of the file kinsol_dense.c , revision 1.5,
 * from the SUNDIALS suite. The file is modified by:
 *
 * Johan Ylikiiskilä - johan.ylikiiskila@gmail.com
 *
 */

/*
 *This is a linear solver using Tikhonov regularization
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "kinsol_jmod_impl.h"
#include "kinpinv.h"
#include "kinsol/kinsol_impl.h"

#include <sundials/sundials_direct.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>

/* Constants */

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* 
 * =================================================================
 * PROTOTYPES FOR PRIVATE FUNCTIONS
 * =================================================================
 */


/* help functions */
void regMatrix(realtype **JTJ_c, realtype **jac, realtype h,int size);

/* KINPinv linit, lsetup, lsolve, and lfree routines */

static int kinPinvInit(KINMem kin_mem);
static int kinPinvSetup(KINMem kin_mem);
static int kinPinvSolve(KINMem kin_mem, N_Vector x, N_Vector b,
                         realtype *res_norm);
static void kinPinvFree(KINMem kin_mem);


/*
 * =================================================================
 * READIBILITY REPLACEMENTS
 * =================================================================
 */

#define lrw1           (kin_mem->kin_lrw1)
#define liw1           (kin_mem->kin_liw1)
#define func           (kin_mem->kin_func)
#define printfl        (kin_mem->kin_printfl)
#define linit          (kin_mem->kin_linit)
#define lsetup         (kin_mem->kin_lsetup)
#define lsolve         (kin_mem->kin_lsolve)
#define lfree          (kin_mem->kin_lfree)
#define lmem           (kin_mem->kin_lmem)
#define inexact_ls     (kin_mem->kin_inexact_ls)
#define uu             (kin_mem->kin_uu)
#define fval           (kin_mem->kin_fval)
#define uscale         (kin_mem->kin_uscale)
#define fscale         (kin_mem->kin_fscale)
#define sqrt_relfunc   (kin_mem->kin_sqrt_relfunc)
#define sJpnorm        (kin_mem->kin_sJpnorm)
#define sfdotJp        (kin_mem->kin_sfdotJp)
#define errfp          (kin_mem->kin_errfp)
#define infofp         (kin_mem->kin_infofp)
#define setupNonNull   (kin_mem->kin_setupNonNull)
#define vtemp1         (kin_mem->kin_vtemp1)
#define vec_tmpl       (kin_mem->kin_vtemp1)
#define vtemp2         (kin_mem->kin_vtemp2)

#define mtype          (kinpinv_mem->d_type)
#define n              (kinpinv_mem->d_n)
#define ml             (kinpinv_mem->d_ml)
#define mu             (kinpinv_mem->d_mu)
#define smu            (kinpinv_mem->d_smu)
#define jacDQ          (kinpinv_mem->d_jacDQ)
#define djac           (kinpinv_mem->d_djac)
#define J              (kinpinv_mem->d_J)
#define pivots         (kinpinv_mem->d_pivots)
#define nje            (kinpinv_mem->d_nje)
#define nfeDQ          (kinpinv_mem->d_nfeDQ)
#define J_data         (kinpinv_mem->d_J_data)
#define last_flag      (kinpinv_mem->d_last_flag)
#define JTJ            (kinpinv_mem->d_JTJ)
#define regularized    (kinpinv_mem->d_regularized)
#define redojac        (kinpinv_mem->d_redojac)
#define beta           (kinpinv_mem->d_beta)
#define reg_param      (kinpinv_mem->d_reg_param)

/* 
 * =================================================================
 * EXPORTED FUNCTIONS
 * =================================================================
 */
             
/*
 * -----------------------------------------------------------------
 * KINPinv
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the dense pseude-inverse linear solver module. 
 * KINPinv sets the kin_linit, kin_lsetup, kin_lsolve, kin_lfree fields 
 * in *kinmem to be kinPinvInit, kinPinvSetup, kinPinvSolve, and 
 * kinPinvFree, respectively.  
 * It allocates memory for a structure of type KINDlsMemRec and sets 
 * the kin_lmem field in *kinmem to the address of this structure.  
 * It sets setupNonNull in *kinmem to TRUE, and the djac field to the 
 * default kinDlsDenseDQJac.
 * Finally, it allocates memory for J and pivots.
 *
 * NOTE: The dense pseudo-inverse linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, KINPinv will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that N_VGetArrayPointer and
 *       N_VSetArrayPointer exist.
 * -----------------------------------------------------------------
 */

int KINPinv(void *kinmem, int N)
{
  KINMem kin_mem;
  KINPinvMem kinpinv_mem;
  /* Check if kinmem is different from NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINPINV_MEM_NULL, "KINPINV", "KINPinv", MSGD_KINMEM_NULL);
    return(KINPINV_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  /* Test if the NVECTOR package is present */
  if (vec_tmpl->ops->nvgetarraypointer == NULL ||
      vec_tmpl->ops->nvsetarraypointer == NULL) {
    KINProcessError(kin_mem, KINPINV_ILL_INPUT, "KINPINV", "KINPinv", MSGD_BAD_NVECTOR);
    return(KINPINV_ILL_INPUT);
  }

  if (lfree !=NULL) lfree(kin_mem);

  /* Set four main function fields in kin_mem */
  linit  = kinPinvInit;
  lsetup = kinPinvSetup;
  lsolve = kinPinvSolve;
  lfree  = kinPinvFree;

  /* Get memory for KINDlsMemRec */
  kinpinv_mem = NULL;
  kinpinv_mem = (KINPinvMem) malloc(sizeof(struct KINPinvMemRec));
  if (kinpinv_mem == NULL) {
    KINProcessError(kin_mem, KINPINV_MEM_FAIL, "KINPINV", "KINPinv", MSGD_MEM_FAIL);
    return(KINPINV_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = SUNDIALS_DENSE;  

  /* Set default Jacobian routine and Jacobian data */
  jacDQ  = TRUE;
  djac   = NULL;
  J_data = NULL;
  last_flag = KINPINV_SUCCESS;

  setupNonNull = TRUE;

  /* Set problem dimension */
  n = N;

  /* Allocate memory for J,RTR and pivots */
  
  J = NULL;
  J = NewDenseMat(N, N);
  if (J == NULL) {
    KINProcessError(kin_mem, KINPINV_MEM_FAIL, "KINPINV", "KINPinv", MSGD_MEM_FAIL);
    free(kinpinv_mem); kinpinv_mem = NULL;
    return(KINPINV_MEM_FAIL);
  }

  JTJ = NULL;
  JTJ = NewDenseMat(n, n);
  if (JTJ == NULL) {
    KINProcessError(kin_mem, KINPINV_MEM_FAIL, "KINPINV", "KINPinv", MSGD_MEM_FAIL);
    DestroyMat(J);
    free(kinpinv_mem); kinpinv_mem = NULL;
    return(KINPINV_MEM_FAIL);
  }

  pivots = NULL;
  pivots = NewIntArray(N);
  if (pivots == NULL) {
    KINProcessError(kin_mem, KINPINV_MEM_FAIL, "KINPINV", "KINPinv", MSGD_MEM_FAIL);
    DestroyMat(J);
    DestroyMat(JTJ);
    free(kinpinv_mem); kinpinv_mem = NULL;
    return(KINPINV_MEM_FAIL);
  }

  beta = NULL;
  beta = NewRealArray(N);
  if (pivots == NULL) {
    KINProcessError(kin_mem, KINPINV_MEM_FAIL, "KINPINV", "KINPinv", MSGD_MEM_FAIL);
    DestroyMat(J);
    DestroyMat(JTJ);
    DestroyArray(pivots);
    free(kinpinv_mem); kinpinv_mem = NULL;
    return(KINPINV_MEM_FAIL);
  }

  /* This is a direct linear solver */
  inexact_ls = FALSE;

  /* Attach linear solver memory to integrator memory */
  lmem = kinpinv_mem;
  /* Set reg_param 'not set' */
  reg_param = 0;
  redojac=0 ;
  regularized = FALSE ;
  nje   = 0;
  nfeDQ = 0;
  return(KINPINV_SUCCESS);
}

/* 
 * =================================================================
 *  PRIVATE FUNCTIONS
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * kinPinvInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the dense
 * linear solver.
 * -----------------------------------------------------------------
 */

static int kinPinvInit(KINMem kin_mem)
{
  KINPinvMem kinpinv_mem;

  kinpinv_mem = (KINPinvMem) lmem;
  
  nje   = 0;
  nfeDQ = 0;
  
  /*
   * Set where to find jacobian data
   * if jacDQ=True it will be calculated by finite differenses
   * otherwise a user-supplied jacobian will be used
   */

  if (jacDQ) {
    djac = kinPinvDQJac;
    J_data = kin_mem;
  } else {
    J_data = kin_mem->kin_user_data;
  }

  /* Set regularization parameter */
  if (reg_param == 0){
    reg_param = 1;
  }

  last_flag = KINPINV_SUCCESS;
  return(0);
}



/* 
 * Function that calculates the regularized matrix for a given h 
 * The result is returned in res while h is the regularization parameter
 * and J is the current prblem jacobian.
 */

void regMatrix(realtype **JTJ_c, realtype **jac, realtype h, int size)
{
  int i,j,k;


  i = 0;
  j = 0;
  k = 0;
  for (i=0;i<size;i++) {
    for (j=0;j<size;j++){
      
      /*Calculate value at RTR(i,j) */
      JTJ_c[j][i] = 0;
      for (k=0;k<size;k++) JTJ_c[j][i] += jac[j][k]*jac[i][k];
      
      /* add the regularization parameter on the diagonal */
      if (i==j)JTJ_c[j][i] += h*h;
    }
  }
}

/*
 * -----------------------------------------------------------------
 * kinPinvSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the linear solver and
 * prepares J transpose J  + h^2 I if necessary for regularization.
 * -----------------------------------------------------------------
 */

static int kinPinvSetup(KINMem kin_mem)
{
  KINPinvMem kinpinv_mem;
  long int ier;
  int retval;
  int i,j;

  realtype **JTJ_c ;
  realtype **jac ;
  realtype rp ;
  double data;
  realtype *b;

  kinpinv_mem = (KINPinvMem) lmem;

  /* Calculate value of jacobian */
  nje++;
  SetToZero(J);
  retval = djac(n, uu, fval, J, J_data, vtemp1, vtemp2);
  if (retval != 0) {
    last_flag = -1;
    return(-1);
  }

  /* Try to do a LU factorization of J */
  ier = DenseGETRF(J, pivots);

  /* If the LU factorization failed, perform regularization */
  regularized = FALSE;
  if (ier > 0) {
    /* Calculate value of jacobian */
    SetToZero(J);

    /* SetToZero(pivots);*/
    for (i=0; i<n; i++) pivots[i] = 0;
    nje++;
    retval = djac(n, uu, fval, J, J_data, vtemp1, vtemp2);
    if (retval != 0) {
      last_flag = -1;
      return(-1);
    }

    /* Calculate J tranpose J */
    SetToZero(JTJ);

    /* Calculate the regularization parameter */

    jac = J->cols;
    b = N_VGetArrayPointer(fval);
    rp = 0;
    for (i = 0;i<n;i++) {
      data = 0;
      for(j = 0;j<n;j++) data += jac[i][j]*b[j];
      rp += data*data ;
    }
    /*printf("rp: %f \n",rp);*/
    if (sqrt(rp)<1) {
      reg_param = sqrt(rp);
    } else {
      reg_param = 1;
    }
    if (printfl > 0) printf("Regparam: %f, \n ",reg_param);
    /* calculate a regularized matrix */
    JTJ_c = JTJ->cols;

    regMatrix(JTJ_c,jac,reg_param,n);
    
    /* LU-factorize the regularized matrix*/
    ier = DenseGETRF(JTJ,pivots);
    regularized = TRUE;
  }

  redojac = FALSE;
  /* Return 0 if the LU was complete; otherwise return -1 */
  last_flag = ier;
  if (ier > 0) return(-1);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * kinPinvSolve
 * -----------------------------------------------------------------
 * This routine handles the solve operation for the dense linear solver
 * by calling the dense backsolve routine.  The returned value is 0.
 * -----------------------------------------------------------------
 */

static int kinPinvSolve(KINMem kin_mem, N_Vector x, N_Vector b, realtype *res_norm)
{
  KINPinvMem kinpinv_mem;
  realtype **jac;
  realtype *xd;
  realtype *bx;
  int i,j;

  kinpinv_mem = (KINPinvMem) lmem;

  if (redojac) {
    return(1);
  }

  if (regularized) {
    if (printfl > 0) printf("Regularized problem\n");
    /* Calculate new right hand side b = J transpose * b */
    bx = N_VGetArrayPointer(b);
    xd = N_VGetArrayPointer(x);
    jac = J->cols;
    i=0;
    j=0;
    for (i=0;i<n;i++){
      xd[i] = 0;
      for (j=0;j<n;j++) xd[i] += jac[i][j]*bx[j];
    }
    /* Back-solve and get solution in x */

    /*N_VScale(ONE,x,b);*/
    

    DenseGETRS(JTJ, pivots, xd);

    regularized = FALSE;
    /* Calculate a fresh jacobian */
    redojac = TRUE;

    
  } else {
    
    /* Copy the right-hand side into x */
    
    N_VScale(ONE, b, x);
    xd = N_VGetArrayPointer(x);
    
    /* Back-solve and get solution in x */
    DenseGETRS(J, pivots, xd);

  }

  /* Compute the terms Jpnorm and sfdotJp for use in the global strategy
     routines and in KINForcingTerm. Both of these terms are subsequently
     corrected if the step is reduced by constraints or the line search.

     sJpnorm is the norm of the scaled product (scaled by fscale) of
     the current Jacobian matrix J and the step vector p.

     sfdotJp is the dot product of the scaled f vector and the scaled
     vector J*p, where the scaling uses fscale. */


  sJpnorm = N_VWL2Norm(b,fscale);
  N_VProd(b, fscale, b);
  N_VProd(b, fscale, b);
  sfdotJp = N_VDotProd(fval, b);
  
  
  last_flag = KINPINV_SUCCESS;
  
  return(0);
}

/*
 * -----------------------------------------------------------------
 * kinPinvFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the dense linear solver.
 * -----------------------------------------------------------------
 */

static void kinPinvFree(KINMem kin_mem)
{
  KINPinvMem  kinpinv_mem;
  kinpinv_mem = (KINPinvMem) lmem;
  
  DestroyMat(J);
  DestroyMat(JTJ);
  DestroyArray(pivots);
  DestroyArray(beta);
  free(kinpinv_mem); kinpinv_mem = NULL;
}


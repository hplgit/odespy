/*
 * Copyright (C) 2010 Modelon AB
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
 * Johan Ylikiiskilä - johan.ylikiiskila@gmail.com
 *
 */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _REG_ROUTINES_H
#define _REG_ROUTINES_H


#include "slu_ddefs.h"

  /*
   * Routines used in regularisation
   */

SuperMatrix* regSparseMatrix( SuperMatrix *jac, double h);

SuperMatrix* getRegRHS( SuperMatrix *jac, SuperMatrix *B);

double getRegParam(SuperMatrix *jac, SuperMatrix *B);
#endif

#ifdef __cplusplus
}
#endif

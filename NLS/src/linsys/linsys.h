/* --------------------------------------------------------------------- **
** Tecgraf / PUC-Rio ( Computer Graphics Technology Group )              ** 
** --------------------------------------------------------------------- ** 
**                                                                       **  
**           Program for solving linear systems of equations             **  
**                                                                       **  
** File: linsys.h   ( Prototypes and Definitions )                    **  
**                                                                       **
** --------------------------------------------------------------------- **
**                                                                       **
** Developed  by:  Anderson Pereira                                      **
**                 Ivan Menezes                                          **
**                 Shun Wang                                             **
**                 Eric de Sturler                                       **
**                 Glaucio H. Paulino                                    **
**                                                                       **
** Version: 27-Jul-09                                                    **
**                                                                       **
** --------------------------------------------------------------------- */

#ifndef _LINSYS_H_
#define _LINSYS_H_


#ifdef __cplusplus
extern "C" {
#endif

/*
 * A sample matrix to explain the storage schemes
 *
 *  A = [2 1 0 0]
 *      |1 2 1 0|
 *      |0 1 2 1|
 *      [0 0 1 1]
 *  
 *  Size of the matrix
 *   nrow = ncol = 4
 *
 * Profile format
 *
 *  Start index of each column
 *   profile = { 0 | 0 | 1 | 2 }
 *  
 *  Nonzero matrix value
 *   val = { 2 | 1 2 | 1 2 | 1 }
 * 
 * Compressed-Sparse Row format (CSR)
 *
 *  Number of nonzeros in the matrix
 *   nnz = 10
 *
 *  Start index of each row
 *   row_ptr[nrow+1] = { 0 | 2 | 5 | 8 | 10 }
 *  
 *  Column index
 *   col_ind[nnz] = { 0 1 | 1 0 2 | 2 1 3 | 3 2 }
 *  
 *  Nonzero matrix value
 *   val[nnz]     = { 2 1 | 2 1 1 | 2 1 1 | 1 1 }
 */

#include "sparsetools.h"

int Crout_Profile(int flag, double **val, double *b, int n, int *profile);
int PCG_Profile(double **val, int *profile,
  int n, double *b, double *x, int *iter, int maxiter, double tol);
int PCG_CSR(int *row_ptr, int *col_ind, double *val,
  int n, double *b, double *x, int *iter, int maxiter, double tol);

#ifdef __cplusplus
}
#endif

#endif
